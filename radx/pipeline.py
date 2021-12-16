"""File containing the main pipeline
"""
import logging
import os
import subprocess
import time
import pandas as pd
from multiprocessing import Process, Queue
from os import chdir, listdir, makedirs, remove, removedirs, system
from os.path import exists, isdir, isfile, join

from radx.settings import PATH_TO_HOSTING, PATH_TO_REFS
from radx.utils import (read_ivar, read_lofreq, merge_calls,
                        filter_merged_calls)

class SRAProcess(Process):
    def __init__(self, fname, in_dir, out_dir, overwrite=False, queue=None):
        super(SRAProcess, self).__init__()
        self.name = fname
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.overwrite = overwrite

    def run(self):
        logging.info("Running job for %s", self.name)
        self.sdir = join(self.out_dir, self.name)
        # The trimmed sorted file is one of the final files used
        # so skip running the process if it exists already
        final_bam_file = join(self.sdir, self.name+".trimmed.sorted.bam")
        if not self.overwrite and exists(final_bam_file):
            logging.info("Found existing file %s. Skipping processing.", final_bam_file)
            logging.info("To process sample again, please use the overwrite flag.")
        # Process the input files
        self.prep()
        self.align()
        self.trim_to_bam()
        self.variants()
        # Analyze files
        self.collect_metrics()
        self.plot()
        # Perform logistics and cleanup
        self.move_files()
        self.cleanup()

    def prep(self):
        logging.info("Prepping for %s", self.name)
        if not isdir(self.sdir):
            makedirs(self.sdir)
        self.log_path = self.name+"_radx.log"
        self.std_out = self.name+"_stdout_radx.log"
        # Source files
        self.sra_r1 = self.name+"_R1.fastq.gz"
        self.sra_r2 = self.name+"_R2.fastq.gz"
        self.ref_fa = "NC_045512.2.fa"
        self.ref_gff = "NC_045512.2.gff"
        self.primer_bed = "swift_primers.bed"
        # some of the important output files to indicate pipeline processing
        self.sort_bam = self.name+".sorted.bam"
        self.trim_sort_bam = self.name+".trimmed.sorted.bam"
        # move files if not fully processed or they don't exist already
        if not exists(join(self.sdir, self.trim_sort_bam)) or self.overwrite:
            for sfile in [self.sra_r1, self.sra_r2]:
                if not isfile(join(self.sdir, sfile)):
                    self.run_cmd(["cp", join(self.in_dir, sfile), self.sdir])
        # TODO: verify if the following are needed for each record
        # or if they can be just done once and referred across SRA jobs
        for sfile in [self.ref_fa, self.ref_gff, self.primer_bed]:
            if not isfile(join(self.sdir, sfile)):
                self.run_cmd(["cp", join(PATH_TO_REFS, sfile), self.sdir])
        # Update paths to job specific directory
        if isdir("radx"):
            chdir(self.sdir)
            logging.info("Changed directory to work dir: %s", self.sdir)
            logging.info("Logging will be continued in: %s", join(self.sdir, self.log_path))
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(format='%(asctime)s: Sample '+self.name+' %(levelname)s:%(message)s',
                            filemode=('w' if self.overwrite else 'a'), filename=self.log_path,
                            level=logging.DEBUG)
        logging.info("Current dir contents: %s", " ".join(listdir()))
        # Verify if all source files exist
        if self.overwrite:
            assert all([isfile(x) for x in [self.sra_r1, self.sra_r2]])
            assert all([isfile(x) for x in [self.ref_fa, self.ref_gff]])
            assert all([isfile(x) for x in [self.primer_bed]])
        # Variables intermediate and result files
        self.sort_bai = self.name+".sorted.bai"
        self.trim = self.name+".trimmed"
        self.trim_bam = self.name+".trimmed.bam"
        self.final_ivar = self.name + "_ivar.tsv"
        self.sort_dep = self.name+".sorted.depth"
        self.trim_sort_dep = self.name+".trimmed.sorted.depth"
        self.trim_sort_bai = self.name+".trimmed.sorted.bai"
        
        self.trim_sort_indelqual = self.name+"_trim_sort_indelqual"
        self.final_lofreq = self.name+"_lofreq.vcf"
        self.variants_merged = self.name+"_variants_merged.tsv"
        self.stats = self.name+".stats"
        self.plot_dir = "plots"
        self.alcov_dir = "alcov"
        self.metrics = self.name+"_metrics.tsv"

    def align(self):
        if not exists(self.sort_bam) or self.overwrite:
            # first align to reference sequence
            self.run_cmd(["bwa", "index", self.ref_fa])
            self.run_cmd(["bwa", "mem", "-t", "32", self.ref_fa, self.sra_r1, self.sra_r2,
                        "|", "samtools", "view", "-b", "-F", "4",
                        "|", "samtools","sort", "-o", self.sort_bam])

    def trim_to_bam(self):
        if not exists(self.trim_sort_bam) or self.overwrite:
            #trimming primers and base quality
            self.run_cmd(["ivar", "trim", "-b", self.primer_bed, "-p", self.trim, "-i", self.sort_bam])
            #checking trimmed vs non trimmed
            self.run_cmd(["samtools", "sort", "-o", self.trim_sort_bam, self.trim_bam])
            # index the sortedbam file
            self.run_cmd(["samtools", "index", self.trim_sort_bam])
            # get depth of the trimmed and sorted bam file for later
            if not exists(self.trim_sort_dep) or self.overwrite:
                self.run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, ">", self.trim_sort_dep])

    def variants(self):
        # Run ivar to get variants
        if not exists(self.final_ivar) or self.overwrite:
            self.run_cmd(["samtools", "mpileup", "-aa", "-A", "-B", "-d", "0", "--reference " +self.ref_fa+ " -Q", "0", self.trim_sort_bam,
                          "|", "ivar", "variants", "-p", self.final_ivar, "-t", "0", "-q", "20", "-m", "10", "-r", self.ref_fa, "-g", self.ref_gff])
        # Run lofreq to get variants
        if not exists(self.final_lofreq) or self.overwrite:
            self.run_cmd(["lofreq", "indelqual", "--dindel", "-f", self.ref_fa, self.trim_sort_bam, "-o", self.trim_sort_indelqual])
            self.run_cmd(["lofreq", "call", "-f", self.ref_fa, "--call-indels", "-o", self.final_lofreq, self.trim_sort_indelqual])
        
        # Read ivar and lofreq
        ivar_calls = read_ivar(self.final_ivar) if exists(self.final_ivar) else pd.DataFrame()
        if len(ivar_calls.index) == 0:
            logging.warning("Empty ivar dataframe :%s", self.final_ivar)
        lofreq_calls = read_lofreq(self.final_lofreq) if exists(self.final_lofreq) else pd.DataFrame()
        if len(lofreq_calls.index) == 0:
            logging.warning("Empty lofreq dataframe :%s", self.final_lofreq)
        merged_calls = merge_calls(ivar_calls, lofreq_calls)
        filtered_merged_calls = filter_merged_calls(merged_calls)
        filtered_merged_calls.to_csv(self.variants_merged, sep="\t", index=False)
        with open("variants.csv", "w") as ofile:
            if "Variant" in filtered_merged_calls.columns.tolist():
                print(",".join(filtered_merged_calls["Variant"].tolist()), end="", file=ofile)
            else:
                print("", file=ofile)

    def collect_metrics(self):
        if not exists(self.trim_sort_bam):
            logging.warning("Cannot find the bam file :%s", self.trim_sort_bam)
        else:
            # "Sample", "Breadth of coverage", "Total read count", "Mean reads"
            breadth, count, mean, variants = "final.bam.breadth", "final.bam.count", "final.bam.mean", "variants.csv"
            self.run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, "|", 
                          "awk '{c++; if($3>10) total+=1} END {print total*100/c}'",
                          ">", breadth], redirect=False)
            self.run_cmd(["samtools", "view", "-c", "-F", "4", self.trim_sort_bam, ">", count], redirect=False)
            self.run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, "|", 
                          "awk '{c++;s+=$3}END{print s/c}'", ">", mean], redirect=False)
            breadth = [y for y in open(breadth)]
            breadth = breadth[0].strip() if breadth else "0" 
            count = [y for y in open(count)]
            count = count[0].strip() if count else "0" 
            mean = [y for y in open(mean)]
            mean = mean[0].strip() if mean else "0" 
            variants = [y for y in open(variants)]
            variants = variants[0].strip() if variants else "" 
            with open(self.metrics, "w") as ofile:
                ofile.write("\t".join([self.name, breadth, count, mean, variants]))

    def move_files(self):
        if not PATH_TO_HOSTING.strip():
            logging.info("Skipping move to hosting directory. PATH_TO_HOSTING not configured.")
        else:
            path_to_move = join(PATH_TO_HOSTING, self.name)
            if not exists(path_to_move):
                makedirs(path_to_move)
            files_to_move = [self.alcov_dir, self.plot_dir]
            for mfile in files_to_move:
                self.run_cmd(["mv", mfile, path_to_move])

    def plot(self):
        if not exists(self.trim_sort_bam):
            logging.warning("Cannot find the bam file :%s", self.trim_sort_bam)
        else:
            # use samtools to get stats and then plot
            self.run_cmd(["samtools", "stats", self.trim_sort_bam, ">", self.stats])
            self.run_cmd(["plot-bamstats", "-p", self.plot_dir+"/", self.stats])
            # run alcov 
            if not exists(self.alcov_dir):
                makedirs(self.alcov_dir)
            self.run_cmd(["alcov", "find_lineages", self.trim_sort_bam])
            self.run_cmd(["alcov", "find_mutants", self.trim_sort_bam])
            self.run_cmd(["alcov", "amplicon_coverage", self.trim_sort_bam])
            self.run_cmd(["alcov", "gc_depth", self.trim_sort_bam])

    def cleanup(self):
        logging.info("Running cleanup for %s", self.name)
        if isdir("radx"):
            chdir(self.sdir)
        filelist = listdir()
        files_to_keep = [
                            self.sort_bam,
                            self.sort_dep,
                            self.trim_sort_bam,
                            self.trim_sort_bai,
                            self.trim_sort_dep,
                            self.plot_dir,
                            self.alcov_dir,
                            self.stats,
                            self.log_path,
                            self.std_out,
                            self.metrics,                   # metrics file
                            self.final_ivar,                # variants from ivar
                            self.final_lofreq,              # variants from lofreq 
                            self.variants_merged            # merged variants
                        ]
        logging.info("Current dir contents: %s", " ".join(filelist))
        files_to_delete = [x for x in filelist if x not in files_to_keep]
        logging.info("Files to delete: %s", " ".join(files_to_delete))
        for dfile in files_to_delete:
            if isfile(dfile):
                remove(dfile)
            else:
                logging.info("Cannot remove directory %s", dfile)
                # can be dangerous to remove directory, so disabled for now
                # removedirs(dfile)

    def run_cmd(self, cmd_list, redirect=True, timeout=None):
        ret = None
        start_file_list = listdir()
        try:
            # TODO: Subprocess doesn't support the pipe command by default
            if "|" not in cmd_list and ">" not in cmd_list:
                logging.info("--- Running SP '%s'", " ".join(cmd_list))
                if redirect and not isdir("radx"): #dont print unless in the working directory
                    with open(self.std_out, "a") as ofile:
                        print("\n".join(["","-"*64," ".join(cmd_list),"-"*64,""]), file=ofile)
                        ret = subprocess.run(cmd_list, check=True, stdout=ofile, stderr=ofile)
                else:
                    ret = subprocess.run(cmd_list, check=True)
                if ret.returncode != 0:
                    logging.info("--- Return code '%s'", ret.returncode)
            else:
                if redirect and ">" not in cmd_list:
                    cmd_list += [">>", self.std_out]
                logging.info("--- Running SYS '%s'", " ".join(cmd_list))
                if redirect and not isdir("radx"): #dont print unless in the working directory
                    with open(self.std_out, "a") as ofile:
                        print("\n".join(["","-"*64," ".join(cmd_list),"-"*64,""]), file=ofile)
                ret = system(" ".join(cmd_list))
                if ret != 0:
                    logging.info("--- Return code '%s'", ret)
            added_file_list = [x for x in listdir() if x not in start_file_list]
            if not isdir("radx"):
                logging.info("--- Added '%s' files '%s'", len(added_file_list), ", ".join(added_file_list))
        except Exception as e:
            logging.info("--- ERROR encountered when processing command '%s'    Message: '%s'", 
                         ", ".join(cmd_list), repr(e))
        return ret
