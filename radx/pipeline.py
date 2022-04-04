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

from radx.settings import PATH_TO_HOSTING, PATH_TO_REFS, PATH_TO_AGGREGATE
from radx.utils import (read_lofreq, completemask)

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
        self.primer_info = "swif_primers_info.tsv"
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
        self.sort_dep = self.name+".sorted.depth"
        self.trim_sort_dep = self.name+".trimmed.sorted.depth"
        self.trim_sort_bai = self.name+".trimmed.sorted.bam.bai"
        self.out_r1 = self.name+"_R1.fastq"
        self.out_r2 = self.name+"_R2.fastq"
        
        self.trim_sort_indelqual = self.name+"_trim_sort_indelqual"
        self.trim_sort_indelqual_readsremoved = self.name+"_trim_sort_indelqual_readsremoved"
        self.lofreq_1 = self.name+"_lofreq_1.vcf"
        self.lofreq_1_filtered = self.name+"_lofreq_1_filtered.vcf"
        self.lofreq_2 = self.name+"_lofreq_2.vcf"
        self.intersect = self.name+"_intersect.vcf"
        self.union = self.name+"_union.vcf"
        self.union_filtered = self.name+"_union_filtered.vcf"
        self.union_snpEff = self.name+"_union_snpEff.vcf"
        self.union_snpEff_tsv = self.name+"_union_snpEff.tsv"
        self.intersect_renamed = self.name+"_intersect_renamed.vcf"
        self.primer_mismatch_indices = self.name+"_primer_mismatch_indices"
        self.primer_mismatch_indices_v2 = self.name+"_primer_mismatch_indices_v2"
        #self.variants_merged = self.name+"_variants_merged.tsv"
        self.stats = self.name+".stats"
        self.plot_dir = "plots"
        self.freyja_variants = self.name+"_freyja_variants.tsv"
        self.freyja_summary = self.name+"_freyja_summary.tsv"
        self.metrics = self.name+"_metrics.tsv"

    def align(self):
        if not exists(self.sort_bam) or self.overwrite:
            # first align to reference sequence
            self.run_cmd(["bwa", "index", self.ref_fa])
            self.run_cmd(["bwa", "mem", "-t", "32","-v","1", self.ref_fa, self.sra_r1, self.sra_r2,
                        "|", "samtools", "view", "-@", "32","-b", "-f", "1", "-F", "268", "-q", "20", "-s","1.0",
                        "|", "samtools","sort", "-@", "32", "-o", self.sort_bam])

    def trim_to_bam(self):
        if not exists(self.trim_sort_bam) or self.overwrite:
            #trimming primers and base quality
            self.run_cmd(["ivar", "trim", "-e", "-m", "1", "-q", "0", "-b", self.primer_bed, "-p", self.trim, "-i", self.sort_bam])
            #realign with lofreq
            self.run_cmd(["lofreq", "viterbi", "--defqual", "2", "-f", self.ref_fa, self.trim, "|"
            "samtools", "sort", "-@", "32", "-o", self.trim_sort_bam])
            # index the sortedbam file

    def remove_reads(self):
        self.run_cmd(["lofreq", "indelqual", "--dindel", "-f", self.ref_fa, self.trim_sort_bam, "-o", self.trim_sort_indelqual])
        self.run_cmd(["lofreq", "call", "--min-cov", "5", "--max-depth", "1000000", "-q", "30", "-Q", "30", "-e", "--min-mq", "20", "--sig", "0.0005", "--bonf", "dynamic",
         "--no-default-filter", "--no-default-filter", "-f", self.ref_fa, "--call-indels", "-o", self.lofreq_1])
        self.run_cmd(["lofreq", "filter", "-V", "0", "-v", "5", "-a", "0.05", "-A", "0.95", "-i", self.lofreq_1, "-o", self.lofreq_1_filtered])
        self.run_cmd(["ivar", "getmasked", "-i", self.lofreq_1_filtered, "-b", self.primer_bed, "-f", self.primer_info, "-p", self.primer_mismatch_indices])
        complete_mask = completemask(self.primer_mismatch_indices, self.primer_info)
        with open(self.primer_mismatch_indices_v2, 'w') as o:
            o.write(complete_mask + '\n')
        self.run_cmd(["ivar", "removereads", "-i", self.trim_sort_indelqual, "-b", self.primer_bed, "-t", self.primer_mismatch_indices_v2, "-p", self.trim_sort_indelqual_readsremoved]) # final bam file
        if not exists(self.trim_sort_bai) or self.overwrite:
            self.run_cmd(["samtools", "index", self.trim_sort_indelqual_readsremoved])
            # get depth of the trimmed and sorted bam file for later
        if not exists(self.trim_sort_dep) or self.overwrite:
            self.run_cmd(["samtools", "depth", "-a", self.trim_sort_indelqual_readsremoved, ">", self.trim_sort_dep])
            # convert to fastq.gz for uploads
        if not exists(self.out_r1) or not exists(self.out_r2) or self.overwrite:
            self.run_cmd(["bedtools", "bamtofastq", "-i", self.trim_sort_indelqual_readsremoved, "-fq", self.out_r1, "-fq2", self.out_r2])

    def variants(self):
        # Run lofreq to get variants
        self.run_cmd(["lofreq", "call", "--min-cov", "5", "--max-depth", "1000000", "-q", "30", "-Q", "30", "-e", "--min-mq", "20", "--sig", "0.0005", "--bonf", "dynamic",
         "--no-default-filter", "--no-default-filter", "-f", self.ref_fa, "--call-indels", "-o", self.lofreq_2, self.trim_sort_indelqual_readsremoved])
        self.run_cmd(["vcfintersect", "-r", self.ref_fa, "-v", "-w", "0", "-i", self.lofreq_2, self.lofreq_1_filtered, ">", self.intersect])
        self.run_cmd(["sed", "-r", "--sandbox", "-e", 's/^(#CHROM.+)$/##FILTER=<ID=AmpliconRemoval,Description="Variant removed upon removal of amplicon">\n\1/g', "-e", 's/(.+\t)PASS(\t.+)/\1AmpliconRemoval\2/g', self.intersect > self.intersect_renamed])
        self.run_cmd(["vcfintersect", "-r", self.ref_fa, "-w", "0", "-u", self.intersect_renamed, self.lofreq_2, ">", self.union])
        self.run_cmd(["lofreq", "filter", "-V", "0", "-v", "0", "-a", "0.0", "-A", "0.0", "-b", "fdr", "-c", "0.001", "--print-all", "-i", self.union, "-o", self.union_filtered])
        self.run_cmd(["snpEff", "eff", "-nodownload", "-dataDir", "/data", "-i", "vcf", "-o", "vcf", "-formatEff", "-classic", "-no-downstream", "-no-intergenic", "-no-upstream", "-ud", "0", "-stats", "stats.html", "-noLog", "NC_045512.2", self.union_filtered, ">", self.union_snpEff])
        # variants are recorded in metrics, so don't run this if it exists
        if not exists(self.metrics) or self.overwrite:
            # Read vcf lofreq final
            lofreq_calls = read_lofreq(self.union_snpEff) if exists(self.union_snpEff) else pd.DataFrame()
            if len(lofreq_calls.index) == 0:
                logging.warning("Empty lofreq dataframe :%s", self.final_lofreq)
            lofreq_calls.to_csv( self.union_snpEff_tsv, sep="\t", index=False)
            with open("variants.csv", "w") as ofile:
                if "Variant" in lofreq_calls.columns.tolist():
                    print(",".join(lofreq_calls["Variant"].tolist()), end="", file=ofile)
                else:
                    print("", file=ofile)

    def collect_metrics(self):
        if not exists(self.metrics) or self.overwrite:
            if not exists(self.trim_sort_bam):
                logging.warning("Cannot find the bam file :%s . Cannot collect metrics.", self.trim_sort_bam)
            else:
                # "Sample", "Breadth of coverage", "Total read count", "Mean reads"
                breadth, count, mean, variants = "final.bam.breadth", "final.bam.count", "final.bam.mean", "variants.csv"
                self.run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, "|", 
                            "awk '{c++; if($3>10) total+=1} END {print total*100/c}'",
                            ">", breadth], redirect=False)
                self.run_cmd(["samtools", "view", "-c", "-F", "4", self.trim_sort_bam, ">", count], redirect=False)
                self.run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, "|", 
                            "awk '{c++;s+=$3}END{print s/c}'", ">", mean], redirect=False)
                if exists(breadth):
                    breadth = [y for y in open(breadth)]
                    breadth = breadth[0].strip() if breadth else "0" 
                else:
                    breadth = "0" 
                if exists(count):
                    count = [y for y in open(count)]
                    count = count[0].strip() if count else "0" 
                else:
                    count = "0" 
                if exists(mean):
                    mean = [y for y in open(mean)]
                    mean = mean[0].strip() if mean else "0" 
                else:
                    mean = "0" 
                if exists(variants):
                    variants = [y for y in open(variants)]
                    variants = variants[0].strip() if variants else "" 
                else:
                    variants = "0" 
                    with open(self.metrics, "w") as ofile:
                        ofile.write("\t".join([self.name, breadth, count, mean, variants]))

    def plot(self):
        if not exists(self.trim_sort_bam):
            logging.warning("Cannot find the bam file :%s", self.trim_sort_bam)
        else:
            # use samtools to get stats and then plot
            # if not exists(self.stats) or self.overwrite: 
                # self.run_cmd(["samtools", "stats", self.trim_sort_bam, ">", self.stats], redirect=False)
                # self.run_cmd(["plot-bamstats", "-p", self.plot_dir+"/", self.stats], redirect=False)

            # run freyja variant generation - this has to be done again for some reason
            if not exists(self.freyja_variants) or self.overwrite:
                self.run_cmd(["samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f", self.ref_fa, self.trim_sort_bam, "|",
                              "tee >(cut -f1-4 > "+self.trim_sort_dep+")", "|",
                              "ivar variants -p", self.freyja_variants, "-q 20 -t 0.0 -r", self.ref_fa])
            # run freyja summary
            if not exists(self.freyja_summary) or self.overwrite:
                self.run_cmd(["freyja", "demix", self.freyja_variants, self.trim_sort_dep, "--output", self.freyja_summary])

    def move_files(self):
        if not PATH_TO_HOSTING.strip():
            logging.info("Skipping move to hosting directory. PATH_TO_HOSTING not configured.")
        else:
            path_to_move = join(PATH_TO_HOSTING, self.name)
            if not exists(path_to_move):
                makedirs(path_to_move)
            files_to_move = [self.plot_dir]
            for mfile in files_to_move:
                self.run_cmd(["mv", mfile, path_to_move])
        if not PATH_TO_AGGREGATE.strip():
            logging.info("Skipping move to hosting directory. PATH_TO_HOSTING not configured.")
        else:
            files_to_copy = [self.freyja_summary]
            for cfile in files_to_copy:
                self.run_cmd(["cp", cfile, PATH_TO_AGGREGATE])

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
                            self.trim_sort_indelqual_readsremoved,
                            self.out_r1,
                            self.out_r2,
                            self.plot_dir,
                            self.stats,
                            self.log_path,
                            self.std_out,
                            self.metrics,                     # metrics file
                            self.union_snpEff,                  # variants from lofreq vcf
                            self.union_snpEff_tsv,              # variants from lofreq final tsv
                            self.freyja_variants,
                            self.freyja_summary
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
                    logging.info("--- Running SP in shell '%s'", " ".join(cmd_list))
                    if redirect and not isdir("radx"): #dont print unless in the working directory
                        with open(self.std_out, "a") as ofile:
                            print("\n".join(["","-"*64," ".join(cmd_list),"-"*64,""]), file=ofile)
                            ret = subprocess.run(" ".join(cmd_list), shell=True,
                                                        executable="/bin/bash",
                                                        check=True, stdout=ofile,
                                                        stderr=ofile)
                    else:
                        ret = subprocess.run(cmd_list, check=True)
                    if ret.returncode != 0:
                        logging.info("--- Return code '%s'", ret.returncode)
            except Exception as e:
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
