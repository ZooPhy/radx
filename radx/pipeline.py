"""File containing the main pipeline
"""
import logging
import os
import subprocess
import time
from multiprocessing import Process, Queue
from os import chdir, listdir, makedirs, remove, removedirs, system
from os.path import exists, isdir, isfile, join

from radx.settings import (PATH_TO_HOSTING, PATH_TO_JOBS, PATH_TO_REFS,
                           PATH_TO_SRA)

class SRAProcess(Process):
    def __init__(self, fname, overwrite=False, queue=None):
        super(SRAProcess, self).__init__()
        self.name = fname
        self.sdir = join(PATH_TO_JOBS, self.name)
        self.overwrite = overwrite

    def run(self):
        logging.info("Running job for %s", self.name)
        # The masked sorted file is one of the final files used
        # so skip running the process if it exists already
        final_bam_file = join(self.sdir, self.name+".masked.sorted.bam")
        if not self.overwrite and exists(final_bam_file):
            logging.info("Found existing file %s. Skipping processing.", final_bam_file)
            logging.info("To process sample again, please use the overwrite flag.")
            return
        # Process the input files
        self.prep()
        self.align()
        self.trim_sra()
        self.get_depths()
        self.consensus()
        # Analyze files
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
        self.primer_fa = "swift_primers.fasta"
        self.primer_tsv = "swift_primers.tsv"
        # move files if they don't exist already
        for sfile in [self.sra_r1, self.sra_r2]:
            if not isfile(join(self.sdir, sfile)):
                self.run_cmd(["cp", join(PATH_TO_SRA, sfile), self.sdir])
        # TODO: verify if the following are needed for each record
        # or if they can be just done once and referred across SRA jobs
        for sfile in [self.ref_fa, self.ref_gff, self.primer_bed, self.primer_fa, self.primer_tsv]:
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
        assert all([isfile(x) for x in [self.sra_r1, self.sra_r2]])
        assert all([isfile(x) for x in [self.ref_fa, self.ref_gff]])
        assert all([isfile(x) for x in [self.primer_bed, self.primer_fa, self.primer_tsv]])
        # Variables intermediate and result files
        self.sort_bam = self.name+".sorted.bam"
        self.sort_bai = self.name+".sorted.bai"
        self.trim = self.name+".trimmed"
        self.trim_sort_bam = self.name+".trimmed.sorted.bam"
        self.trim_bam = self.name+".trimmed.bam"
        self.final = self.name + "_final"
        self.sort_dep = self.name+".sorted.depth"
        self.trim_sort_dep = self.name+".trimmed.sorted.depth"
        self.trim_cons = self.name+".trimmed.consensus"
        self.trim_cons_fa = self.name+".trimmed.consensus.fa"
        self.trim_cons_tsv = self.name+".trimmed.consensus.tsv"
        self.sw_prim_cons_bam = self.name+"_swift_primers_consensus.bam"
        self.sw_prim_cons_bed = self.name+"_swift_primers_consensus.bed"
        self.prim_mis_ind = "primer_mismatchers_indices"
        self.mask_bam = self.name+".masked.bam"
        self.mask_sort_bam = self.name+".masked.sorted.bam"
        self.mask_sort_bai = self.name+".masked.sorted.bai"
        self.mask_sort_dep = self.name+".masked.sorted.depth"
        self.final_mask_0freq = self.name+"_final_masked_0freq"
        self.stats = self.name+".stats"
        self.plot_dir = "plots"
        self.alcov_dir = "alcov"

    def align(self):
        # index the reference sequences
        self.run_cmd(["bwa", "index", self.ref_fa])
        self.run_cmd(["samtools", "faidx", self.ref_fa])
        # first align to reference sequence
        self.run_cmd(["bwa", "mem", "-t", "32", self.ref_fa, self.sra_r1, self.sra_r2,
                      "|", "samtools", "view", "-b", "-F", "4",
                      "|", "samtools","sort", "-o", self.sort_bam])

    def trim_sra(self):
        #trimming primers and base quality
        self.run_cmd(["bwa", "index", self.sort_bam])
        self.run_cmd(["ivar", "trim", "-b", self.primer_bed, "-p", self.trim, "-i", self.sort_bam])
        #checking trimmed vs non trimmed
        self.run_cmd(["samtools", "sort", "-o", self.trim_sort_bam, self.trim_bam])
        # TODO: Check if there's a better name for the following 
        if not exists(self.final) or self.overwrite:
            self.run_cmd(["samtools", "mpileup" "-aa" "-A" "-B" "-d" "0" "--reference" +self.ref_fa+ "-Q" "0", self.trim_sort_bam,
                          "|", "ivar", "variants", "-p", self.final, "-t", "0.3", "-q", "20", "-m", "10", "-r", self.ref_fa, "-g", self.ref_gff])
        # index the sortedbam file TODO: The following should be above?
        self.run_cmd(["samtools", "index", self.trim_sort_bam])

    def get_depths(self):
        # get depth of the trimmed and sorted sorted bam file for later
        if not exists(self.trim_sort_dep) or self.overwrite:
            self.run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, ">", self.trim_sort_dep])
        if not exists(self.sort_dep) or self.overwrite:
            self.run_cmd(["samtools", "depth", "-a", self.sort_bam, ">", self.sort_dep])

    def consensus(self):
        #getting consensus and info for masking
        if not exists(self.trim_cons) or self.overwrite:
            self.run_cmd(["samtools", "mpileup", "-A", "-d", "0", "-Q", "0", self.trim_sort_bam, 
                          "|", "ivar", "consensus", "-p", self.trim_cons, "-n", "N"])
            self.run_cmd(["bwa", "index", "-p", self.trim_cons, self.trim_cons_fa])

        if not exists(self.sw_prim_cons_bam) or self.overwrite:
            self.run_cmd(["bwa", "mem", "-k", "5", "-T", "16", self.trim_cons, self.primer_fa, 
                          "|", "samtools", "view", "-bS", "-F", "4", 
                          "|", "samtools", "sort", "-o", self.sw_prim_cons_bam])

        self.run_cmd(["samtools", "mpileup", "-A", "-d", "0", "--reference", self.trim_cons_fa, "-Q", "0", self.sw_prim_cons_bam, 
                      "|", "ivar", "variants", "-p", self.trim_cons, "-t", "0.3"])

        if not exists(self.sw_prim_cons_bed) or self.overwrite:
            self.run_cmd(["bedtools", "bamtobed", "-i", self.sw_prim_cons_bam, ">", self.sw_prim_cons_bed])

        if not exists(self.prim_mis_ind) or self.overwrite:
            self.run_cmd(["ivar", "getmasked", "-i", self.trim_cons_tsv, "-b", self.sw_prim_cons_bed, "-f", self.primer_tsv, "-p", self.prim_mis_ind])

        #final analysis and masking
        if not exists(self.mask_bam) or self.overwrite:
            self.run_cmd(["ivar", "removereads", "-i", self.trim_sort_bam, "-p", self.mask_bam, "-t", self.prim_mis_ind+".txt", "-b", self.primer_bed])

        if not exists(self.mask_sort_bam) or self.overwrite:
            self.run_cmd(["samtools", "sort", "-o", self.mask_sort_bam, self.mask_bam])
            # index the masked bam file for alcov analysis
            self.run_cmd(["samtools", "index", self.mask_sort_bam])

        if not exists(self.mask_sort_dep) or self.overwrite:
            self.run_cmd(["samtools", "depth", "-a", self.mask_sort_bam, ">", self.mask_sort_dep])

        if not exists(self.final_mask_0freq) or self.overwrite:
            self.run_cmd(["samtools", "mpileup", "-aa", "-A", "-B", "-d", "0", "--reference", self.ref_fa, "-Q", "0", self.mask_sort_bam, 
                          "|", "ivar", "variants", "-p", self.final_mask_0freq, "-t", "0", "-q", "20", "-m", "20", "-r", self.ref_fa, "-g", self.ref_gff])

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
        if not exists(self.mask_sort_bam):
            logging.warning("Cannot find the bam file :%s", self.mask_sort_bam)
        else:
            # use samtools to get stats and then plot
            # samtools stats B1.masked.sorted.bam > B1.masked.sorted.bam.stats
            self.run_cmd(["samtools", "stats", self.mask_sort_bam, ">", self.stats])
            # plot-bamstats -p B1.masked.sorted.bam.plot B1.masked.sorted.bam.stats
            self.run_cmd(["plot-bamstats", "-p", self.plot_dir+"/", self.stats])
            # run alcov 
            if not exists(self.alcov_dir):
                makedirs(self.alcov_dir)
            self.run_cmd(["alcov", "find_lineages", self.mask_sort_bam])
            self.run_cmd(["alcov", "find_mutants", self.mask_sort_bam])
            self.run_cmd(["alcov", "amplicon_coverage", self.mask_sort_bam])
            self.run_cmd(["alcov", "gc_depth", self.mask_sort_bam])

    def cleanup(self):
        logging.info("Running cleanup for %s", self.name)
        if isdir("radx"):
            chdir(self.sdir)
        filelist = listdir()
        files_to_keep = [
                            # self.sort_bam,
                            # self.trim_sort_bam,
                            self.mask_sort_bam, 
                            self.mask_sort_bai,
                            self.mask_sort_dep,
                            self.final_mask_0freq,
                            self.final+".tsv",
                            self.plot_dir,
                            self.alcov_dir,
                            self.stats,
                            self.log_path,
                            self.std_out
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

    def run_cmd(self, cmd_list, timeout=None):
        ret = None
        # TODO: Subprocess doesn't support the pipe command by default
        if "|" not in cmd_list and ">" not in cmd_list:
            logging.info("--- Running SP '%s'", " ".join(cmd_list))
            if not isdir("radx"): #dont print unless in the working directory
                with open(self.std_out, "a") as ofile:
                    print("\n".join(["","-"*64," ".join(cmd_list),"-"*64,""]), file=ofile)
                    ret = subprocess.run(cmd_list, check=True, stdout=ofile, stderr=ofile)
            else:
                ret = subprocess.run(cmd_list, check=True)
            if ret.returncode != 0:
                logging.info("--- Return code '%s'", ret.returncode)
        else:
            if ">" not in cmd_list:
                cmd_list += [">>", self.std_out]
            logging.info("--- Running SYS '%s'", " ".join(cmd_list))
            if not isdir("radx"): #dont print unless in the working directory
                with open(self.std_out, "a") as ofile:
                    print("\n".join(["","-"*64," ".join(cmd_list),"-"*64,""]), file=ofile)
            ret = system(" ".join(cmd_list))
            if ret != 0:
                logging.info("--- Return code '%s'", ret)
        return ret
