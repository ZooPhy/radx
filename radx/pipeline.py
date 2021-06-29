"""File containing the main pipeline
"""
import logging
import subprocess
from os import chdir, listdir, makedirs, system
from os.path import isdir, isfile, join

from radx.settings import PATH_TO_REFS, PATH_TO_SRA
from radx.utils import run_cmd


class SRARecord(object):
    def __init__(self, fname):
        self.name = fname
        self.sdir = join(PATH_TO_SRA, self.name)
        self.sra_r1 = join(PATH_TO_SRA, self.name+"_R1.fastq.gz")
        self.sra_r2 = join(PATH_TO_SRA, self.name+"_R2.fastq.gz")
        self.ref_fa = join(PATH_TO_REFS, "NC_045512.2.fa")
        self.ref_gff = join(PATH_TO_REFS, "NC_045512.2.gff")
        self.primer_bed = join(PATH_TO_REFS, "swift_primers.bed")
        self.primer_fa = join(PATH_TO_REFS, "swift_primers.fasta")
        self.primer_tsv = join(PATH_TO_REFS, "swift_primers.tsv")

    def process(self):
        self.prep()
        self.align()
        self.trim_sra()
        self.consensus()

    def prep(self):
        if not isdir(self.sdir):
            makedirs(self.sdir)
        run_cmd(["cp", self.sra_r1, self.sdir])
        run_cmd(["cp", self.sra_r2, self.sdir])
        # TODO: verify if the following are needed for each record
        # or if they can be just done once and referred
        run_cmd(["cp", self.ref_fa, self.sdir])
        run_cmd(["cp", self.ref_gff, self.sdir])
        run_cmd(["cp", self.primer_bed, self.sdir])
        run_cmd(["cp", self.primer_fa, self.sdir])
        run_cmd(["cp", self.primer_tsv, self.sdir])
        # Update paths
        if isdir("radx"):
            chdir(self.sdir)
            logging.info("Current dir contents: %s", " ".join(listdir()))
        self.ref_fa = "NC_045512.2.fa"
        self.ref_gff = "NC_045512.2.gff"
        self.sra_r1 = self.name+"_R1.fastq.gz"
        self.sra_r2 = self.name+"_R2.fastq.gz"
        self.primer_bed = "swift_primers.bed"
        self.primer_fa = "swift_primers.fasta"
        self.primer_tsv = "swift_primers.tsv"
        # Verify if all files exist
        assert all([isfile(x) for x in [self.sra_r1, self.sra_r2]])
        assert all([isfile(x) for x in [self.ref_fa, self.ref_gff]])
        assert all([isfile(x) for x in [self.primer_bed, self.primer_fa, self.primer_tsv]])

    def align(self):
        # index the reference sequences
        run_cmd(["bwa", "index", self.ref_fa])
        run_cmd(["samtools", "faidx", self.ref_fa])
        self.sort_bam = self.name+".sorted.bam"
        # first align to reference sequence
        run_cmd(["bwa", "mem", "-t", "32", self.ref_fa, self.sra_r1, self.sra_r2,
                 "|", "samtools", "view", "-b", "-F", "4",
                 "|", "samtools","sort", "-o", self.sort_bam])

    def trim_sra(self):
        #trimming primers and base quality
        run_cmd(["bwa", "index", self.sort_bam])
        self.trim = self.name+".trimmed"
        # TODO: Following doesn't seem to work on Mac. Try earlier version.
        run_cmd(["ivar", "trim", "-b", self.primer_bed, "-p", self.trim, "-i", self.sort_bam])
        #checking trimmed vs non trimmed
        self.trim_sort_bam = self.name+".trimmed.sorted.bam"
        self.trim_bam = self.name+".trimmed.bam"
        run_cmd(["samtools", "sort", "-o", self.trim_sort_bam, self.trim_bam])
        # TODO: Check if there's a better name for the following 
        self.final = self.name + "_final"
        run_cmd(["samtools", "mpileup" "-aa" "-A" "-B" "-d" "0" "--reference" "+self.ref_fa+" "-Q" "0", self.trim_sort_bam,
                 "|", "ivar", "variants", "-p", self.final, "-t", "0.3", "-q", "20", "-m", "10", "-r", self.ref_fa, "-g", self.ref_gff])
        # index the sortedbam file
        run_cmd(["samtools", "index", self.trim_sort_bam])
        # get depth of the trimmed and sorted sorted bam file for later
        self.sort_dep = self.name+".sorted.depth"
        self.trim_sort_dep = self.name+".trimmed.sorted.depth"
        run_cmd(["samtools", "depth", "-a", self.trim_sort_bam, ">", self.trim_sort_dep])
        run_cmd(["samtools", "depth", "-a", self.sort_bam, ">", self.sort_dep])

    def consensus(self):
        #getting consensus and info for masking
        self.trim_cons = self.name+".trimmed.consensus"
        self.trim_cons_fa = self.name+".trimmed.consensus.fa"
        self.trim_cons_tsv = self.name+".trimmed.consensus.tsv"
        self.sw_prim_cons_bam = self.name+"_swift_primers_consensus.bam"
        self.sw_prim_cons_bed = self.name+"_swift_primers_consensus.bed"
        self.prim_mis_ind = "primer_mismatchers_indices"
        run_cmd(["samtools", "mpileup", "-A", "-d", "0", "-Q", "0", self.trim_sort_bam, 
                 "|", "ivar", "consensus", "-p", self.trim_cons, "-n", "N"])
        run_cmd(["bwa", "index", "-p", self.trim_cons, self.trim_cons_fa])

        run_cmd(["bwa", "mem", "-k", "5", "-T", "16", self.trim_cons, self.primer_fa, 
                 "|", "samtools", "view", "-bS", "-F", "4", 
                 "|", "samtools", "sort", "-o", self.sw_prim_cons_bam])

        run_cmd(["samtools", "mpileup", "-A", "-d", "0", "--reference", self.trim_cons_fa, "-Q", "0", self.sw_prim_cons_bam, 
                 "|", "ivar", "variants", "-p", self.trim_cons, "-t", "0.3"])

        run_cmd(["bedtools", "bamtobed", "-i", self.sw_prim_cons_bam, ">", self.sw_prim_cons_bed])

        run_cmd(["ivar", "getmasked", "-i", self.trim_cons_tsv, "-b", self.sw_prim_cons_bed, "-f", self.primer_tsv, "-p", self.prim_mis_ind])

        #final analysis and masking
        self.mask_bam = self.name+".masked.bam"
        self.mask_sort_bam = self.name+".masked.sorted.bam"
        self.mask_sort_dep = self.name+".masked.sorted.depth"
        self.final_mask_0freq = self.name+"_final_masked_0freq"
        run_cmd(["ivar", "removereads", "-i", self.trim_sort_bam, "-p", self.mask_bam, "-t", self.prim_mis_ind+".txt", "-b", self.primer_bed])

        run_cmd(["samtools", "sort", "-o", self.mask_sort_bam, self.mask_bam])

        run_cmd(["samtools", "depth", "-a", self.mask_sort_bam, ">", self.mask_sort_dep])

        run_cmd(["samtools", "mpileup", "-aa", "-A", "-B", "-d", "0", "--reference", self.ref_fa, "-Q", "0", self.mask_sort_bam, 
                 "|", "ivar", "variants", "-p", self.final_mask_0freq, "-t", "0", "-q", "20", "-m", "20", "-r", self.ref_fa, "-g", self.ref_gff])
