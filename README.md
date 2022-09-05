# RADx-rad v0.2
RADx-rad pipeline for metagenomic data and analysis of SARS-CoV-2 from wastewater samples.

## Please note: This project is currently under development.

Currently, the architecture and information flow in this pipeline is designed as shown below.

![arch](resources/radx-arch.png)

# Requirements
Please install the following 
* [Python 3.7+](https://www.python.org) - Our code runs on Python
* [BWA 0.7.17](https://github.com/lh3/bwa) - Can be installed using bioconda channel
* [htslib 1.14](http://www.htslib.org/download/) - Can be installed using bioconda channel
* [samtools 1.14](http://www.htslib.org/download/) - Can be installed using bioconda channel
* [iVar 1.3](https://github.com/andersen-lab/ivar) - Mac users, please install iVar from source and not from conda
* [bedtools 2.30.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)
* [Freyja 1.3.1](https://github.com/andersen-lab/Freyja) - Requires installation of [Usher](https://usher-wiki.readthedocs.io/en/latest/Installation.html) from source for regular updates of global trees.
* [lofreq 2.1.5](https://github.com/CSB5/lofreq) - Mac users, please install from bioconda
* [pyvcf 0.6.8](https://pyvcf.readthedocs.io/en/latest/)

If you installed BWA or other tools at a custom location, you may have to add the executable to the environment PATH table. You may also have to run them when you restart your computer.  
```
export PATH="$PATH:PATH-TO-CUSTOM-LOCATION-CONTAINING-BINARIES"
```
e.g.
```
export PATH="$PATH:/home/arjun/pyspace/radx/resources/bwa"
export PATH="$PATH:/home/arjun/pyspace/radx/resources/htslib/bin"
export PATH="$PATH:/home/arjun/pyspace/radx/resources/samtools/bin"
```

# Installation instructions
Copy ```radx/settings_template.py``` into a new file and rename the new file to ```radx/settings.py``` and change paths if necessary.

Create a directory ```refs``` under ```resources``` and download/move the following files into them.
* NC_045512.2.fa
* NC_045512.2.gff
* swift_primers.bed
* swift_primers.fasta
* swift_primers.tsv

For now, the project can be run in standalone mode using the following command:
```
python radx.py INPUT_DIRECTORY OUTPUT_DIRECTORY
```

INPUT_DIRECTORY contains all sample fastq.gz files in the format below:
* sample01_R1.fastq.gz
* sample01_R2.fastq.gz
* sample02_R1.fastq.gz
* sample02_R2.fastq.gz
* sample03_R1.fastq.gz
* sample03_R2.fastq.gz


On successful processing, OUTPUT_DIRECTORY will contain subdirectories for each sample with output files in them. 
* sample01
* sample02
* sample03

There are some options when it comes to running radx, you may view these by running:
```
python radx.py --help
```
