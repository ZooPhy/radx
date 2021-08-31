# Radx-rad
RADx-rad pipeline for metagenomic data and analysis of SARS-CoV-2 from wastewater samples.

# Please note: This project is currently under development and not functional.

# Requirements

* [Python 3.7+](https://www.python.org)
* [BWA](https://github.com/lh3/bwa)
* [htslib](http://www.htslib.org/download/)
* [samtools](http://www.htslib.org/download/)
* [iVar](https://github.com/andersen-lab/ivar)
* [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

If you installed BWA or other tools at a custom location, you may have to add the executable to the environment PATH table. 
```
export PATH="$PATH:PATH-TO-CUSTOM-LOCATION-CONTAINING-BINARIES"
```
e.g.
```
export PATH="$PATH:/Users/amagge/pyspace/radx/resources/bwa"
export PATH="$PATH:/Users/amagge/pyspace/radx/resources/htslib/bin"
export PATH="$PATH:/Users/amagge/pyspace/radx/resources/samtools/bin"
```

# Installation instructions
Copy ```radx/settings_template.py``` into a new file and rename the new file to ```radx/settings.py``` and change paths if necessary.

For now, the project can be run in standalone mode using the following command:
```
python radx.py
```


