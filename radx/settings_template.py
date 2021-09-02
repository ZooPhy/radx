# ******************** MAIN SETTINGS FILE ********************
# This holds global variables used locally for system specific settings
# Copy settings_template.py into a new file and rename the new file to to settings.py  
# Do not commit settings.py to GitHub repository as it may contain passwords

### We use the Selenium library to automatically download GISAID FASTA and metadata files
SELENIUM_BROWSER = "chrome"

### GISAID username and password
GISAID_USERNAME = ""
GISAID_PASSWORD = ""

### Important paths
# path where reference sequences (fasta, primers etc.) are held
PATH_TO_REFS = "resources/refs/"
# path where GISAID sequences are stored
PATH_TO_DOWNLOADS = ""
PATH_TO_GISAID = "resources/gisaid/"
# path where local fastq.gz files are stored
PATH_TO_SRA = "resources/sra/"
# path local parent working directory
PATH_TO_JOBS = "resources/jobs/"
# path to hosting directory for sharing
PATH_TO_HOSTING = ""