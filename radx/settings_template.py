# ******************** MAIN SETTINGS FILE ********************
# This holds global variables used locally for system specific settings
# Copy settings_template.py into a new file and rename the new file to to settings.py  
# Do not commit settings.py to GitHub repository as it may contain passwords

### Important paths - REQUIRED
# path where reference sequences (fasta, primers etc.) are held
PATH_TO_REFS = "resources/refs/"

## Paths for secondary tasks - NOT REQUIRED
### We use the Selenium library to automatically download GISAID FASTA and metadata files
SELENIUM_BROWSER = "chrome"

### GISAID username and password
GISAID_USERNAME = ""
GISAID_PASSWORD = ""

# path where mutation info needs to be stored - 
PATH_TO_MUTATIONS = "resources/mutations/"
# path where alcov is available
PATH_TO_ALCOV = "../alcov/"
# path where GISAID sequences are stored
PATH_TO_DOWNLOADS = ""
PATH_TO_GISAID = "resources/gisaid/"
# path to hosting directory for sharing
PATH_TO_HOSTING = ""
# path to aggregate directory for collection (absolute path)
PATH_TO_AGGREGATE = "/Users/amagge/pyspace/radx/resources/aggregate_b05/"
