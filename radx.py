import argparse
import logging
import os

from radx import PATH_TO_SRA, SRARecord

logging.basicConfig(format='%(asctime)s: %(levelname)s:%(message)s',
                    filemode='w', filename='logs/server.log',
                    level=logging.DEBUG)

def main():
    '''Main method : parse input arguments and train'''
    # TODO: Process arguments
    # parser = argparse.ArgumentParser()
    # parser.add_argument('mode', type=str,
    #                     help='Run standalone/server mode')
    # args = parser.parse_args()
    # Get all valid files
    sra_files = [x for x in os.listdir(PATH_TO_SRA) if x[-8:]=="fastq.gz"]
    logging.info("Found %s fastq.gz files.", len(sra_files))
    sra_files = sorted(list(set([x.split(".")[0][:-3] for x in sra_files])))
    for sra_file in sra_files[:1]:
        logging.info("Processing %s", sra_file)
        sra_rec = SRARecord(sra_file)
        sra_rec.process()

if __name__ == '__main__':
    main()
