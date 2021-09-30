import argparse
import logging
import os
import sys
from multiprocessing import Queue

from radx import PATH_TO_SRA, SRAProcess

logging.basicConfig(format='%(asctime)s: %(levelname)s:%(message)s',
                    filemode='a', filename='logs/radx.log',
                    level=logging.DEBUG)

def main():
    '''Main method : parse input arguments and train'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true',
                        help='Redo all steps and overwrite preexisting files')
    parser.add_argument('--include', type=str, default=None,
                        help='Run pipeline on the subset of files (comma separated)')
    parser.add_argument('--multiproc', action='store_false',
                        help='Use multiprocessing to run jobs in parallel')
    parser.add_argument('--maxproc', type=int, default=4,
                        help='Max processes to run in parallel')
    args = parser.parse_args()
    
    # Get all valid fastq.gz files
    sra_files = [x for x in os.listdir(PATH_TO_SRA) if x[-8:]=="fastq.gz"]
    # Store them without the extensions
    sra_files = [x.split(".")[0][:-3] for x in sra_files]
    logging.info("Found %s fastq.gz files from %s samples", len(sra_files), len(set(sra_files)))
    if args.include:
        included_files = [x.strip() for x in args.include.split(",")]
        logging.info("Processing input subsample %s", len(included_files))
        sra_files = [x for x in sra_files if x in included_files]
    sra_files = sorted(list(set(sra_files)))
    logging.info("Processing %s samples", len(sra_files))
    # If multiprocessing is enabled launch processes in a batch
    if args.multiproc:
        for i in range(0, len(sra_files), args.maxproc):
            batch = sra_files[i:i+args.maxproc]
            logging.info("Starting batch %s", batch)
            processes = []
            q = Queue()
            for sra_file in batch[:]:
                logging.info("Processing %s", sra_file)
                # TODO : Probably better to pass overwrite with start
                sra_proc = SRAProcess(sra_file,
                                      overwrite=args.overwrite,
                                      queue=q)
                sra_proc.start()
                processes.append(sra_proc)
            [proc.join() for proc in processes]
            logging.info("Finished batch %s", batch)
        logging.info("Done")
    else:
        for sra_file in sra_files[:]:
            logging.info("Processing %s", sra_file)
            sra_proc = SRAProcess(sra_file)
            sra_proc.start()

if __name__ == '__main__':
    main()
