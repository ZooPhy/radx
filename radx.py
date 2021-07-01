import argparse
import logging
import os
from multiprocessing import Queue

from radx import PATH_TO_SRA, SRAProcess

logging.basicConfig(format='%(asctime)s: %(levelname)s:%(message)s',
                    filemode='w', filename='logs/radx.log',
                    level=logging.DEBUG)

def main():
    '''Main method : parse input arguments and train'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--multiproc', action='store_false',
                        help='Use multiprocessing to run jobs in parallel')
    parser.add_argument('--maxproc', type=int, default=4,
                        help='Max processes to run in parallel')
    args = parser.parse_args()
    
    # Get all valid fastq.gz files
    sra_files = [x for x in os.listdir(PATH_TO_SRA) if x[-8:]=="fastq.gz"]
    sra_files = sorted(list(set([x.split(".")[0][:-3] for x in sra_files])))
    logging.info("Found %s fastq.gz files.", len(sra_files))
    # If multiprocessing is enabled launch processes in a batch
    if args.multiproc:
        for i in range(0, len(sra_files), args.maxproc):
            batch = sra_files[i:i+args.maxproc]
            logging.info("Starting batch %s", batch)    
            processes = []
            q = Queue()
            for sra_file in batch[:]:
                logging.info("Processing %s", sra_file)
                sra_proc = SRAProcess(sra_file, queue=q)
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
