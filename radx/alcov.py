import logging

class BAMAnalyzer(object):
    def __init__(self, name, bam_path):
        self.name = name
        self.path = bam_path

    def run(self):
        logging.info("Running analysis for %s", self.name)
 