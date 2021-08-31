"""File containing class declarations
"""

class GISAIDRecord(object):
    def __init__(self, gdict):
        self.virus_name = gdict['Virus name']
        self.virus_type = gdict['Type']
        self.accession_id = gdict['Accession ID']
        self.collection_date = gdict['Collection date']
        self.location = gdict['Location']
        self.addn_loc_info = gdict['Additional location information']
        self.sequence_length = gdict['Sequence length']
        self.host = gdict['Host']
        self.patient_age = gdict['Patient age']
        self.gender = gdict['Gender']
        self.clade = gdict['Clade']
        self.pango_lineage = gdict['Pango lineage']
        self.pango_version = gdict['Pangolin version']
        self.variant = gdict['Variant']
        self.aa_subst = gdict['AA Substitutions']
        self.submission_date = gdict['Submission date']
        self.is_reference = gdict['Is reference?']
        self.is_complete = gdict['Is complete?']
        self.is_high_coverage = gdict['Is high coverage?']
        self.is_low_coverage = gdict['Is low coverage?']
        self.n_content = gdict['N-Content']

    