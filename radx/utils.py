import logging
from os.path import exists
import vcf
import pandas as pd
import numpy as np

GENE_MAP = {'ORF1a': [266, 13468],
            'ORF1b': [13468, 21555],
            'S': [21563, 25384],
            'ORF3a': [25393, 26220],
            'E': [26245, 26472],
            'M': [26523, 27191],
            'ORF6': [27202, 27387],
            'ORF7a': [27394, 27759],
            'ORF7b': [27756, 27887],
            'ORF8': [27894, 28259],
            'N': [28274, 29533],
            'ORF10': [29558, 29674]}

def read_lofreq(filename):
    lofreq_calls = pd.DataFrame(columns=["REGION", "POS", "REF", "ALT", "QUAL", "FILTER", 
                                         "REF_DP", "REF_RV", "ALT_DP", "ALT_RV",
                                         "ALT_FREQ", "TOTAL_DP", "STRAND-BIAS", "EFFECT"])
    vcf_reader = vcf.Reader(filename=filename)
    for row in vcf_reader:
        if len(row.FILTER) == 0 and "EFF" not in row.INFO:
            lofreq_calls = lofreq_calls.append({"REGION": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": str(row.REF),
                                            "ALT": str(row.ALT[0]),
                                            "QUAL": row.QUAL, 
                                            "FILTER": "",
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
                                            "STRAND-BIAS":row.INFO["SB"],
                                            "EFFECT":""
                                           }, 
                                           ignore_index=True)
        if len(row.FILTER) > 0 and "EFF" in row.INFO:
            lofreq_calls = lofreq_calls.append({"REGION": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": str(row.REF),
                                            "ALT": str(row.ALT[0]),
                                            "QUAL": row.QUAL, 
                                            "FILTER": row.FILTER[0],
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
                                            "STRAND-BIAS":row.INFO["SB"],
                                            "EFFECT":row.INFO["EFF"]
                                           }, 
                                           ignore_index=True)
        if len(row.FILTER) == 0 and "EFF" in row.INFO:
            lofreq_calls = lofreq_calls.append({"REGION": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": str(row.REF),
                                            "ALT": str(row.ALT[0]),
                                            "QUAL": row.QUAL, 
                                            "FILTER": "",
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
                                            "STRAND-BIAS":row.INFO["SB"],
                                            "EFFECT":row.INFO["EFF"]
                                           }, 
                                           ignore_index=True)
        if len(row.FILTER) > 0 and "EFF" not in row.INFO:
            lofreq_calls = lofreq_calls.append({"REGION": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": str(row.REF),
                                            "ALT": str(row.ALT[0]),
                                            "QUAL": row.QUAL, 
                                            "FILTER": row.FILTER[0],
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
                                            "STRAND-BIAS":row.INFO["SB"],
                                            "EFFECT":""
                                           }, 
                                           ignore_index=True)
    if lofreq_calls.empty:
        return lofreq_calls
    lofreq_calls["Variant"] = lofreq_calls.apply(lambda row:
                                         str(row["POS"]) + ":" + \
                                         str(row["ALT"][1:]) if len(row["ALT"]) > 1 else \
                                         (str(row["POS"]) + "-" + \
                                         str(row["POS"] + len(row["REF"][1:]) + 1) if len(row["REF"]) > 1 else \
                                         str(row["REF"]) + str(row["POS"]) + \
                                         str(row["ALT"])), axis=1)
    #lofreq_calls = lofreq_calls[lofreq_calls["FILTER"] != "AmpliconRemoval"]
    return lofreq_calls

def completemask(file, primerinfo):
    with open(file) as i:
        getmasked_output = i.readline().strip()
    if not getmasked_output:
        logging.info('No affected primer binding sites found!')
        result = ""
    else:
        masked_primers = getmasked_output.split('\t')
        with open(primerinfo) as i:
            amplicon_data = [line.strip().split('\t') for line in i]
        masked_complete = []
        for primer in masked_primers:
            for amplicon in amplicon_data:
                if primer in amplicon:
                    masked_complete += amplicon
        result = '\t'.join(sorted(set(masked_complete)))
        logging.info('Removing reads primed with any of:')
    return result

