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

def read_ivar(filename):
    ivar_calls = pd.read_csv(filename, sep='\t')
    if ivar_calls.empty:
        return ivar_calls
    ivar_calls["Variant"] = ivar_calls.apply(lambda row:
                                    np.nan if row["ALT"][0] == "-" else \
                                    (str(row["POS"]) +  ":" + str(row["ALT"][1:]) \
                                     if row["ALT"][0] == "+" else \
                                     str(row["REF"]) + str(row["POS"]) +  str(row["ALT"])), axis=1)                                     
    ivar_calls = ivar_calls.dropna(subset=["Variant"]).drop_duplicates(subset=["Variant"])
    return ivar_calls

def read_lofreq(filename):
    lofreq_calls = pd.DataFrame(columns=["REGION", "POS", "REF", "ALT", "QUAL", 
                                         "REF_DP", "REF_RV", "ALT_DP", "ALT_RV",
                                         "ALT_FREQ", "TOTAL_DP"])
    vcf_reader = vcf.Reader(filename=filename)
    for row in vcf_reader:
        lofreq_calls = lofreq_calls.append({"REGION": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": str(row.REF),
                                            "ALT": str(row.ALT[0]),
                                            "QUAL": row.QUAL, 
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
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
    return lofreq_calls

def annotate_mutations(row):
    mutations = []
    gene = next((gene for gene, (first, last) in GENE_MAP.items() if first <= row["POS"] <= last), "")
    print(row)
    if gene != "" and row["PASS"] == True and "+" not in row["ALT"]:
        aa_pos = int((row["POS"] - GENE_MAP[gene][0]) // 3 + 1)
        mutations.append(f"{gene}:{row['REF_AA']}{str(aa_pos)}{row['ALT_AA']}")
    return mutations

def merge_calls(ivar, lofreq):
    if ivar.empty:
        return lofreq
    elif lofreq.empty:
        return ivar
    merged = pd.concat([lofreq, ivar[ivar['ALT_FREQ'] > 0.01]], ignore_index=True)
    merged = merged.drop_duplicates(subset=["Variant"], keep='last') # keeps the ivar information for amino acid anotation
    merged = merged.sort_values(by=["POS"])
    if merged.empty:
        return merged
    merged["Mutation"] = merged.apply(lambda x: annotate_mutations(x), axis=1)
    return merged

def filter_merged_calls(merged, min_af=0.05):
    if "ALT_FREQ" in merged:
        filtered_merged = merged[~((merged["ALT_FREQ"] < min_af))]
        if "PASS" in filtered_merged:
            filtered_merged = filtered_merged[~filtered_merged['PASS'].isin([False])]
        return filtered_merged
    else:
        return merged


