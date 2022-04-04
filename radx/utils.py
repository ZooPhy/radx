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

# def read_ivar(filename):
#     ivar_calls = pd.read_csv(filename, sep='\t')
#     if ivar_calls.empty:
#         return ivar_calls
#     ivar_calls["Variant"] = ivar_calls.apply(lambda row:
#                                     np.nan if row["ALT"][0] == "-" else \
#                                     (str(row["POS"]) +  ":" + str(row["ALT"][1:]) \
#                                      if row["ALT"][0] == "+" else \
#                                      str(row["REF"]) + str(row["POS"]) +  str(row["ALT"])), axis=1)
#     ivar_calls = ivar_calls.dropna(subset=["Variant"]).drop_duplicates(subset=["Variant"])
#     ivar_calls["SOURCE"] = "ivar"
#     return ivar_calls

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
        print()
        print('No affected primer binding sites found!')
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
        print()
        print('Removing reads primed with any of:')
    return result


# def annotate_mutations(row):
#     mutations = []
#     gene = next((gene for gene, (first, last) in GENE_MAP.items() if first <= row["POS"] <= last), "")
#     if gene != "" and row["PASS"] == True and "+" not in row["ALT"]:
#         aa_pos = int((row["POS"] - GENE_MAP[gene][0]) // 3 + 1)
#         mutations.append(f"{gene}:{row['REF_AA']}{str(aa_pos)}{row['ALT_AA']}")
#     return mutations

# def merge_calls(ivar, lofreq):
#     if ivar.empty:
#         return ivar
#     elif lofreq.empty:
#         return lofreq
#     merged = pd.concat([lofreq, ivar[ivar['ALT_FREQ'] > 0.01]], ignore_index=True)
#     # We want to keep mutations that are common in both ivar and lofreq
#     ivar_snps = set([x for x in merged[merged["SOURCE"]=="ivar"]["Variant"] if "-" not in x and ":" not in x])
#     lofreq_snps = set([x for x in merged[merged["SOURCE"]=="lofreq"]["Variant"] if "-" not in x and ":" not in x])
#     # remove mutations that were not found in both ivar and lofreq
#     outer_snps = ivar_snps.symmetric_difference(lofreq_snps)
#     before = len(merged.index)
#     merged = merged[~merged['Variant'].isin(outer_snps)]
#     logging.info("Removed %s symmetric difference mutations on ivar and lofreq", (before-len(merged.index)))
#     # from redundant mutations, keeps the ivar information for amino acid anotation
#     merged = merged.drop_duplicates(subset=["Variant"], keep='last')
#     merged = merged.sort_values(by=["POS"])
#     if merged.empty:
#         return merged
#     merged["Mutation"] = merged.apply(lambda x: annotate_mutations(x), axis=1)
#     return merged

# def filter_merged_calls(merged, min_af=0.05):
#     if "ALT_FREQ" in merged:
#         filtered_merged = merged[~((merged["ALT_FREQ"] < min_af))]
#         if "PASS" in filtered_merged:
#             filtered_merged = filtered_merged[~filtered_merged['PASS'].isin([False])]
#         return filtered_merged
#     else:
#         return merged
