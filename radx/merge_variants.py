import argparse
import vcf
import pandas as pd

genes = {'ORF1a': [266, 13468], 'ORF1b': [13468, 21555], 'S': [21563, 25384], 'ORF3a': [25393, 26220], 'E': [26245, 26472], 'M': [26523, 27191], 'ORF6': [27202, 27387], 'ORF7a': [27394, 27759], 'ORF7b': [27756, 27887], 'ORF8': [27894, 28259], 'N': [28274, 29533], 'ORF10': [29558, 29674]}

def read_ivar(filename):
    ivar_calls = pd.read_csv(filename, sep='\t')
    if ivar_calls.empty:
        return ivar_calls
    ivar_calls["Variant"] = ivar_calls.apply(lambda row:
                                     next if row["ALT"][0] == "-" else \
                                    (str(row["POS"]) + \
                                     ":" + \
                                     str(row["ALT"][1:]) if row["ALT"][0] == "+" \
                                     else \
                                     str(row["REF"]) + \
                                     str(row["POS"]) + \
                                     str(row["ALT"])), axis=1)                                     
    ivar_calls = ivar_calls.drop_duplicates(subset=["Variant"])
    return ivar_calls

def read_lofreq(filename):
    lofreq_calls = pd.DataFrame(columns=["REGION", "POS", "REF", "ALT", "QUAL", 
                                     "REF_DP", "REF_RV", "ALT_DP", "ALT_RV",
                                     "ALT_FREQ", "TOTAL_DP"])
    try:
        vcf_reader = vcf.Reader(filename=filename)
    except FileNotFoundError:
        return None
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
                                         str(row["POS"]) + \
                                         ":" + \
                                         str(row["ALT"][1:]) if len(row["ALT"]) > 1 else \
                                         (str(row["POS"]) + \
                                         "-" + \
                                         str(row["POS"] + len(row["REF"][1:]) + 1) if len(row["REF"]) > 1 else \
                                         str(row["REF"]) + \
                                         str(row["POS"]) + \
                                         str(row["ALT"])), axis=1)
    return lofreq_calls

def annotate_mutations(row):
    mutations = []
    gene = next((gene for gene, (first, last) in genes.items() if first <= row["POS"] <= last), "")
    if gene != "" and row["PASS"] == True:
        aa_pos = int((row["POS"] - genes[gene][0]) // 3 + 1)
        mutations.append(f"{gene}:{row['REF_AA']}{str(aa_pos)}{row['ALT_AA']}")
    return mutations

def merge_calls(ivar, lofreq):
    if ivar.empty:
        return ivar
    elif lofreq.empty:
        return lofreq
    merged = pd.concat([lofreq, ivar[ivar['ALT_FREQ'] > 0.01]], ignore_index=True)
    merged = merged.drop_duplicates(subset=["Variant"], keep='last') # keeps the ivar information for amino acid anotation
    merged = merged.sort_values(by=["POS"])
    if merged.empty:
        return merged
    merged["Mutation"] = merged.apply(lambda x: annotate_mutations(x), axis=1)
    return merged


def filter_merged_calls(merged, min_af):
    filtered_merged = merged[~((merged["ALT_FREQ"] < min_af))]
    filtered_merged = filtered_merged[~filtered_merged['PASS'].isin([False])]
    return filtered_merged


def main():
    parser = argparse.ArgumentParser(description="Merge LoFreq and iVar outputs")
    parser.add_argument("-m", "--min-af", type=float, default=0.05, 
        help="Minimum allele frequency of the variants (between 0 and 1)")
    parser.add_argument("-i","--ivar_input", type=str,
        help="Path to the iVar output .tsv", required=True)
    parser.add_argument("-l", "--lofreq_input", type=str,
        help="Path to the LoFreq output .vcf", required=True)
    parser.add_argument("-o", "--output", type=str, default="output.tsv",
        help="Name of the output file for merged calls")
    args = parser.parse_args()

    ivar_calls = read_ivar(args.ivar_input)
    lofreq_calls = read_lofreq(args.lofreq_input)
    merged_calls = merge_calls(ivar_calls, lofreq_calls)
    if merged_calls.empty:
        exit(0)
    filtered_merged_calls = filter_merged_calls(merged_calls, args.min_af)
    filtered_merged_calls.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
