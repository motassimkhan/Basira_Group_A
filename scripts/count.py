import pandas as pd, glob, re

# 1. Read all masked MAFs
files = glob.glob("mutation/**/*.maf.gz", recursive=True)
maf = pd.concat([pd.read_csv(f, sep="\t", comment="#", low_memory=False,
                             usecols=["Hugo_Symbol","Variant_Classification",
                                      "Tumor_Sample_Barcode"])
                 for f in files])
maf["patient"] = maf["Tumor_Sample_Barcode"].str[:12]

# 2. Keep non-silent coding variants
keep = ["Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
        "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site",
        "Translation_Start_Site","Nonstop_Mutation"]
ns = maf[maf["Variant_Classification"].isin(keep)]

# 3. Cadherin genes
cad_re = r"^(CDH\d+|PCDH.+|FAT[1-4]|DSC[1-3]|DSG[1-4]|CELSR[1-3])$"
ns = ns.assign(is_cad=ns["Hugo_Symbol"].str.match(cad_re))

# 4. Patient universe: mutation patients that also have RNA data
mut_ss = pd.read_csv("mutation/gdc_sample_sheet.2026-10-11.tsv", sep="\t")
rna_ss = pd.read_csv("database/gdc_sample_sheet.2026-10-11.tsv", sep="\t")
patients = set(mut_ss["Sample ID"].str[:12]) & set(rna_ss["Sample ID"].str[:12])

# 5. Flag and count
tp53 = set(ns.loc[ns["Hugo_Symbol"] == "TP53", "patient"])
cad  = set(ns.loc[ns["is_cad"], "patient"])

def grp(p):
    t, c = p in tp53, p in cad
    return ("TP53 + cadherin" if t and c else
            "TP53 only" if t else
            "Cadherin only" if c else "Wild-type (neither)")

print(pd.Series([grp(p) for p in patients]).value_counts())