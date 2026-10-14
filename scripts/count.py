import pandas as pd, glob

# 1. Read all masked MAFs
files = glob.glob("mutation/**/*.maf.gz", recursive=True)
maf = pd.concat(
    [pd.read_csv(f, sep="\t", comment="#", low_memory=False,
                 usecols=["Hugo_Symbol", "Variant_Classification",
                          "Tumor_Sample_Barcode"])
     for f in files],
    ignore_index=True,
)

# 2. Patient and sample-type IDs from the barcode
maf["patient"] = maf["Tumor_Sample_Barcode"].str[:12]
maf["sample_type"] = maf["Tumor_Sample_Barcode"].str[13:15]

# 3. Keep one tumor type per patient (primary solid tumor = 01).
PRIMARY = {"01"}
maf = maf[maf["sample_type"].isin(PRIMARY)]

# 4. Patient count
patients = set(maf["patient"].unique())

# 5. Non-silent coding variants
keep = ["Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
        "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
        "Translation_Start_Site", "Nonstop_Mutation"]
ns = maf[maf["Variant_Classification"].isin(keep)]

# 6. Cadherin genes
cad_re = r"^(CDH\d+|PCDH.+|FAT[1-4]|DSC[1-3]|DSG[1-4]|CELSR[1-3])$"
is_cad = ns["Hugo_Symbol"].str.match(cad_re, na=False)

tp53 = set(ns.loc[ns["Hugo_Symbol"] == "TP53", "patient"])
cad  = set(ns.loc[is_cad, "patient"])

# 7. Assign each patient to exactly one group
def grp(p):
    t, c = p in tp53, p in cad
    return ("TP53 + cadherin" if t and c else
            "TP53 only"      if t else
            "Cadherin only"  if c else
            "Wild-type (neither)")

groups = pd.Series({p: grp(p) for p in patients}, name="group")
counts = groups.value_counts()
print(counts)

# 8. Sanity checks
assert counts.sum() == len(patients) == maf["patient"].nunique()
assert tp53 <= patients and cad <= patients
print("Total patients:", len(patients))
