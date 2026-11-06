
import sys
from pathlib import Path
import pandas as pd

# ============================================================
# STEP 1: DATA PREPARATION
# Dataset: GSE248821
# Input: Processed RSEM-CPM Excel file
# ============================================================

if len(sys.argv) != 2:
    print("Usage: python step1_data_preparation.py <input_excel>")
    sys.exit(1)

input_excel = Path(sys.argv[1])

print(f"Input file: {input_excel}")
print("Reading GEO Excel file...")

# 1. Read GEO processed expression data
df = pd.read_excel(
    input_excel,
    sheet_name="gene_RSEM_CPM"
)

print(f"Original dataset: {df.shape[0]} rows × {df.shape[1]} columns")

# 2. Define the 13 biological samples
sample_columns = [
    "control, rep1",
    "control, rep2",
    "control, rep3",
    "control, rep4",
    "Brd4 KO, rep1",
    "Brd4 KO, rep2",
    "Brd4 KO, rep3",
    "control, IFNg, rep1",
    "control, IFNg, rep2",
    "control, IFNg, rep3",
    "Brd4 KO, IFNg, rep1",
    "Brd4 KO, IFNg, rep2",
    "Brd4 KO, IFNg, rep3"
]

missing_samples = [
    sample for sample in sample_columns
    if sample not in df.columns
]

if missing_samples:
    raise ValueError(
        f"Missing expected sample columns: {missing_samples}"
    )

# 3. Gene annotation
gene_info = df[["GeneID", "TranscriptID"]].copy()

# 4. Expression matrix
expression = df[sample_columns].copy()

expression = expression.apply(
    pd.to_numeric,
    errors="coerce"
)

# 5. Validate original expression data BEFORE filtering
missing_values = int(expression.isna().sum().sum())

if missing_values > 0:
    raise ValueError(
        f"Found {missing_values} missing/non-numeric expression values."
    )

if (expression < 0).any().any():
    raise ValueError("Negative expression values detected.")

duplicate_genes = int(
    gene_info["GeneID"].duplicated().sum()
)

if duplicate_genes > 0:
    raise ValueError(
        f"Found {duplicate_genes} duplicate GeneIDs."
    )

# 6. Original dataset statistics
original_gene_count = expression.shape[0]

all_zero_genes = int(
    (expression == 0).all(axis=1).sum()
)

expressed_genes = int(
    (expression > 0).any(axis=1).sum()
)

# 7. Filter low-expression genes
# Retain genes with CPM >= 1 in at least 3 samples

keep_genes = (
    expression.ge(1).sum(axis=1) >= 3
)

expression = expression.loc[keep_genes].copy()
gene_info = gene_info.loc[keep_genes].copy()

filtered_gene_count = expression.shape[0]

print("Genes retained after filtering:", filtered_gene_count)

# 8. Sample metadata
metadata = []

for sample in sample_columns:

    brd4 = "KO" if "Brd4 KO" in sample else "Control"
    ifng = "IFNg" if "IFNg" in sample else "No_IFNg"

    metadata.append({
        "Sample": sample,
        "BRD4": brd4,
        "IFNg": ifng,
        "Condition": f"{brd4}_{ifng}"
    })

metadata = pd.DataFrame(metadata)

# 9. Final consistency checks
if expression.shape[0] != gene_info.shape[0]:
    raise ValueError("Expression and annotation row counts do not match.")

if expression.columns.tolist() != sample_columns:
    raise ValueError("Expression sample columns are not in expected order.")

unique_genes = int(gene_info["GeneID"].nunique())

if unique_genes != filtered_gene_count:
    raise ValueError("GeneID uniqueness check failed after filtering.")

# 10. Save outputs
gene_info.to_csv("gene_annotation.csv", index=False)
expression.to_csv("expression_CPM.csv", index=False)
metadata.to_csv("sample_metadata.csv", index=False)

# 11. Summary
summary = f"""GSE248821 — STEP 1: DATA PREPARATION
=====================================

Original dataset:
{df.shape[0]} rows × {df.shape[1]} columns

Original number of genes:
{original_gene_count}

Number of samples:
{len(sample_columns)}

All-zero genes in original dataset:
{all_zero_genes}

Genes with expression > 0 in original dataset:
{expressed_genes}

Filtering criterion:
CPM >= 1 in at least 3 samples

Genes retained after filtering:
{filtered_gene_count}

Unique GeneIDs after filtering:
{unique_genes}

Duplicate GeneIDs:
{duplicate_genes}

Missing expression values:
{missing_values}

Sample groups:
Control = 4
BRD4 KO = 3
Control + IFNg = 3
BRD4 KO + IFNg = 3

Expression data:
Processed RSEM-CPM
"""

with open("step1_summary.txt", "w") as file:
    file.write(summary)

print("\n====================================")
print("STEP 1 COMPLETED SUCCESSFULLY")
print("====================================")
print(summary)

print("Output files:")
print("gene_annotation.csv")
print("expression_CPM.csv")
print("sample_metadata.csv")
print("step1_summary.txt")
