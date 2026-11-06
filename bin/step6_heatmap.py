
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

EXPRESSION_FILE = "02_Data_Preparation/expression_CPM.csv"
ANNOTATION_FILE = "02_Data_Preparation/gene_annotation.csv"

INTERACTION_FILE = (
    "04_Differential_Expression/BRD4_IFNg_Interaction/"
    "BRD4_IFNg_interaction_results.csv"
)

OUTPUT_DIR = "07_Visualizations/Heatmaps"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Read input files
expression = pd.read_csv(EXPRESSION_FILE)
annotation = pd.read_csv(ANNOTATION_FILE)
interaction = pd.read_csv(INTERACTION_FILE)

# Validate inputs
if len(expression) != len(annotation):
    raise ValueError("Expression and annotation row counts do not match.")

if "GeneID" not in annotation.columns:
    raise ValueError("GeneID missing from annotation file.")

required = ["GeneID", "logFC", "t", "adj.P.Val"]
missing = [col for col in required if col not in interaction.columns]

if missing:
    raise ValueError(f"Missing interaction columns: {missing}")

if expression.columns.duplicated().any():
    raise ValueError("Duplicate sample columns detected.")

# Add GeneID while preserving STEP1 row order
expression.insert(0, "GeneID", annotation["GeneID"].values)

sample_columns = expression.columns[1:].tolist()

if len(sample_columns) != 13:
    raise ValueError(f"Expected 13 samples, found {len(sample_columns)}.")

# Validate numeric expression values
expression[sample_columns] = expression[sample_columns].apply(
    pd.to_numeric, errors="coerce"
)

if expression[sample_columns].isna().any().any():
    raise ValueError("Missing or non-numeric expression values detected.")

if (expression[sample_columns] < 0).any().any():
    raise ValueError("Negative CPM values detected.")

# Select significant interaction DEGs
interaction["logFC"] = pd.to_numeric(
    interaction["logFC"], errors="coerce"
)
interaction["t"] = pd.to_numeric(
    interaction["t"], errors="coerce"
)
interaction["adj.P.Val"] = pd.to_numeric(
    interaction["adj.P.Val"], errors="coerce"
)

sig = interaction[
    (interaction["adj.P.Val"] < 0.05)
    & (interaction["logFC"].abs() >= 1)
    & interaction["t"].notna()
].copy()

if sig.empty:
    raise ValueError("No significant interaction DEGs found.")

# Rank genes by absolute t-statistic
sig["abs_t"] = sig["t"].abs()

top_genes = (
    sig.sort_values("abs_t", ascending=False)
    .drop_duplicates("GeneID")
    .head(50)["GeneID"]
    .tolist()
)

# Extract expression in ranked gene order
expression = expression.set_index("GeneID")

available_genes = [
    gene for gene in top_genes if gene in expression.index
]

if not available_genes:
    raise ValueError("No selected interaction genes found in expression matrix.")

heatmap_data = expression.loc[available_genes, sample_columns].copy()

# Log transformation
heatmap_data = np.log2(heatmap_data + 1)

# Row-wise Z-score
row_sd = heatmap_data.std(axis=1).replace(0, np.nan)

heatmap_z = heatmap_data.sub(
    heatmap_data.mean(axis=1), axis=0
).div(row_sd, axis=0).fillna(0)

# Use actual sample names as labels
plt.figure(figsize=(15, 13))

sns.heatmap(
    heatmap_z,
    cmap="vlag",
    center=0,
    xticklabels=sample_columns,
    yticklabels=True,
    cbar_kws={"label": "Row Z-score"}
)

plt.title(
    "Top BRD4 × IFNγ Interaction DEGs",
    fontsize=15,
    pad=15
)
plt.xlabel("Samples")
plt.ylabel("Gene ID")
plt.xticks(rotation=45, ha="right", fontsize=8)
plt.yticks(fontsize=7)
plt.tight_layout()

heatmap_file = os.path.join(
    OUTPUT_DIR,
    "BRD4_IFNg_Interaction_Top50_Heatmap.png"
)

plt.savefig(heatmap_file, dpi=300, bbox_inches="tight")
plt.close()

# Save selected gene list
selected = sig[sig["GeneID"].isin(available_genes)].copy()

selected.sort_values(
    "abs_t", ascending=False
).drop(columns=["abs_t"]).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "BRD4_IFNg_Interaction_Top50_Genes.csv"
    ),
    index=False
)

# Save summary
with open(
    os.path.join(OUTPUT_DIR, "Heatmap_summary.txt"), "w"
) as f:
    f.write("BRD4 × IFNγ Interaction Heatmap\n")
    f.write("================================\n")
    f.write(f"Significant interaction DEGs: {len(sig)}\n")
    f.write(f"Genes displayed: {len(available_genes)}\n")
    f.write("Selection: Top 50 by absolute limma t-statistic\n")
    f.write("Transformation: log2(CPM + 1)\n")
    f.write("Scaling: Row-wise Z-score\n")

print("Heatmap completed successfully.")
print(f"Genes displayed: {len(available_genes)}")
print(f"Output: {heatmap_file}")
