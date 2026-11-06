import pandas as pd
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


# ============================================================
# STEP 2: SAMPLE QUALITY CONTROL
# Dataset: GSE248821
# Input: Processed RSEM-CPM expression matrix
# ============================================================

expression_file = Path("expression_CPM.csv")
metadata_file = Path("sample_metadata.csv")

# ------------------------------------------------------------
# 1. Create output directories
# ------------------------------------------------------------

Path("PCA").mkdir(exist_ok=True)
Path("Correlation").mkdir(exist_ok=True)
Path("Clustering").mkdir(exist_ok=True)
Path("Expression_Distribution").mkdir(exist_ok=True)

# ------------------------------------------------------------
# 2. Read data
# ------------------------------------------------------------

expression = pd.read_csv(expression_file)
metadata = pd.read_csv(metadata_file)

print("Reading expression matrix...")
print(f"Expression matrix: {expression.shape}")

# Samples are columns
samples = expression.columns.tolist()

if set(metadata["Sample"]) != set(samples):
    raise ValueError("Sample names in expression matrix and metadata do not match.")

# ------------------------------------------------------------
# 3. Log2 transformation
# ------------------------------------------------------------

log_cpm = np.log2(expression + 1)

print("Applied log2(CPM + 1) transformation.")
# ------------------------------------------------------------
# 4. Remove genes with zero variance
# ------------------------------------------------------------

gene_variance = log_cpm.var(axis=1)

log_cpm = log_cpm.loc[
    gene_variance > 0,
    :
]

print(
    f"Genes retained after variance filtering: "
    f"{log_cpm.shape[0]}"
)

# ------------------------------------------------------------
# 5. Select most variable genes for sample QC
# ------------------------------------------------------------

n_variable = min(1000, log_cpm.shape[0])

variance = log_cpm.var(axis=1)

top_genes = variance.nlargest(n_variable).index

qc_matrix = log_cpm.loc[
    top_genes,
    :
]

# Samples × genes
sample_matrix = qc_matrix.T

print(
    f"Top variable genes used for QC: "
    f"{n_variable}"
)
# ------------------------------------------------------------
# 6. PCA
# ------------------------------------------------------------

pca = PCA(n_components=2)

pca_coordinates = pca.fit_transform(sample_matrix)

pca_df = pd.DataFrame(
    pca_coordinates,
    columns=["PC1", "PC2"],
    index=sample_matrix.index
)

pca_df = pca_df.merge(
    metadata,
    left_index=True,
    right_on="Sample"
)

pc1_variance = pca.explained_variance_ratio_[0] * 100
pc2_variance = pca.explained_variance_ratio_[1] * 100

plt.figure(figsize=(9, 7))

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="Condition",
    s=120
)

for _, row in pca_df.iterrows():
    plt.text(
        row["PC1"],
        row["PC2"],
        row["Sample"],
        fontsize=8
    )

plt.xlabel(f"PC1 ({pc1_variance:.2f}% variance)")
plt.ylabel(f"PC2 ({pc2_variance:.2f}% variance)")
plt.title("PCA of GSE248821 Samples")
plt.tight_layout()

plt.savefig(
    "PCA/PCA_plot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

pca_df.to_csv(
    "PCA/PCA_coordinates.csv",
    index=False
)

# ------------------------------------------------------------
# 7. Sample correlation
# ------------------------------------------------------------

correlation = sample_matrix.T.corr(method="pearson")

plt.figure(figsize=(11, 9))

sns.heatmap(
    correlation,
    annot=True,
    fmt=".2f",
    cmap="viridis",
    square=True
)

plt.title("Sample Correlation Heatmap")
plt.tight_layout()

plt.savefig(
    "Correlation/sample_correlation_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

correlation.to_csv(
    "Correlation/sample_correlation.csv"
)

# ------------------------------------------------------------
# 8. Hierarchical clustering
# ------------------------------------------------------------

distance_matrix = pdist(
    sample_matrix,
    metric="euclidean"
)

linkage_matrix = linkage(
    distance_matrix,
    method="complete"
)

plt.figure(figsize=(12, 7))

dendrogram(
    linkage_matrix,
    labels=sample_matrix.index,
    leaf_rotation=90
)

plt.title("Hierarchical Clustering of Samples")
plt.xlabel("Samples")
plt.ylabel("Euclidean distance")
plt.tight_layout()

plt.savefig(
    "Clustering/sample_dendrogram.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# 9. Expression distribution
# ------------------------------------------------------------

plot_df = log_cpm.copy()
plot_df["Gene"] = plot_df.index

long_df = plot_df.melt(
    id_vars="Gene",
    var_name="Sample",
    value_name="log2_CPM"
)

plt.figure(figsize=(14, 7))

sns.boxplot(
    data=long_df,
    x="Sample",
    y="log2_CPM"
)

plt.xticks(rotation=90)
plt.xlabel("Sample")
plt.ylabel("log2(CPM + 1)")
plt.title("Expression Distribution Across Samples")
plt.tight_layout()

plt.savefig(
    "Expression_Distribution/expression_boxplot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# 10. QC summary
# ------------------------------------------------------------

with open("QC_summary.txt", "w") as f:

    f.write("GSE248821 — STEP 2: SAMPLE QUALITY CONTROL\n")
    f.write("===========================================\n\n")

    f.write(f"Number of samples: {len(samples)}\n")
    f.write(f"Genes in expression matrix: {expression.shape[0]}\n")
    f.write(f"Genes used for QC after variance filtering: {log_cpm.shape[0]}\n")
    f.write(f"Top variable genes used for PCA/clustering: {n_variable}\n\n")

    f.write(
        f"PC1 variance explained: {pc1_variance:.2f}%\n"
    )
    f.write(
        f"PC2 variance explained: {pc2_variance:.2f}%\n"
    )

print("\n====================================")
print("STEP 2 QC COMPLETED SUCCESSFULLY")
print("====================================")