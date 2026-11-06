import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# STEP 4: DEG PROCESSING + VOLCANO PLOTS
# Dataset: GSE248821
# Input: limma differential-expression result tables
# ============================================================

if len(sys.argv) != 5:
    print(
        "Usage: python step4_deg_processing.py "
        "<BRD4_results> <IFNg_results> "
        "<KO_IFNg_results> <Interaction_results>"
    )
    sys.exit(1)


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

input_files = {
    "BRD4_Effect": Path(sys.argv[1]),
    "IFNg_Effect": Path(sys.argv[2]),
    "IFNg_Response_KO": Path(sys.argv[3]),
    "BRD4_IFNg_Interaction": Path(sys.argv[4]),
}


# ------------------------------------------------------------
# 2. Create output directories
# ------------------------------------------------------------

Path("DEG_Results/Significant").mkdir(
    parents=True,
    exist_ok=True
)

Path("DEG_Results/Upregulated").mkdir(
    parents=True,
    exist_ok=True
)

Path("DEG_Results/Downregulated").mkdir(
    parents=True,
    exist_ok=True
)

Path("Volcano").mkdir(
    parents=True,
    exist_ok=True
)


# ------------------------------------------------------------
# 3. Processing function
# ------------------------------------------------------------

def process_contrast(name, input_file):

    print(f"\nProcessing: {name}")
    print(f"Input: {input_file}")

    df = pd.read_csv(input_file)

    required_columns = [
        "logFC",
        "adj.P.Val",
        "Significant"
    ]

    missing = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{name}: missing required columns: {missing}"
        )

    # --------------------------------------------------------
    # Significant genes
    # --------------------------------------------------------

    significant = df[
        (df["Significant"] == "YES")
    ].copy()

    # --------------------------------------------------------
    # Upregulated genes
    # logFC >= 1
    # --------------------------------------------------------

    upregulated = significant[
        significant["logFC"] >= 1
    ].copy()

    # --------------------------------------------------------
    # Downregulated genes
    # logFC <= -1
    # --------------------------------------------------------

    downregulated = significant[
        significant["logFC"] <= -1
    ].copy()

    # --------------------------------------------------------
    # Save DEG tables
    # --------------------------------------------------------

    significant.to_csv(
        f"DEG_Results/Significant/{name}_significant.csv",
        index=False
    )

    upregulated.to_csv(
        f"DEG_Results/Upregulated/{name}_upregulated.csv",
        index=False
    )

    downregulated.to_csv(
        f"DEG_Results/Downregulated/{name}_downregulated.csv",
        index=False
    )

    # --------------------------------------------------------
    # Volcano plot preparation
    # --------------------------------------------------------

    plot_df = df.copy()

    # Avoid log10(0)
    plot_df["adj.P.Val"] = plot_df["adj.P.Val"].clip(
        lower=np.finfo(float).tiny
    )

    plot_df["neg_log10_FDR"] = -np.log10(
        plot_df["adj.P.Val"]
    )

    plot_df["Category"] = "Not Significant"

    plot_df.loc[
        (plot_df["Significant"] == "YES") &
        (plot_df["logFC"] >= 1),
        "Category"
    ] = "Upregulated"

    plot_df.loc[
        (plot_df["Significant"] == "YES") &
        (plot_df["logFC"] <= -1),
        "Category"
    ] = "Downregulated"

    # --------------------------------------------------------
    # Volcano plot
    # --------------------------------------------------------

    plt.figure(figsize=(10, 7))

    sns.scatterplot(
        data=plot_df,
        x="logFC",
        y="neg_log10_FDR",
        hue="Category",
        s=35,
        alpha=0.7
    )

    plt.axvline(
        x=1,
        linestyle="--"
    )

    plt.axvline(
        x=-1,
        linestyle="--"
    )

    plt.axhline(
        y=-np.log10(0.05),
        linestyle="--"
    )

    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(FDR)")
    plt.title(f"Volcano Plot — {name}")

    plt.tight_layout()

    plt.savefig(
        f"Volcano/{name}_volcano.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    #--------------------------------------------------------
    # MA plot
    # --------------------------------------------------------

    if "AveExpr" not in df.columns:
        raise ValueError(
            f"{name}: AveExpr column missing from limma results."
        )

    plt.figure(figsize=(10, 7))

    ma_df = df.copy()

    ma_df["Category"] = "Not Significant"

    ma_df.loc[
        (ma_df["Significant"] == "YES") &
        (ma_df["logFC"] >= 1),
        "Category"
    ] = "Upregulated"

    ma_df.loc[
        (ma_df["Significant"] == "YES") &
        (ma_df["logFC"] <= -1),
        "Category"
    ] = "Downregulated"

    sns.scatterplot(
        data=ma_df,
        x="AveExpr",
        y="logFC",
        hue="Category",
        s=30,
        alpha=0.7
    )

    plt.axhline(
        y=0,
        linestyle="--"
    )

    plt.xlabel("Average Expression (log2)")
    plt.ylabel("log2 Fold Change")
    plt.title(f"MA Plot — {name}")

    plt.tight_layout()

    plt.savefig(
        f"Volcano/{name}_MA.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    return {
        "Contrast": name,
        "Total genes": len(df),
        "Significant genes": len(significant),
        "Upregulated": len(upregulated),
        "Downregulated": len(downregulated)
    }


# ------------------------------------------------------------
# 4. Process all contrasts
# ------------------------------------------------------------

summary = []

for name, input_file in input_files.items():

    summary.append(
        process_contrast(
            name,
            input_file
        )
    )


# ------------------------------------------------------------
# 5. Save summary
# ------------------------------------------------------------

summary_df = pd.DataFrame(summary)

summary_df.to_csv(
    "DEG_Results/DEG_summary.csv",
    index=False
)

# ------------------------------------------------------------
# 6. Final report
# ------------------------------------------------------------

with open(
    "DEG_Results/DEG_processing_summary.txt",
    "w"
) as f:

    f.write(
        "GSE248821 — STEP 4: DEG PROCESSING\n"
        "===================================\n\n"
    )

    f.write(
        "Significance threshold:\n"
        "FDR < 0.05\n"
        "|logFC| >= 1\n\n"
    )

    f.write(
        summary_df.to_string(index=False)
    )

print("\n====================================")
print("STEP 4 COMPLETED SUCCESSFULLY")
print("====================================\n")

print(summary_df)