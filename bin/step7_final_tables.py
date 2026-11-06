
import os
import shutil
import pandas as pd

OUT = "08_Final_Tables"

DEG_BASE = "05_DEG_Results/DEG_Results"
PATHWAY_BASE = "06_Pathway_Analysis/Pathway_Analysis"

CONTRASTS = [
    "BRD4_Effect",
    "IFNg_Effect",
    "IFNg_Response_KO",
    "BRD4_IFNg_Interaction"
]

CATEGORIES = [
    "Significant",
    "Upregulated",
    "Downregulated"
]

PATHWAY_METHODS = [
    "GO",
    "KEGG",
    "Reactome",
    "GSEA"
]


def copy_required(source, destination):
    if not os.path.isfile(source):
        raise FileNotFoundError(
            f"Required input file not found: {source}"
        )

    os.makedirs(os.path.dirname(destination), exist_ok=True)
    shutil.copy2(source, destination)


# ---------------------------------------------------------
# 1. Create final output directories
# ---------------------------------------------------------
os.makedirs(OUT, exist_ok=True)

for folder in ["DEG", "Interaction", "Pathways", "GSEA"]:
    os.makedirs(os.path.join(OUT, folder), exist_ok=True)


# ---------------------------------------------------------
# 2. Copy DEG summaries
# ---------------------------------------------------------
copy_required(
    f"{DEG_BASE}/DEG_summary.csv",
    f"{OUT}/DEG/DEG_summary.csv"
)

copy_required(
    f"{DEG_BASE}/DEG_processing_summary.txt",
    f"{OUT}/DEG/DEG_processing_summary.txt"
)


# ---------------------------------------------------------
# 3. Copy DEG tables for all four contrasts
# ---------------------------------------------------------
for contrast in CONTRASTS:

    for category in CATEGORIES:

        filename = f"{contrast}_{category.lower()}.csv"

        source = os.path.join(
            DEG_BASE,
            category,
            filename
        )

        destination = os.path.join(
            OUT,
            "DEG",
            contrast,
            filename
        )

        copy_required(source, destination)


# ---------------------------------------------------------
# 4. Copy main interaction DEG tables
# ---------------------------------------------------------
for category in CATEGORIES:

    filename = (
        f"BRD4_IFNg_Interaction_{category.lower()}.csv"
    )

    source = os.path.join(
        DEG_BASE,
        category,
        filename
    )

    destination = os.path.join(
        OUT,
        "Interaction",
        filename
    )

    copy_required(source, destination)


# ---------------------------------------------------------
# 5. Copy pathway and GSEA results for all contrasts
# ---------------------------------------------------------
for method in PATHWAY_METHODS:

    for contrast in CONTRASTS:

        source_dir = os.path.join(
            PATHWAY_BASE,
            method,
            contrast
        )

        if not os.path.isdir(source_dir):
            raise FileNotFoundError(
                f"Pathway output directory not found: {source_dir}"
            )

        if method == "GSEA":
            destination_dir = os.path.join(
                OUT, "GSEA", contrast
            )
        else:
            destination_dir = os.path.join(
                OUT, "Pathways", method, contrast
            )

        shutil.copytree(
            source_dir,
            destination_dir,
            dirs_exist_ok=True
        )


# ---------------------------------------------------------
# 6. Generate final results summary
# ---------------------------------------------------------
summary_rows = []

# DEG counts
for contrast in CONTRASTS:

    filename = f"{contrast}_significant.csv"

    filepath = os.path.join(
        OUT, "DEG", contrast, filename
    )

    df = pd.read_csv(filepath)

    summary_rows.append({
        "Contrast": contrast,
        "Analysis": "Differential Expression",
        "Result": "Significant DEGs",
        "Count": len(df)
    })


# Pathway and GSEA result counts
for method in PATHWAY_METHODS:

    for contrast in CONTRASTS:

        if method == "GSEA":
            base_dir = os.path.join(
                OUT, "GSEA", contrast
            )
        else:
            base_dir = os.path.join(
                OUT, "Pathways", method, contrast
            )

        for root, _, files in os.walk(base_dir):

            for filename in files:

                if not filename.endswith(".csv"):
                    continue

                filepath = os.path.join(root, filename)
                df = pd.read_csv(filepath)

                summary_rows.append({
                    "Contrast": contrast,
                    "Analysis": method,
                    "Result": filename,
                    "Count": len(df)
                })


pd.DataFrame(summary_rows).to_csv(
    os.path.join(
        OUT, "final_results_summary.csv"
    ),
    index=False
)


# ---------------------------------------------------------
# 7. Save methodology note
# ---------------------------------------------------------
with open(
    os.path.join(OUT, "final_tables_method.txt"),
    "w"
) as f:

    f.write(
        "Final tables generated from differential expression "
        "and pathway analysis outputs.\n\n"
    )

    f.write("Differential expression contrasts:\n")

    for contrast in CONTRASTS:
        f.write(f"- {contrast}\n")

    f.write("\nPathway analysis methods:\n")

    for method in PATHWAY_METHODS:
        f.write(f"- {method}\n")

    f.write(
        "\nDEG threshold: FDR < 0.05 and |logFC| >= 1.\n"
    )

    f.write(
        "GSEA ranking metric: limma t-statistic.\n"
    )

    f.write(
        "Primary biological focus: BRD4 x IFNg interaction.\n"
    )


print("Final tables generated successfully.")
print("All four contrasts included.")
print("GO, KEGG, Reactome and GSEA results copied.")
print(f"Output directory: {OUT}")
