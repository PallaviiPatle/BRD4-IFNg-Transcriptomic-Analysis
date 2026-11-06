
# ============================================================
# STEP 5: PATHWAY AND FUNCTIONAL ENRICHMENT ANALYSIS
# Dataset: GSE248821
# Organism: Mus musculus (Mouse)
#
# Input: limma differential-expression result tables
# ORA:
#   1. GO Biological Process
#   2. KEGG
#   3. Reactome
#
# GSEA:
#   1. GO Biological Process
#   2. KEGG
#   3. Reactome
#
# All four experimental contrasts are analyzed separately.
# ============================================================

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({

    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(ReactomePA)
    library(AnnotationDbi)
    library(dplyr)
    library(ggplot2)

})

options(stringsAsFactors = FALSE)
set.seed(1234)

# ------------------------------------------------------------
# 1. Input arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Expected four differential-expression CSV files.")
}

names(args) <- c(
    "BRD4_Effect",
    "IFNg_Effect",
    "IFNg_Response_KO",
    "BRD4_IFNg_Interaction"
)

# ------------------------------------------------------------
# 2. Validate input files
# ------------------------------------------------------------

for (file in args) {

    if (!file.exists(file)) {
        stop("Input file not found: ", file)
    }

}

# ------------------------------------------------------------
# 3. Output directories
# ------------------------------------------------------------

base_out <- "Pathway_Analysis"

dir.create(
    base_out,
    recursive = TRUE,
    showWarnings = FALSE
)

for (analysis in names(args)) {

    for (database in c("GO", "KEGG", "Reactome", "GSEA")) {

        dir.create(
            file.path(base_out, database, analysis),
            recursive = TRUE,
            showWarnings = FALSE
        )

    }

}

# ------------------------------------------------------------
# 4. Helper function to save enrichment results
# ------------------------------------------------------------

save_enrichment <- function(
    result,
    csv_file,
    plot_file,
    title
) {

    if (is.null(result)) {

        write.csv(
            data.frame(),
            csv_file,
            row.names = FALSE
        )

        message("No enrichment result: ", title)

        return(invisible(NULL))
    }

    result_df <- as.data.frame(result)

    write.csv(
        result_df,
        csv_file,
        row.names = FALSE
    )

    if (nrow(result_df) == 0) {

        message("No significant terms: ", title)

        return(invisible(NULL))
    }

    p <- dotplot(
        result,
        showCategory = min(15, nrow(result_df)),
        title = title
    )

    ggsave(
        filename = plot_file,
        plot = p,
        width = 9,
        height = 7,
        dpi = 300
    )

}

# ------------------------------------------------------------
# 5. Helper function for GSEA results
# ------------------------------------------------------------

save_gsea <- function(
    result,
    csv_file,
    plot_file,
    title
) {

    if (is.null(result)) {

        write.csv(
            data.frame(),
            csv_file,
            row.names = FALSE
        )

        message("GSEA returned no result: ", title)

        return(invisible(NULL))
    }

    result_df <- as.data.frame(result)

    write.csv(
        result_df,
        csv_file,
        row.names = FALSE
    )

    if (nrow(result_df) == 0) {

        message("No significant GSEA pathways: ", title)

        return(invisible(NULL))
    }

    p <- dotplot(
        result,
        showCategory = min(15, nrow(result_df)),
        title = title
    )

    ggsave(
        filename = plot_file,
        plot = p,
        width = 9,
        height = 7,
        dpi = 300
    )

}

# ------------------------------------------------------------
# 6. Process all four contrasts
# ------------------------------------------------------------

summary_list <- list()

for (analysis in names(args)) {

    input_file <- args[[analysis]]

    message("\n========================================")
    message("Processing contrast: ", analysis)
    message("Input: ", input_file)
    message("========================================")

    # --------------------------------------------------------
    # 6A. Read differential-expression results
    # --------------------------------------------------------

    df <- read.csv(
        input_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    required_cols <- c(
        "GeneID",
        "logFC",
        "t",
        "adj.P.Val",
        "Significant"
    )

    missing_cols <- setdiff(
        required_cols,
        colnames(df)
    )

    if (length(missing_cols) > 0) {

        stop(
            "Missing columns in ",
            input_file,
            ": ",
            paste(missing_cols, collapse = ", ")
        )

    }

    if (nrow(df) == 0) {
        stop("No genes found in input file: ", input_file)
    }

    if (anyDuplicated(df$GeneID)) {
        stop("Duplicate GeneIDs found in: ", input_file)
    }

    # --------------------------------------------------------
    # 6B. Validate and normalize Ensembl gene IDs
    # --------------------------------------------------------

    df$GeneID <- as.character(df$GeneID)

    # Remove Ensembl version suffix if present.
    # Example: ENSMUSG00000000001.5
    # becomes: ENSMUSG00000000001

    df$ENSEMBL <- sub(
        "\\.[0-9]+$",
        "",
        df$GeneID
    )

    if (any(is.na(df$ENSEMBL)) ||
        any(df$ENSEMBL == "")) {

        stop("Missing or invalid Ensembl gene IDs.")
    }

    # --------------------------------------------------------
    # 6C. Ensembl to Entrez and gene symbol mapping
    # --------------------------------------------------------

    mapping <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = unique(df$ENSEMBL),
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "ENSEMBL"
    )

    mapping <- mapping %>%
        filter(
            !is.na(ENSEMBL),
            !is.na(ENTREZID),
            ENTREZID != ""
        ) %>%
        distinct(
            ENSEMBL,
            ENTREZID,
            .keep_all = TRUE
        )

    # Remove ambiguous Ensembl-to-Entrez mappings.
    # Retain only Ensembl IDs mapping to one unique Entrez ID.

    mapping_counts <- mapping %>%
        count(ENSEMBL, name = "mapping_count")

    valid_ensembl <- mapping_counts$ENSEMBL[
        mapping_counts$mapping_count == 1
    ]

    mapping <- mapping %>%
        filter(ENSEMBL %in% valid_ensembl)

    if (nrow(mapping) == 0) {
        stop(
            "No Ensembl IDs could be mapped to Entrez IDs. ",
            "Verify the GeneID format and organism annotation."
        )
    }

    # --------------------------------------------------------
    # 6D. Merge expression results with gene annotation
    # --------------------------------------------------------

    df_map <- merge(
        df,
        mapping,
        by = "ENSEMBL",
        all = FALSE,
        sort = FALSE
    )

    if (nrow(df_map) == 0) {
        stop("No genes remained after ID mapping.")
    }

    # --------------------------------------------------------
    # 6E. Define ORA background universe
    # --------------------------------------------------------

    # Background includes all tested genes successfully
    # mapped to Entrez IDs.

    universe <- unique(
        as.character(df_map$ENTREZID)
    )

    # --------------------------------------------------------
    # 6F. Identify significant DEGs
    # --------------------------------------------------------

    sig_df <- df_map[
        !is.na(df_map$adj.P.Val) &
        !is.na(df_map$logFC) &
        df_map$adj.P.Val < 0.05 &
        abs(df_map$logFC) >= 1,
    ]

    sig_genes <- unique(
        as.character(sig_df$ENTREZID)
    )

    up_genes <- unique(
        as.character(
            sig_df$ENTREZID[sig_df$logFC > 0]
        )
    )

    down_genes <- unique(
        as.character(
            sig_df$ENTREZID[sig_df$logFC < 0]
        )
    )

    message("Tested genes: ", nrow(df))
    message("Mapped genes: ", length(universe))
    message("Significant DEGs: ", length(sig_genes))

    # --------------------------------------------------------
    # 7. GO Biological Process ORA
    # --------------------------------------------------------

    go_result <- tryCatch({

        if (length(sig_genes) == 0) {
            return(NULL)
        }

        enrichGO(
            gene = sig_genes,
            universe = universe,
            OrgDb = org.Mm.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.20,
            minGSSize = 10,
            maxGSSize = 500,
            readable = TRUE
        )

    }, error = function(e) {

        message("GO ORA failed: ", e$message)

        NULL

    })

    save_enrichment(
        go_result,
        file.path(
            base_out,
            "GO",
            analysis,
            "GO_BP_enrichment.csv"
        ),
        file.path(
            base_out,
            "GO",
            analysis,
            "GO_BP_dotplot.png"
        ),
        paste0(analysis, " - GO Biological Process")
    )

    # --------------------------------------------------------
    # 8. KEGG ORA
    # --------------------------------------------------------

    kegg_result <- tryCatch({

        if (length(sig_genes) == 0) {
            return(NULL)
        }

        enrichKEGG(
            gene = sig_genes,
            universe = universe,
            organism = "mmu",
            keyType = "ncbi-geneid",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.20,
            minGSSize = 10,
            maxGSSize = 500
        )

    }, error = function(e) {

        message("KEGG ORA failed: ", e$message)

        NULL

    })

    save_enrichment(
        kegg_result,
        file.path(
            base_out,
            "KEGG",
            analysis,
            "KEGG_enrichment.csv"
        ),
        file.path(
            base_out,
            "KEGG",
            analysis,
            "KEGG_dotplot.png"
        ),
        paste0(analysis, " - KEGG")
    )

    # --------------------------------------------------------
    # 9. Reactome ORA
    # --------------------------------------------------------

    reactome_result <- tryCatch({

        if (length(sig_genes) == 0) {
            return(NULL)
        }

        enrichPathway(
            gene = sig_genes,
            universe = universe,
            organism = "mouse",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.20,
            minGSSize = 10,
            maxGSSize = 500,
            readable = TRUE
        )

    }, error = function(e) {

        message("Reactome ORA failed: ", e$message)

        NULL

    })

    save_enrichment(
        reactome_result,
        file.path(
            base_out,
            "Reactome",
            analysis,
            "Reactome_enrichment.csv"
        ),
        file.path(
            base_out,
            "Reactome",
            analysis,
            "Reactome_dotplot.png"
        ),
        paste0(analysis, " - Reactome")
    )

    # --------------------------------------------------------
    # 10. Prepare ranked gene list for GSEA
    # --------------------------------------------------------

    # Use ALL tested and mapped genes.
    # Ranking metric: limma moderated t-statistic.

    rank_df <- df_map[
        !is.na(df_map$t) &
        is.finite(df_map$t) &
        !is.na(df_map$ENTREZID),
    ]

    # Retain one row per Entrez ID.
    # For duplicate Entrez IDs, retain the largest
    # absolute moderated t-statistic.

    rank_df <- rank_df[
        order(
            abs(rank_df$t),
            decreasing = TRUE
        ),
    ]

    rank_df <- rank_df[
        !duplicated(rank_df$ENTREZID),
    ]

    gene_list <- rank_df$t

    names(gene_list) <- as.character(
        rank_df$ENTREZID
    )

    gene_list <- sort(
        gene_list,
        decreasing = TRUE
    )

    if (length(gene_list) < 10) {
        stop(
            "Insufficient mapped genes for GSEA in ",
            analysis
        )
    }

    # --------------------------------------------------------
    # 11. GSEA - GO Biological Process
    # --------------------------------------------------------

    gsea_go <- tryCatch({

        gseGO(
            geneList = gene_list,
            OrgDb = org.Mm.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            verbose = FALSE,
            seed = TRUE
        )

    }, error = function(e) {

        message("GSEA GO failed: ", e$message)

        NULL

    })

    save_gsea(
        gsea_go,
        file.path(
            base_out,
            "GSEA",
            analysis,
            "GSEA_GO_BP.csv"
        ),
        file.path(
            base_out,
            "GSEA",
            analysis,
            "GSEA_GO_BP_dotplot.png"
        ),
        paste0(
            analysis,
            " - GSEA GO Biological Process"
        )
    )

    # --------------------------------------------------------
    # 12. GSEA - KEGG
    # --------------------------------------------------------

    gsea_kegg <- tryCatch({

        gseKEGG(
            geneList = gene_list,
            organism = "mmu",
            keyType = "ncbi-geneid",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            verbose = FALSE,
            seed = TRUE
        )

    }, error = function(e) {

        message("GSEA KEGG failed: ", e$message)

        NULL

    })

    save_gsea(
        gsea_kegg,
        file.path(
            base_out,
            "GSEA",
            analysis,
            "GSEA_KEGG.csv"
        ),
        file.path(
            base_out,
            "GSEA",
            analysis,
            "GSEA_KEGG_dotplot.png"
        ),
        paste0(analysis, " - GSEA KEGG")
    )

    # --------------------------------------------------------
    # 13. GSEA - Reactome
    # --------------------------------------------------------

    gsea_reactome <- tryCatch({

        gsePathway(
            geneList = gene_list,
            organism = "mouse",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            verbose = FALSE,
            seed = TRUE
        )

    }, error = function(e) {

        message("GSEA Reactome failed: ", e$message)

        NULL

    })

    save_gsea(
        gsea_reactome,
        file.path(
            base_out,
            "GSEA",
            analysis,
            "GSEA_Reactome.csv"
        ),
        file.path(
            base_out,
            "GSEA",
            analysis,
            "GSEA_Reactome_dotplot.png"
        ),
        paste0(analysis, " - GSEA Reactome")
    )

    # --------------------------------------------------------
    # 14. Contrast summary
    # --------------------------------------------------------

    summary_list[[analysis]] <- data.frame(

        Contrast = analysis,

        Tested_genes = nrow(df),

        Mapped_genes = length(universe),

        Mapping_percentage = round(
            100 * length(universe) / nrow(df),
            2
        ),

        Significant_DEGs = length(sig_genes),

        Upregulated_DEGs = length(up_genes),

        Downregulated_DEGs = length(down_genes),

        GSEA_ranked_genes = length(gene_list)

    )

}

# ------------------------------------------------------------
# 15. Save final pathway analysis summary
# ------------------------------------------------------------

summary_df <- bind_rows(summary_list)

write.csv(
    summary_df,
    file.path(
        base_out,
        "pathway_analysis_summary.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 16. Save methodology report
# ------------------------------------------------------------

methodology <- c(

    "GSE248821 - STEP 5: PATHWAY ANALYSIS",

    "======================================",

    "",

    "Organism: Mus musculus",

    "Gene identifier mapping: Ensembl to Entrez",

    "Annotation database: org.Mm.eg.db",

    "",

    "OVER-REPRESENTATION ANALYSIS (ORA)",

    "GO Biological Process",

    "KEGG",

    "Reactome",

    "",

    "GENE SET ENRICHMENT ANALYSIS (GSEA)",

    "GO Biological Process",

    "KEGG",

    "Reactome",

    "",

    "ORA significance criteria:",

    "FDR < 0.05 and |logFC| >= 1",

    "",

    "ORA background:",

    "All tested genes successfully mapped to Entrez IDs",

    "",

    "GSEA ranking:",

    "All mapped tested genes ranked by limma moderated t-statistic",

    "",

    "Multiple testing correction: Benjamini-Hochberg",

    "",

    "All four differential-expression contrasts were analyzed separately."

)

writeLines(
    methodology,
    file.path(
        base_out,
        "Pathway_analysis_method.txt"
    )
)

# ------------------------------------------------------------
# 17. Final status
# ------------------------------------------------------------

cat("\n========================================\n")
cat("STEP 5 PATHWAY ANALYSIS COMPLETED\n")
cat("========================================\n\n")

print(summary_df)

message("\nAll four contrasts processed.")

