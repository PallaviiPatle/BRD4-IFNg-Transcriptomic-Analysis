
# ============================================================
# STEP 3: DIFFERENTIAL EXPRESSION ANALYSIS
# Dataset: GSE248821
# Input: Processed RSEM-CPM expression matrix
# Statistical method: limma-trend
# Design: BRD4 * IFNg (2 x 2 factorial design)
# ============================================================

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({
    library(limma)
})

# ------------------------------------------------------------
# 1. Read input files
# ------------------------------------------------------------

expression <- read.csv(
    "expression_CPM.csv",
    check.names = FALSE
)

metadata <- read.csv(
    "sample_metadata.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

annotation <- read.csv(
    "gene_annotation.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat("Expression matrix:", nrow(expression), "genes x",
    ncol(expression), "samples\n")

cat("Metadata:", nrow(metadata), "samples\n")
cat("Annotation:", nrow(annotation), "genes\n")

# ------------------------------------------------------------
# 2. Validate input dimensions and sample matching
# ------------------------------------------------------------

if (ncol(expression) != 13) {
    stop("Expected 13 samples in expression matrix.")
}

if (nrow(expression) != nrow(annotation)) {
    stop("Expression and annotation row counts do not match.")
}

sample_names <- colnames(expression)

if (anyDuplicated(sample_names)) {
    stop("Duplicate sample names in expression matrix.")
}

if (anyDuplicated(metadata$Sample)) {
    stop("Duplicate sample names in metadata.")
}

if (!setequal(metadata$Sample, sample_names)) {
    stop("Sample names in expression matrix and metadata do not match.")
}

metadata <- metadata[
    match(sample_names, metadata$Sample),
    ,
    drop = FALSE
]

if (any(is.na(metadata$Sample))) {
    stop("Metadata matching failed.")
}

if (anyDuplicated(annotation$GeneID)) {
    stop("Duplicate GeneIDs found in annotation.")
}

if (any(is.na(annotation$GeneID))) {
    stop("Missing GeneIDs found in annotation.")
}

# STEP 1 filters expression and annotation using the same
# row mask. Their original row order must therefore match.

cat("Input validation passed.\n")

# ------------------------------------------------------------
# 3. Convert expression data into numeric matrix
# ------------------------------------------------------------

expr_matrix <- as.matrix(expression)

mode(expr_matrix) <- "numeric"

rownames(expr_matrix) <- as.character(annotation$GeneID)

if (any(!is.finite(expr_matrix))) {
    stop("Non-finite expression values detected.")
}

if (any(expr_matrix < 0)) {
    stop("Negative expression values detected.")
}

# ------------------------------------------------------------
# 4. Validate expression filtering
# ------------------------------------------------------------

genes_before_filtering <- nrow(expr_matrix)

keep <- rowSums(expr_matrix >= 1) >= 3

expr_matrix <- expr_matrix[keep, , drop = FALSE]

annotation <- annotation[keep, , drop = FALSE]

genes_after_filtering <- nrow(expr_matrix)

cat("Genes before filtering:", genes_before_filtering, "\n")
cat("Genes after filtering:", genes_after_filtering, "\n")

if (genes_after_filtering == 0) {
    stop("No genes remain after expression filtering.")
}

if (!identical(
    rownames(expr_matrix),
    as.character(annotation$GeneID)
)) {
    stop("Gene annotation order does not match expression matrix.")
}

# ------------------------------------------------------------
# 5. Log2 transformation
# ------------------------------------------------------------

# Input contains processed RSEM-CPM values,
# not raw integer counts.

logCPM <- log2(expr_matrix + 1)

if (any(!is.finite(logCPM))) {
    stop("Non-finite values detected after log transformation.")
}

cat("Applied log2(CPM + 1) transformation.\n")

# ------------------------------------------------------------
# 6. Define experimental factors
# ------------------------------------------------------------

metadata$BRD4 <- factor(
    metadata$BRD4,
    levels = c("Control", "KO")
)

metadata$IFNg <- factor(
    metadata$IFNg,
    levels = c("No_IFNg", "IFNg")
)

if (any(is.na(metadata$BRD4)) ||
    any(is.na(metadata$IFNg))) {
    stop("Unexpected BRD4 or IFNg metadata values.")
}

# Check experimental group sizes
group_counts <- table(
    metadata$BRD4,
    metadata$IFNg
)

cat("\nExperimental group sizes:\n")
print(group_counts)

# ------------------------------------------------------------
# 7. Create factorial design matrix
# ------------------------------------------------------------

design <- model.matrix(
    ~ BRD4 * IFNg,
    data = metadata
)

colnames(design) <- make.names(colnames(design))

cat("\nDesign matrix:\n")
print(design)

if (qr(design)$rank != ncol(design)) {
    stop("Design matrix is not full rank.")
}

if (nrow(design) <= ncol(design)) {
    stop("Insufficient residual degrees of freedom.")
}

# ------------------------------------------------------------
# 8. Fit limma model
# ------------------------------------------------------------

fit <- lmFit(
    logCPM,
    design
)

# ------------------------------------------------------------
# 9. Define four biological contrasts
# ------------------------------------------------------------

contrast_matrix <- makeContrasts(

    # BRD4 KO vs Control without IFNg
    BRD4_Effect =
        BRD4KO,

    # IFNg vs No_IFNg in Control cells
    IFNg_Effect =
        IFNgIFNg,

    # IFNg vs No_IFNg in BRD4 KO cells
    IFNg_Response_KO =
        IFNgIFNg + BRD4KO.IFNgIFNg,

    # Difference between IFNg responses
    # in Control and BRD4 KO cells
    BRD4_IFNg_Interaction =
        BRD4KO.IFNgIFNg,

    levels = design
)

cat("\nContrast matrix:\n")
print(contrast_matrix)

# ------------------------------------------------------------
# 10. Apply contrasts and empirical Bayes moderation
# ------------------------------------------------------------

fit2 <- contrasts.fit(
    fit,
    contrast_matrix
)

fit2 <- eBayes(
    fit2,
    trend = TRUE
)

# ------------------------------------------------------------
# 11. Create output directories
# ------------------------------------------------------------

output_directories <- c(
    "BRD4_Effect",
    "IFNg_Effect",
    "IFNg_Response_KO",
    "BRD4_IFNg_Interaction"
)

for (directory in output_directories) {
    dir.create(
        directory,
        showWarnings = FALSE,
        recursive = TRUE
    )
}

# ------------------------------------------------------------
# 12. Function to generate result tables
# ------------------------------------------------------------

create_results <- function(coef_name) {

    results <- topTable(
        fit2,
        coef = coef_name,
        number = Inf,
        sort.by = "P",
        adjust.method = "BH"
    )

    results$GeneID <- rownames(results)

    results$TranscriptID <- annotation$TranscriptID[
        match(
            results$GeneID,
            annotation$GeneID
        )
    ]

    results$Significant <- ifelse(
        !is.na(results$adj.P.Val) &
        results$adj.P.Val < 0.05 &
        abs(results$logFC) >= 1,
        "YES",
        "NO"
    )

    results <- results[
        ,
        c(
            "GeneID",
            "TranscriptID",
            "logFC",
            "AveExpr",
            "t",
            "P.Value",
            "adj.P.Val",
            "B",
            "Significant"
        )
    ]

    results
}

# ------------------------------------------------------------
# 13. Generate results for all four contrasts
# ------------------------------------------------------------

brd4_results <- create_results(
    "BRD4_Effect"
)

ifng_results <- create_results(
    "IFNg_Effect"
)

ko_ifng_results <- create_results(
    "IFNg_Response_KO"
)

interaction_results <- create_results(
    "BRD4_IFNg_Interaction"
)

# ------------------------------------------------------------
# 14. Save complete differential expression tables
# ------------------------------------------------------------

write.csv(
    brd4_results,
    "BRD4_Effect/BRD4_effect_results.csv",
    row.names = FALSE
)

write.csv(
    ifng_results,
    "IFNg_Effect/IFNg_effect_results.csv",
    row.names = FALSE
)

write.csv(
    ko_ifng_results,
    "IFNg_Response_KO/IFNg_response_KO_results.csv",
    row.names = FALSE
)

write.csv(
    interaction_results,
    "BRD4_IFNg_Interaction/BRD4_IFNg_interaction_results.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 15. Save significant gene tables
# ------------------------------------------------------------

write.csv(
    subset(brd4_results, Significant == "YES"),
    "BRD4_Effect/BRD4_effect_significant.csv",
    row.names = FALSE
)

write.csv(
    subset(ifng_results, Significant == "YES"),
    "IFNg_Effect/IFNg_effect_significant.csv",
    row.names = FALSE
)

write.csv(
    subset(ko_ifng_results, Significant == "YES"),
    "IFNg_Response_KO/IFNg_response_KO_significant.csv",
    row.names = FALSE
)

write.csv(
    subset(interaction_results, Significant == "YES"),
    "BRD4_IFNg_Interaction/BRD4_IFNg_interaction_significant.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 16. Calculate significant gene counts
# ------------------------------------------------------------

count_significant <- function(results) {

    sum(
        results$Significant == "YES",
        na.rm = TRUE
    )
}

brd4_sig <- count_significant(brd4_results)
ifng_sig <- count_significant(ifng_results)
ko_ifng_sig <- count_significant(ko_ifng_results)
interaction_sig <- count_significant(interaction_results)

# ------------------------------------------------------------
# 17. Generate differential expression summary
# ------------------------------------------------------------

summary_text <- paste0(

    "GSE248821 — STEP 3: DIFFERENTIAL EXPRESSION\n",
    "==============================================\n\n",

    "Input: Processed RSEM-CPM expression matrix\n",
    "Transformation: log2(CPM + 1)\n",
    "Method: limma-trend\n",
    "Design: ~ BRD4 * IFNg\n",
    "Multiple testing correction: Benjamini-Hochberg\n",
    "Significance: FDR < 0.05 and |logFC| >= 1\n\n",

    "Original genes in STEP 1 output: ",
    genes_before_filtering, "\n",

    "Genes used in differential expression: ",
    genes_after_filtering, "\n",

    "Number of samples: ",
    ncol(logCPM), "\n\n",

    "Experimental groups:\n",
    "Control: 4\n",
    "BRD4 KO: 3\n",
    "Control + IFNg: 3\n",
    "BRD4 KO + IFNg: 3\n\n",

    "Significant genes:\n",

    "BRD4 Effect: ",
    brd4_sig, "\n",

    "IFNg Effect: ",
    ifng_sig, "\n",

    "IFNg Response in BRD4 KO: ",
    ko_ifng_sig, "\n",

    "BRD4 x IFNg Interaction: ",
    interaction_sig, "\n"
)

writeLines(
    summary_text,
    "DE_summary.txt"
)

# ------------------------------------------------------------
# 18. Final status
# ------------------------------------------------------------

cat("\n====================================\n")
cat("STEP 3 COMPLETED SUCCESSFULLY\n")
cat("====================================\n\n")

cat(summary_text)

cat("\nGenerated output directories:\n")

for (directory in output_directories) {
    cat("-", directory, "\n")
}
