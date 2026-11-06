
nextflow.enable.dsl=2

// ==========================================
// MODULE IMPORTS
// ==========================================

include { STEP1_DATA_PREPARATION } from './modules/data_preparation'
include { SAMPLE_QC } from './modules/sample_qc'
include { DIFFERENTIAL_EXPRESSION } from './modules/differential_expression'
include { DEG_PROCESSING } from './modules/deg_processing'
include { PATHWAY_ANALYSIS } from './modules/pathway_analysis'
include { HEATMAP_ANALYSIS } from './modules/heatmap_analysis'
include { FINAL_TABLES } from './modules/final_tables'


// ==========================================
// MAIN WORKFLOW
// ==========================================

workflow {

    // Input Excel file
    input_excel = file(params.input_file)

    // STEP 1: Data Preparation
    step1 = STEP1_DATA_PREPARATION(input_excel)

    // STEP 2: Sample Quality Control
    SAMPLE_QC(
        step1.expression_CPM,
        step1.sample_metadata
    )

    // STEP 3: Differential Expression Analysis
    step3 = DIFFERENTIAL_EXPRESSION(
        step1.expression_CPM,
        step1.sample_metadata,
        step1.gene_annotation
    )

    // STEP 4: DEG Processing
    deg = DEG_PROCESSING(
        step3.brd4_results,
        step3.ifng_results,
        step3.ko_ifng_results,
        step3.interaction_results
    )

    // STEP 5: Pathway Analysis
    pathways = PATHWAY_ANALYSIS(
        step3.brd4_results,
        step3.ifng_results,
        step3.ko_ifng_results,
        step3.interaction_results
    )

    // STEP 6: Interaction Heatmap
    HEATMAP_ANALYSIS(
        step1.expression_CPM,
        step1.gene_annotation,
        step3.interaction_results
    )

    // STEP 7: Final Tables
    FINAL_TABLES(
        deg.deg_results,
        pathways.pathway_results
    )
}
