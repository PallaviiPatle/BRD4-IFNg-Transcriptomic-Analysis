process DIFFERENTIAL_EXPRESSION {

    tag "GSE248821 Differential Expression"

    publishDir "${params.output_dir}/04_Differential_Expression",
        mode: 'copy',
        overwrite: true

    input:
    path expression_file
    path metadata_file
    path annotation_file

    output:
    path "BRD4_Effect", emit: brd4_results
    path "IFNg_Effect", emit: ifng_results
    path "IFNg_Response_KO", emit: ko_ifng_results
    path "BRD4_IFNg_Interaction", emit: interaction_results
    path "DE_summary.txt"

    script:
    """
    Rscript ${projectDir}/bin/step3_differential_expression.R
    """
}
