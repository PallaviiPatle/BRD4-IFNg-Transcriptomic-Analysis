process PATHWAY_ANALYSIS {

    tag "GSE248821 Pathway Analysis"

    publishDir "${params.output_dir}/06_Pathway_Analysis",
        mode: 'copy',
        overwrite: true

    input:
    path brd4_results
    path ifng_results
    path ko_ifng_results
    path interaction_results

    output:
    path "Pathway_Analysis", emit: pathway_results

    script:
    """
    echo "Running actual pathway analysis v2"

    Rscript ${projectDir}/bin/step5_pathway_analysis.R \
        "${brd4_results}/BRD4_effect_results.csv" \
        "${ifng_results}/IFNg_effect_results.csv" \
        "${ko_ifng_results}/IFNg_response_KO_results.csv" \
        "${interaction_results}/BRD4_IFNg_interaction_results.csv"
    """
}
