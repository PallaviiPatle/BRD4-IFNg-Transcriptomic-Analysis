process DEG_PROCESSING {

    tag "GSE248821 DEG Processing"

    publishDir "${params.output_dir}/05_DEG_Results",
        mode: 'copy',
        overwrite: true

    publishDir "${params.output_dir}/07_Visualizations",
        pattern: "Volcano/*.png",
        mode: 'copy',
        overwrite: true

    input:
    path brd4_results
    path ifng_results
    path ko_ifng_results
    path interaction_results

    output:
    path "DEG_Results", emit: deg_results
    path "Volcano/*.png", emit: volcano_results

    script:
    """
    python ${projectDir}/bin/step4_deg_processing.py \
        "${brd4_results}/BRD4_effect_results.csv" \
        "${ifng_results}/IFNg_effect_results.csv" \
        "${ko_ifng_results}/IFNg_response_KO_results.csv" \
        "${interaction_results}/BRD4_IFNg_interaction_results.csv"
    """
}
