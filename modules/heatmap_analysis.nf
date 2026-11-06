process HEATMAP_ANALYSIS {

    tag "GSE248821 Interaction Heatmap"

    publishDir "${params.output_dir}",
        mode: 'copy',
        overwrite: true

    input:
    path expression_file
    path annotation_file
    path interaction_file

    output:
    path "07_Visualizations/Heatmaps"

    script:
    """
    mkdir -p 02_Data_Preparation
    mkdir -p 04_Differential_Expression/BRD4_IFNg_Interaction

    cp "${expression_file}" \
        02_Data_Preparation/expression_CPM.csv

    cp "${annotation_file}" \
        02_Data_Preparation/gene_annotation.csv

    cp "${interaction_file}/BRD4_IFNg_interaction_results.csv" \
        04_Differential_Expression/BRD4_IFNg_Interaction/BRD4_IFNg_interaction_results.csv

    python ${projectDir}/bin/step6_heatmap.py
    """
}
