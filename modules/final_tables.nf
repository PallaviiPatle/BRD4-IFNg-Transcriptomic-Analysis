process FINAL_TABLES {

    tag "GSE248821 Final Tables"

    publishDir "${params.output_dir}",
        mode: 'copy',
        overwrite: true

    input:
    path deg_results
    path pathway_results

    output:
    path "08_Final_Tables"

    script:
    """
    mkdir -p 05_DEG_Results
    mkdir -p 06_Pathway_Analysis

    cp -r "${deg_results}" 05_DEG_Results/DEG_Results
    cp -r "${pathway_results}" 06_Pathway_Analysis/Pathway_Analysis

    python ${projectDir}/bin/step7_final_tables.py
    """
}
