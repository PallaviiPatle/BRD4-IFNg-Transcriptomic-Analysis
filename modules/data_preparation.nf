process STEP1_DATA_PREPARATION {

    tag "GSE248821 Data Preparation"

    input:
    path input_excel

    output:

    path "gene_annotation.csv", emit: gene_annotation
    path "expression_CPM.csv", emit: expression_CPM
    path "sample_metadata.csv", emit: sample_metadata
    path "step1_summary.txt", emit: summary

    script:
    """
    python ${projectDir}/bin/step1_data_preparation.py "${input_excel}"
    """
}
