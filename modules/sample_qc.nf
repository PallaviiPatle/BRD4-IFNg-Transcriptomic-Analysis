process SAMPLE_QC {

    tag "GSE248821 Sample QC"

    publishDir "${params.output_dir}/03_Sample_QC",
        mode: 'copy',
        overwrite: true

    input:
    path expression_file
    path metadata_file

    output:
    path "PCA"
    path "Correlation"
    path "Clustering"
    path "Expression_Distribution"
    path "QC_summary.txt"

    script:
    """
    python ${projectDir}/bin/step2_sample_qc.py
    """
}