process MAJIQ_SUMMARIZE {

    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/pip_openpyxl_pandas:6a8c6a3118605483'

    input:
    tuple val(meta), path(tsv_files)
    path report_script

    output:
    path 'MAJIQ_summary.xlsx',
        emit: workbook

    tuple val(meta), path('MAJIQ_gene_summary.tsv'),
        emit: gene_summary

    script:
    """
    python3 ${report_script} .
    """
}