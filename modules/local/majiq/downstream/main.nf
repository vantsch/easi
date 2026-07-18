process MAJIQ_SUMMARIZE {

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/pip_openpyxl_pandas:6a8c6a3118605483'

    input:
    path tsv_files
    path report_script

    output:
    path 'MAJIQ_summary.xlsx', emit: workbook

    script:
    """
    python3 ${report_script} .
    """
}