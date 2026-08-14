process RSEQC_STRANDNESS {

    tag 'strandness'
    label 'process_single'

    container 'community.wave.seqera.io/library/pip_openpyxl_pandas:6a8c6a3118605483'

    input:
    path reports

    output:
    stdout emit: strandness

    script:
    """
    python ${moduleDir}/bin/parse_rseqc_strandness.py \
        ${reports.join(' ')}
    """

    stub:
    """
    echo "REVERSE"
    """
}