process RSEQC_STRANDNESS {

    tag 'strandness'
    label 'process_single'

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