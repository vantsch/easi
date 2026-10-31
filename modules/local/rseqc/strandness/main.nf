process RSEQC_STRANDNESS {

    tag "${meta.id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/pip_openpyxl_pandas:6a8c6a3118605483'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), stdout, emit: strandness

    script:
    """
    python ${moduleDir}/bin/parse_rseqc_strandness.py \
        ${report}
    """

    stub:
    """
    echo "REVERSE"
    """
}