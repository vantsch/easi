process MAJIQ_QUANTIFY {

    tag "${condition}"
    label 'process_low'

    input:
    tuple val(condition), path(psicov), val(meta_splicegraph), path(splicegraph)

    output:
    tuple val(condition), path("quantify/*.tsv"), emit: psi_tsv
    tuple val("${task.process}"),
          val('majiq'),
          eval('majiq-v3 --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${condition}"

    """
    mkdir -p quantify

    majiq-v3 \
        quantify \
        --nthreads ${task.cpus} \
        --output-tsv quantify/${prefix}.tsv \
        ${psicov} \
        --splicegraph ${splicegraph} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${condition}"

    """
    mkdir -p quantify
    echo -e "gene_id\\tpsi_value\\nGENE1\\t0.85\\nGENE2\\t0.92" \
        > quantify/${prefix}.tsv
    """
}