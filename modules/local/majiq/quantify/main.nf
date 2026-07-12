process MAJIQ_QUANTIFY {

    tag "${condition}"
    label 'process_low'
    secret 'MAJIQ_LICENSE'

    input:
    tuple val(condition), path(psicov), val(meta_splicegraph), path(splicegraph)

    output:
    tuple val(condition), path("quantify/*.tsv"), emit: psi_tsv
    tuple val("${task.process}"),
          val('majiq'),
          eval('/Users/ucl/micromamba/envs/majiq/bin/majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${condition}"

    """
    export MAJIQ_LICENSE_FILE=\$(mktemp)
    printf '%s' "\$MAJIQ_LICENSE" >| "\$MAJIQ_LICENSE_FILE"

    mkdir -p quantify

    /Users/ucl/micromamba/envs/majiq/bin/majiq-v3 \\
        quantify \\
        --nthreads ${task.cpus} \\
        --output-tsv quantify/${prefix}.tsv \\
        ${psicov} \\
        --splicegraph ${splicegraph} \\
        ${args}
    """

    stub:
    """
    mkdir -p quantify
    echo -e "gene_id\\tpsi_value\\nGENE1\\t0.85\\nGENE2\\t0.92" \
        > quantify/${condition}.tsv
    """
}