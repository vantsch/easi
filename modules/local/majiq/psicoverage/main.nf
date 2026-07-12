process MAJIQ_PSICOVERAGE {

    tag "${condition}"
    label 'process_single'
    secret 'MAJIQ_LICENSE'

    input:
    tuple val(condition), path(sj), val(meta_splicegraph), path(splicegraph)

    output:
    tuple val(condition), path("psi-coverage/*.psicov"), emit: psi_coverage
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

    mkdir -p psi-coverage

    /Users/ucl/micromamba/envs/majiq/bin/majiq \\
        psi-coverage \\
        --nthreads ${task.cpus} \\
        ${splicegraph} \\
        psi-coverage/${prefix}.psicov \\
        ${sj} \\
        ${args}
    """

    stub:
    """
    mkdir -p psi-coverage
    touch psi-coverage/${condition}.psicov
    """
}