process MAJIQ_PSICOVERAGE {

    tag "${condition}"
    label 'process_single'

    input:
    tuple val(condition), path(sj), val(meta_splicegraph), path(splicegraph)

    output:
    tuple val(condition), path("psi-coverage/*.psicov"), emit: psi_coverage
    tuple val("${task.process}"),
          val('majiq'),
          eval('majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${condition}"

    """
    mkdir -p psi-coverage

    majiq \
        psi-coverage \
        --nthreads ${task.cpus} \
        ${splicegraph} \
        psi-coverage/${prefix}.psicov \
        ${sj} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${condition}"

    """
    mkdir -p psi-coverage
    touch psi-coverage/${prefix}.psicov
    """
}