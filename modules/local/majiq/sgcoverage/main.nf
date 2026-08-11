process MAJIQ_SGCOVERAGE {

    tag "${condition}"
    label 'process_single'

    input:
    tuple val(condition),
          path(sj_files),
          val(meta_splicegraph),
          path(splicegraph)

    output:
    tuple val(condition), path("*.sgc"), emit: sgc_files
    tuple val("${task.process}"),
          val('majiq'),
          eval('majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    majiq-v3 \
        sg-coverage \
        ${splicegraph} \
        ${condition}.sgc \
        ${sj_files} \
        --nthreads ${task.cpus} \
        ${args}
    """

    stub:
    """
    echo "This is a stub .sgc file for ${condition}" \
        > ${condition}.sgc
    """
}