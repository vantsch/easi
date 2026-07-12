process MAJIQ_SGCOVERAGE {

    tag "${condition}"
    label 'process_single'
    secret 'MAJIQ_LICENSE'

    input:
    tuple val(condition),
          path(sj_files),
          val(meta_splicegraph),
          path(splicegraph)

    output:
    tuple val(condition), path("*.sgc"), emit: sgc_files
    tuple val("${task.process}"),
          val('majiq'),
          eval('/Users/ucl/micromamba/envs/majiq/bin/majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    export MAJIQ_LICENSE_FILE=\$(mktemp)
    printf '%s' "\$MAJIQ_LICENSE" >| "\$MAJIQ_LICENSE_FILE"

    /Users/ucl/micromamba/envs/majiq/bin/majiq-v3 \\
        sg-coverage \\
        ${splicegraph} \\
        ${condition}.sgc \\
        ${sj_files} \\
        --nthreads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    echo "This is a stub .sgc file for ${condition}" > ${condition}.sgc
    """
}