process MAJIQ_MODULIZE {

    label 'process_high'
    secret 'MAJIQ_LICENSE'

    input:
    tuple path(voila_files), val(meta_splicegraph), path(splicegraph)

    output:
    path("modulize"), emit: modulize_files
    tuple val("${task.process}"),
          val('voila'),
          eval('/Users/ucl/micromamba/envs/majiq/bin/voila --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    export MAJIQ_LICENSE_FILE=\$(mktemp)
    printf '%s' "\$MAJIQ_LICENSE" >| "\$MAJIQ_LICENSE_FILE"

    mkdir -p modulize

    /Users/ucl/micromamba/envs/majiq/bin/voila \\
        modulize \\
        ${splicegraph} \\
        ${voila_files} \\
        --logger module.logger.txt \\
        --directory modulize \\
        --nproc ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mkdir -p modulize
    touch modulize/dummy_output.txt
    echo "This is a stub log file" > modulize/module.logger.txt
    """
}