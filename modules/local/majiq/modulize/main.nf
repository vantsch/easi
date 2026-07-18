process MAJIQ_MODULIZE {

    label 'process_high'
    secret 'MAJIQ_LICENSE'

    input:
    tuple path(voila_files), val(meta_splicegraph), path(splicegraph)

    output:
    path "${prefix}", emit: modulize
    path "${prefix}/*.tsv", emit: tsv_files

    tuple val("${task.process}"), val('voila'),
        eval('/Users/ucl/micromamba/envs/majiq/bin/voila --version | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+"'),
        emit: versions_voila,
        topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'modulize'

    """
    export MAJIQ_LICENSE_FILE=\$(mktemp)
    printf '%s' "\$MAJIQ_LICENSE" >| "\$MAJIQ_LICENSE_FILE"

    mkdir -p ${prefix}

    /Users/ucl/micromamba/envs/majiq/bin/voila \\
        modulize \\
        ${splicegraph} \\
        ${voila_files} \\
        --logger ${prefix}/voila.log \\
        --directory ${prefix} \\
        --nproc ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'modulize'

    """
    mkdir -p ${prefix}
    touch ${prefix}/dummy_output.txt
    echo "This is a stub log file" > ${prefix}/voila.log
    """
}