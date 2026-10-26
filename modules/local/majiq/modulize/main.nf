process MAJIQ_MODULIZE {

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(voila_files), path(splicegraph)

    output:
    tuple val(meta), path("${prefix}"),       emit: modulize
    tuple val(meta), path("${prefix}/*.tsv"), emit: tsv_files

    tuple val("${task.process}"),
          val('voila'),
          eval('voila --version | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+"'),
          emit: versions_voila,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'modulize'

    """
    mkdir -p ${prefix}

    voila \
        modulize \
        ${splicegraph} \
        ${voila_files} \
        --logger ${prefix}/voila.log \
        --directory ${prefix} \
        --nproc ${task.cpus} \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'modulize'

    """
    mkdir -p ${prefix}
    touch ${prefix}/dummy_output.txt
    echo -e "module_id\\tstatus\\nMODULE1\\tstub" \
        > ${prefix}/dummy_output.tsv
    echo "This is a stub log file" \
        > ${prefix}/voila.log
    """
}