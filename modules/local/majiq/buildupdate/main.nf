process MAJIQ_BUILDUPDATE {

    tag "${meta.id ?: 'majiq_update'}"
    label 'process_low'

    input:
    path sj
    tuple val(meta), path(splicegraph)
    val list_conditions

    output:
    tuple val(meta), path('splicegraph_updated.zarr'), emit: splicegraph
    tuple val("${task.process}"),
          val('majiq'),
          eval('majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def config_rows = list_conditions
        .collect { pair ->
            "printf '%s\\t%s\\n' '${pair[0]}' '${pair[1]}' >> config.tsv"
        }
        .join('\n')

    """
    printf 'sj\\tgroup\\n' > config.tsv
    ${config_rows}

    echo "MAJIQ group configuration:"
    cat config.tsv

    majiq-v3 build \
        ${splicegraph} \
        splicegraph_updated.zarr \
        --groups-tsv config.tsv \
        --nthreads ${task.cpus} \
        ${args}
    """

    stub:
    """
    mkdir -p splicegraph_updated.zarr
    echo "Stub updated splicegraph content" \
        > splicegraph_updated.zarr/dummy_file.txt
    """
}