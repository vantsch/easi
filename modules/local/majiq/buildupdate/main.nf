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
        .collect { pair -> "${pair[0]}\t${pair[1]}" }
        .join('\n')

    """
    cat > config.tsv <<'EOF'
    sj	group
    ${config_rows}
    EOF

    echo "MAJIQ group configuration:"
    cat config.tsv

    majiq-build \
        update \
        --nthreads ${task.cpus} \
        --groups-tsv config.tsv \
        ${splicegraph} \
        splicegraph_updated.zarr \
        ${args}
    """

    stub:
    """
    mkdir -p splicegraph_updated.zarr
    echo "Stub updated splicegraph content" \
        > splicegraph_updated.zarr/dummy_file.txt
    """
}