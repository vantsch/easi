process LEAFCUTTER_CLUSTER {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/vantsch/easi-leafcutter:1.2'

    input:
    tuple val(meta), path(junc_files), val(groups)

    output:
    tuple val(meta), path("*_perind.counts.gz"), emit: counts
    tuple val(meta), path("*_perind_numers.counts.gz"), emit: numerators
    tuple val(meta), path("juncfiles.txt"), emit: junction_list
    tuple val(meta), path("groups_file.txt"), emit: groups_file
    tuple val("${task.process}"), val('leafcutter'),
        eval('printf "%s\n" "$LEAFCUTTER_COMMIT"'),
        emit: versions_leafcutter,
        topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def junc_args = junc_files
        .collect { "'${it}'" }
        .join(' ')

    def group_args = groups
        .collect { sample, condition ->
            "'${sample}\\t${condition}'"
        }
        .join(' ')

    """
    printf '%s\\n' ${junc_args} > juncfiles.txt
    printf '%b\\n' ${group_args} > groups_file.txt

    leafcutter_cluster_regtools.py \
        -j juncfiles.txt \
        -o ${prefix} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_perind.counts.gz
    touch ${prefix}_perind_numers.counts.gz
    touch juncfiles.txt
    touch groups_file.txt
    """
}