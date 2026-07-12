process LEAFCUTTER_DS {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/kuechleo/leafcutter_image:2.2.0'

    input:
    tuple val(meta), path(counts), path(groups)

    output:
    tuple val(meta), path("*_cluster_significance.txt"), emit: cluster_significance
    tuple val(meta), path("*_effect_sizes.txt"), emit: effect_sizes
    tuple val(meta), path("*.pdf"), emit: plots, optional: true
    tuple val("${task.process}"), val('leafcutter'),
        eval('Rscript -e \'cat(as.character(packageVersion("leafcutter")))\''),
        emit: versions_leafcutter,
        topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript /leafcutter/scripts/leafcutter_ds.R \
        --output_prefix ${prefix} \
        --num_threads ${task.cpus} \
        ${args} \
        ${counts} \
        ${groups}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_cluster_significance.txt
    touch ${prefix}_effect_sizes.txt
    """
}