process MAJIQ_BUILDSJ {

    tag "${meta.id}"
    label 'process_single'

    input:
    tuple(
        val(meta),
        path(bam),
        val(meta_splicegraph),
        path(splicegraph)
    )

    output:
    tuple val(meta), path("build/sj/*.sj"), emit: sj

    tuple val("${task.process}"),
          val('majiq'),
          eval('majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p build/sj

    majiq-build \
        sj \
        --nthreads ${task.cpus} \
        ${bam} \
        ${splicegraph} \
        build/sj/${prefix}.sj \
        --strandness AUTO \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p build/sj
    echo "Stub splice junction content" > build/sj/${prefix}.sj
    """
}