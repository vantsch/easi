process MAJIQ_BUILDGFF3 {

    tag "${meta.id ?: 'annotation'}"
    label 'process_single'
    secret 'MAJIQ_LICENSE'

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path('splicegraph.zarr'), emit: splicegraph
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

    /Users/ucl/micromamba/envs/majiq/bin/majiq-build \\
        gff3 \\
        ${gff} \\
        splicegraph.zarr \\
        ${args}
    """

    stub:
    """
    mkdir -p splicegraph.zarr
    echo "Stub splicegraph content" > splicegraph.zarr/dummy_file.txt
    """
}