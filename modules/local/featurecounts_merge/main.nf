process FEATURECOUNTS_MERGE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    path counts

    output:
    path "featureCounts.matrix.tsv", emit: matrix
    tuple val("${task.process}"), val('python'),
          eval('python --version | sed "s/Python //"'),
          emit: versions_python,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${moduleDir}/bin/merge_featurecounts.py \\
        --counts ${counts.join(' ')} \\
        --out featureCounts.matrix.tsv
    """

    stub:
    """
    touch featureCounts.matrix.tsv
    """
}