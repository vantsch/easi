process FEATURECOUNTS_MERGE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    
    container 'community.wave.seqera.io/library/pip_openpyxl_pandas:6a8c6a3118605483'

    input:
    path counts
    path merge_script

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
    python ${merge_script} \\
        --counts ${counts.join(' ')} \\
        --out featureCounts.matrix.tsv
    """

    stub:
    """
    touch featureCounts.matrix.tsv
    """
}