process GENE_LEVEL_INTEGRATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'community.wave.seqera.io/library/pip_pyyaml_numpy_openpyxl_pandas:a884a9e713c9a2a8'

    input:
    tuple val(meta),
          path(deseq2_gene_summary),
          path(majiq_gene_summary),
          path(rmats_gene_summary),
          path(leafcutter_gene_summary),
          path(isa_gene_summary)

    path ranking_config

    output:
    tuple val(meta),
          path("${prefix}.gene_integration_short.tsv"),
          emit: short

    tuple val(meta),
          path("${prefix}.gene_integration_long.tsv"),
          emit: long

    tuple val(meta),
          path("${prefix}.gene_ranking_parameters.tsv"),
          emit: ranking_parameters

    tuple val(meta),
          path("${prefix}.gene_ranking_parameters.json"),
          emit: ranking_parameters_json

    tuple val("${task.process}"),
          val('python'),
          eval('python --version 2>&1 | sed "s/Python //"'),
          emit: versions_python,
          topic: versions

    script:
    prefix = task.ext.prefix ?: meta.id

    """
    integrate_gene_results.py \\
        --deseq2 ${deseq2_gene_summary} \\
        --majiq ${majiq_gene_summary} \\
        --rmats ${rmats_gene_summary} \\
        --leafcutter ${leafcutter_gene_summary} \\
        --isa ${isa_gene_summary} \\
        --config ${ranking_config} \\
        --prefix ${prefix} \\
        --outdir .
    """

    stub:
    prefix = task.ext.prefix ?: meta.id

    """
    touch ${prefix}.gene_integration_short.tsv
    touch ${prefix}.gene_integration_long.tsv
    touch ${prefix}.gene_ranking_parameters.tsv
    touch ${prefix}.gene_ranking_parameters.json
    """
}