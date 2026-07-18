process LEAFCUTTER_DOWNSTREAM {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/easi-deseq2:1.1'

    input:
    tuple val(meta),
          path(cluster_significance),
          path(effect_sizes),
          path(counts),
          path(groups),
          path(gtf)

    output:
    tuple val(meta), path("*.leafcutter.clusters.annotated.tsv"), emit: clusters_annotated
    tuple val(meta), path("*.leafcutter.clusters.significant.tsv"), emit: clusters_significant
    tuple val(meta), path("*.leafcutter.introns.annotated.tsv"), emit: introns_annotated
    tuple val(meta), path("*.leafcutter.introns.significant.tsv"), emit: introns_significant
    tuple val(meta), path("*.leafcutter.counts.annotated.tsv"), emit: counts_annotated
    tuple val(meta), path("*.leafcutter.sample_groups.tsv"), emit: sample_groups
    tuple val(meta), path("*.leafcutter.enrichment_gene_sets.tsv"), emit: enrichment_gene_sets
    tuple val(meta), path("*.leafcutter.*_ora.tsv"), emit: enrichment_tables
    tuple val(meta), path("*.leafcutter.*_ora.*plot.pdf"), emit: enrichment_plots_pdf, optional: true
    tuple val(meta), path("*.leafcutter.*_ora.*plot.png"), emit: enrichment_plots_png, optional: true
    tuple val(meta), path("*.leafcutter.summary.tsv"), emit: summary
    tuple val(meta), path("*.leafcutter.R_sessionInfo.log"), emit: session_info
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript ${moduleDir}/bin/leafcutter_downstream.R \
        --cluster_file ${cluster_significance} \
        --effect_file ${effect_sizes} \
        --counts_file ${counts} \
        --groups_file ${groups} \
        --gtf_file ${gtf} \
        --output_prefix ${prefix} \
        ${args}
    """ 

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.leafcutter.clusters.annotated.tsv
    touch ${prefix}.leafcutter.clusters.significant.tsv
    touch ${prefix}.leafcutter.introns.annotated.tsv
    touch ${prefix}.leafcutter.introns.significant.tsv
    touch ${prefix}.leafcutter.counts.annotated.tsv
    touch ${prefix}.leafcutter.sample_groups.tsv
    touch ${prefix}.leafcutter.summary.tsv
    touch ${prefix}.leafcutter.R_sessionInfo.log
    touch ${prefix}.leafcutter.enrichment_gene_sets.tsv
    touch ${prefix}.leafcutter.significant.hallmark_ora.tsv
    touch ${prefix}.leafcutter.significant_effect_filtered.hallmark_ora.tsv
    touch ${prefix}.leafcutter.significant.go_bp_ora.tsv
    touch ${prefix}.leafcutter.significant.go_mf_ora.tsv
    touch ${prefix}.leafcutter.significant.go_cc_ora.tsv
    touch ${prefix}.leafcutter.significant_effect_filtered.go_bp_ora.tsv
    touch ${prefix}.leafcutter.significant_effect_filtered.go_mf_ora.tsv
    touch ${prefix}.leafcutter.significant_effect_filtered.go_cc_ora.tsv
    touch ${prefix}.leafcutter.significant.kegg_ora.tsv
    touch ${prefix}.leafcutter.significant_effect_filtered.kegg_ora.tsv
    """
}