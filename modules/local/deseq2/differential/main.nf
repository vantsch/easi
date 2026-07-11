process DESEQ2_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/library/easi-deseq2:1.1'
    
    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target), val(formula), val(comparison)
    tuple val(meta2), path(samplesheet), path(counts)
    tuple val(meta3), path(control_genes_file)
    tuple val(meta4), path(transcript_lengths_file)
    path report_template, stageAs: "report_template.Rmd"

    output:
    tuple val(meta), path("*.deseq2.results.tsv")              , emit: results
    tuple val(meta), path("*.deseq2.dispersion.png")           , emit: dispersion_plot_png
    tuple val(meta), path("*.deseq2.dispersion.pdf")           , emit: dispersion_plot_pdf
    tuple val(meta), path("*.dds.rld.rds")                     , emit: rdata
    tuple val(meta), path("*.deseq2.sizefactors.tsv")          , emit: size_factors
    tuple val(meta), path("*.normalised_counts.tsv")           , emit: normalised_counts
    tuple val(meta), path("*.rlog.tsv")                        , optional: true, emit: rlog_counts
    tuple val(meta), path("*.vst.tsv")                         , optional: true, emit: vst_counts
    tuple val(meta), path("*.deseq2.PCA.labelled.pdf")         , emit: pca_pdf
    tuple val(meta), path("*.deseq2.PCA.labelled.png")         , emit: pca_png
    tuple val(meta), path("*.deseq2.sig_genes.log2FC0.5.tsv")  , emit: sig_genes_lfc_0_5
    tuple val(meta), path("*.deseq2.sig_genes.log2FC1.tsv")    , emit: sig_genes_lfc_1
    tuple val(meta), path("*.deseq2.sig_genes.log2FC1.5.tsv")  , emit: sig_genes_lfc_1_5
    tuple val(meta), path("*.deseq2.volcano.pdf")              , emit: volcano_pdf
    tuple val(meta), path("*.deseq2.volcano.png")              , emit: volcano_png
    tuple val(meta), path("*.deseq2.MA.pdf")                   , emit: ma_pdf
    tuple val(meta), path("*.deseq2.MA.png")                   , emit: ma_png
    tuple val(meta), path("*.hallmark_gsea.tsv")               , emit: hallmark_gsea
    tuple val(meta), path("*.hallmark_gsea.dotplot.pdf")       , optional: true, emit: hallmark_gsea_dotplot_pdf
    tuple val(meta), path("*.hallmark_gsea.dotplot.png")       , optional: true, emit: hallmark_gsea_dotplot_png
    tuple val(meta), path("*.hallmark_gsea.barplot.pdf")       , optional: true, emit: hallmark_gsea_barplot_pdf
    tuple val(meta), path("*.hallmark_gsea.barplot.png")       , optional: true, emit: hallmark_gsea_barplot_png
    tuple val(meta), path("*.go_*_ora.tsv")                    , emit: go_ora
    tuple val(meta), path("*.go_*_ora.dotplot.pdf")            , optional: true, emit: go_ora_dotplot_pdf
    tuple val(meta), path("*.go_*_ora.dotplot.png")            , optional: true, emit: go_ora_dotplot_png
    tuple val(meta), path("*.go_*_ora.barplot.pdf")            , optional: true, emit: go_ora_barplot_pdf
    tuple val(meta), path("*.go_*_ora.barplot.png")            , optional: true, emit: go_ora_barplot_png
    tuple val(meta), path("*.deseq2.report.html")              , emit: report_html
    tuple val(meta), path("*.deseq2.model.txt")                , emit: model
    tuple val(meta), path("*.R_sessionInfo.log")               , emit: session_info
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq2_differential.R'

    stub:
    """
    touch ${meta.id}.deseq2.results.tsv
    touch ${meta.id}.deseq2.dispersion.png
    touch ${meta.id}.deseq2.dispersion.pdf
    touch ${meta.id}.dds.rld.rds
    touch ${meta.id}.deseq2.sizefactors.tsv
    touch ${meta.id}.normalised_counts.tsv
    touch ${meta.id}.rlog.tsv
    touch ${meta.id}.deseq2.model.txt
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
