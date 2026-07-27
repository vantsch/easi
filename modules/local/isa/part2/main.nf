process ISA_PART2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bioconductor-isoformswitchanalyzer==2.2.0"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0'}"

    input:
    tuple val(meta), path(switchlist), path(pfam_isa)
    path gtf

    output:
    tuple val(meta), path("switchlist_part2.rds"),
        emit: switchlist

    tuple val(meta), path("isa_part2_domain_analysis.tsv"),
        emit: domain_analysis

    tuple val(meta), path("isa_part2_isoform_features.tsv"),
        emit: isoform_features

    tuple val(meta), path("switch_plots"),
        emit: switch_plots

    tuple val(meta), path("common_switch_consequences.pdf"),
        emit: consequence_summary

    tuple val(meta), path("isa_part2.gene_summary.tsv"),
        emit: gene_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    isa_part2.R \
        ${switchlist} \
        ${pfam_isa} \
        ${gtf} \
        ${args}
    """

    stub:
    """
    touch switchlist_part2.rds
    touch isa_part2_domain_analysis.tsv
    touch isa_part2_isoform_features.tsv
    touch isa_part2.gene_summary.tsv

    mkdir -p switch_plots
    touch switch_plots/switchplot_stub.pdf

    touch common_switch_consequences.pdf
    """
}