process ISA_PART1 {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bioconductor-isoformswitchanalyzer==2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0' }"

    input:
    tuple val(meta), path(salmon_output)
    path gtf
    path transcript_sequences
    path samplesheet
    path contrastsheet
    val alpha
    val dIF

    output:
    tuple val(meta), path("switchlist_part1.rds"),
        emit: switchlist

    tuple val(meta), path("isa_part1_switches.tsv"),
        emit: switches

    tuple val(meta), path("isa_part1_isoform_features.tsv"),
        emit: isoform_features

    tuple val(meta), path("isa_part1_orf_analysis.tsv"),
        emit: orf_analysis,
        optional: true

    tuple val(meta), path("results/isoformSwitchAnalyzeR_isoform_AA.fasta"),
        emit: aa_fasta

    tuple val(meta), path("results/isoformSwitchAnalyzeR_isoform_nt.fasta"),
        emit: nt_fasta

    tuple val(meta), path("results"),
        emit: sequences

    tuple val(meta), path("isa_part1_sequence_files.txt"),
        emit: sequence_manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # ISA contrast orientation fix v1
    mkdir -p results

    Rscript ${moduleDir}/bin/isa_part1.R \
        ${gtf} \
        ${transcript_sequences} \
        ${samplesheet} \
        ${contrastsheet} \
        ${alpha} \
        ${dIF} \
        ${args}
    """
}