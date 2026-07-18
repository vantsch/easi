process PFAM_SCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_2' :
        'quay.io/biocontainers/hmmer:3.4--hdbdd923_2'}"

    input:
    tuple val(meta), path(aa_fasta)

    tuple path(pfam_hmm),
          path(pfam_h3f),
          path(pfam_h3i),
          path(pfam_h3m),
          path(pfam_h3p)

    output:
    tuple val(meta), path("*.pfam.domtblout"),
        emit: domtblout

    tuple val(meta), path("*.pfam.tblout"),
        emit: tblout

    tuple val("${task.process}"),
        val('hmmer'),
        eval("hmmscan -h | grep -m 1 '^# HMMER' | sed 's/^# HMMER //'"),
        topic: versions,
        emit: versions_hmmer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    hmmscan \
        --cpu ${task.cpus} \
        --cut_ga \
        --noali \
        --domtblout ${prefix}.pfam.domtblout \
        --tblout ${prefix}.pfam.tblout \
        ${args} \
        ${pfam_hmm} \
        ${aa_fasta}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.pfam.domtblout
    touch ${prefix}.pfam.tblout
    """
}