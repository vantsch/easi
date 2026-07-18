process PFAM_FORMAT_ISA {
    tag "${meta.id}"
    label 'process_single'

     conda "${moduleDir}/environment.yml"

    container 'community.wave.seqera.io/library/pip_openpyxl_pandas:6a8c6a3118605483'

    input:
    tuple val(meta), path(domtblout)
    path pfam_clans

    output:
    tuple val(meta), path("*.pfam_isa.tsv"),
        emit: pfam_isa

    tuple val("${task.process}"),
        val('python'),
        eval("python --version 2>&1 | sed 's/Python //'"),
        topic: versions,
        emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python ${moduleDir}/bin/format_pfam_for_isa.py \
        --domtblout ${domtblout} \
        --clans ${pfam_clans} \
        --output ${prefix}.pfam_isa.tsv \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.pfam_isa.tsv
    """
}