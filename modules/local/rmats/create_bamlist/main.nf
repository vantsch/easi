process CREATE_BAMLIST {
    label 'process_single'

    conda "conda-forge::sed=4.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'biocontainers/sed:4.7.0' }"

    input:
    tuple val(contrast), val(cond1), path(bam1), val(cond2),path(bam2)

    output:
    tuple val(contrast), path("${cond1}_bamlist.txt"), emit: bam_list1
    tuple val(contrast), path("${cond2}_bamlist.txt"), emit: bam_list2, optional: true
    tuple val("${task.process}"), val('sed'), eval('sed --version 2>&1 | head -n1 | sed "s/.* //"'), topic: versions, emit: versions_sed

    when:
    task.ext.when == null || task.ext.when

    script:
    def b1_list = bam1 instanceof List ? bam1.join(' ') : bam1
    def b2_list = bam2 instanceof List ? bam2.join(' ') : bam2
    def bam2_cmd = (cond2 && b2_list) ? "echo ${b2_list} | sed 's: :,:g' > ${cond2}_bamlist.txt" : ''
    """
    echo ${b1_list} | sed 's: :,:g' > ${cond1}_bamlist.txt
    ${bam2_cmd}
    """
}
