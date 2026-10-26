process MAJIQ_DELTAPSI {

    tag "${contrast}"
    label 'process_medium'

    input:
    tuple val(contrast),
          val(treatment),
          val(control),
          path(treatment_files),
          path(control_files),
          path(splicegraph)

    output:
    path("deltapsi/${contrast}.deltapsi.tsv"),        emit: deltapsi_tsv
    tuple val(contrast), path("deltapsi/${contrast}.dpsicov"), emit: dpsicov
    path("deltapsi/${contrast}.deltapsi.logger.txt"), emit: logger
    tuple val("${task.process}"),
          val('majiq'),
          eval('majiq --version'),
          emit: versions_majiq,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def changing_threshold = task.ext.changing_threshold ?: 0.1
    def args = task.ext.args ?: ''

    """
    mkdir -p deltapsi

    majiq \
        deltapsi \
        --nthreads ${task.cpus} \
        --splicegraph ${splicegraph} \
        -psi1 ${control_files} \
        -psi2 ${treatment_files} \
        --names ${control} ${treatment} \
        --changing-threshold ${changing_threshold} \
        --output-tsv deltapsi/${contrast}.deltapsi.tsv \
        --logger deltapsi/${contrast}.deltapsi.logger.txt \
        --output-voila deltapsi/${contrast}.dpsicov \
        --debug \
        ${args}
    """

    stub:
    """
    mkdir -p deltapsi
    echo -e "gene_id\\tdeltapsi_value\\nGENE1\\t0.5\\nGENE2\\t-0.3" \
        > deltapsi/${contrast}.deltapsi.tsv
    touch deltapsi/${contrast}.dpsicov
    echo "This is a stub log file for ${contrast}" \
        > deltapsi/${contrast}.deltapsi.logger.txt
    """
}