process SPLICE_REGION_INTEGRATION {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'community.wave.seqera.io/library/pip_pyyaml_numpy_openpyxl_pandas:a884a9e713c9a2a8'

    input:
    tuple val(meta),
          path(majiq_junctions),
          path(rmats_jc_tables),
          path(leafcutter_introns)

    output:
    tuple val(meta),
          path("splice_region_integration/${meta.id}.significant_region_overlaps.tsv"),
          emit: significant_region_overlaps

    tuple val(meta),
          path("splice_region_integration/${meta.id}.splice_region_integration.xlsx"),
          emit: workbook

    tuple val(meta),
          path("splice_region_integration/${meta.id}.ranked_splicing_events.tsv"),
          emit: ranked_splicing_events

    tuple val("${task.process}"),
          val('python'),
          eval('python --version 2>&1 | sed "s/Python //"'),
          emit: versions_python,
          topic: versions

    script:
    def args = task.ext.args ?: ''

    def rmats_jc_inputs = (
        rmats_jc_tables instanceof List
            ? rmats_jc_tables.join(' ')
            : rmats_jc_tables
    )

    """
    mkdir -p splice_region_integration

    integrate_splice_regions.py \\
        --majiq-junctions ${majiq_junctions} \\
        --rmats-jc ${rmats_jc_inputs} \\
        --leafcutter-introns ${leafcutter_introns} \\
        --prefix ${meta.id} \\
        --outdir splice_region_integration \\
        ${args}
    """

    stub:
    """
    mkdir -p splice_region_integration

    touch splice_region_integration/${meta.id}.significant_region_overlaps.tsv
    touch splice_region_integration/${meta.id}.splice_region_integration.xlsx
    touch splice_region_integration/${meta.id}.ranked_splicing_events.tsv
    """
}