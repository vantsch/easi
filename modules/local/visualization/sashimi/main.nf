process SASHIMI_PLOT {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'community.wave.seqera.io/library/easi-sashimi:9693973230f8a30e'

    input:
    tuple val(meta),
          path(ranked_events),
          val(sample_meta),
          path(bams),
          path(bais)

    path gtf

    output:
    tuple val(meta),
          path("sashimi/sashimi_manifest.tsv"),
          emit: manifest

    tuple val(meta),
          path("sashimi/selected_ranked_events.tsv"),
          emit: selected_events

    tuple val(meta),
          path("sashimi/sashimi_bam_manifest.tsv"),
          emit: bam_manifest

    tuple val(meta),
          path("sashimi/rank_*.pdf"),
          optional: true,
          emit: plots

    tuple val(meta),
          path("sashimi/rank_*.log"),
          optional: true,
          emit: logs

    tuple val("${task.process}"),
          val('python'),
          eval('python --version 2>&1 | sed "s/Python //"'),
          emit: versions_python,
          topic: versions

    tuple val("${task.process}"),
          val('ggsashimi'),
          eval("ggsashimi.py --version 2>&1 | sed 's/ggsashimi v//'"),
          emit: versions_ggsashimi,
          topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def bam_list = (
        bams instanceof List
            ? bams
            : [bams]
    )

    def bai_list = (
        bais instanceof List
            ? bais
            : [bais]
    )

    def metadata_list = (
        sample_meta instanceof List
            ? sample_meta
            : [sample_meta]
    )

    if (
        bam_list.size()
        != bai_list.size()
    ) {
        error """
        SASHIMI_PLOT received different numbers of BAM and BAI files.

        BAM files: ${bam_list.size()}
        BAI files: ${bai_list.size()}
        """.stripIndent()
    }

    if (
        bam_list.size()
        != metadata_list.size()
    ) {
        error """
        SASHIMI_PLOT received different numbers of BAM files and
        sample metadata records.

        BAM files: ${bam_list.size()}
        Metadata records: ${metadata_list.size()}
        """.stripIndent()
    }

    def manifest_rows = (
        0..<bam_list.size()
    ).collect { index ->

        def sample = metadata_list[index]

        def sample_id = (
            sample.id
            ?: sample.sample
            ?: "sample_${index + 1}"
        )

        def condition = (
            sample.condition
            ?: sample.group
            ?: sample.phenotype
            ?: 'UNSPECIFIED'
        )

        def label = (
            sample.label
            ?: sample_id
        )

        [
            sample_id,
            bam_list[index],
            condition,
            bai_list[index],
            label,
        ].join('\t')

    }.join('\n')

    def manifest_content =
        "sample_id\tbam\tcondition\tbai\tlabel\n${manifest_rows}\n"

    def manifest_base64 = manifest_content
        .getBytes('UTF-8')
        .encodeBase64()
        .toString()

    """
    python -c "import base64; open('sashimi_input_bams.tsv', 'wb').write(base64.b64decode('${manifest_base64}'))"

    plot_ranked_sashimi.py \
        --ranked-events ${ranked_events} \
        --bam-manifest sashimi_input_bams.tsv \
        --gtf ${gtf} \
        --outdir sashimi \
        --workdir .sashimi_work \
        --ggsashimi ggsashimi.py \
        --top-n 30 \
        --flank 500 \
        --min-coverage 3 \
        --aggregation mean_j \
        --output-format pdf \
        --height 2.0 \
        --annotation-height 4.5 \
        --strand-mode NONE ${args}
    """

    stub:
    """
    mkdir -p sashimi

    cat > sashimi/sashimi_manifest.tsv <<-'EOF'
    rank	gene_names	chromosome	event_start	event_end	status
    1	STUB	chr1	100	200	SUCCESS
    EOF

    cat > sashimi/selected_ranked_events.tsv <<-'EOF'
    rank	gene_names	chromosome	region_start	region_end
    1	STUB	chr1	100	200
    EOF

    cat > sashimi/sashimi_bam_manifest.tsv <<-'EOF'
    sample_id	bam	condition	bai	label
    STUB	stub.bam	CONDITION	stub.bam.bai	STUB
    EOF

    touch sashimi/rank_001_STUB_chr1_100_200.pdf
    touch sashimi/rank_001_STUB_chr1_100_200.log
    """
}
