process RMATS_SUMMARIZE {
    tag { contrast?.id ?: (cond2 ? "${cond1}_vs_${cond2}" : "${cond1}") }
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/4697fb2978d58ad7209a2d7b267bc02cf56cb8c256d787767b482c97346330d6/data'
        : 'community.wave.seqera.io/library/r-pairadise_rmats:4a32921359285283'}"

    input:
    path gtf
    tuple val(contrast),
          val(cond1),
          val(meta1),
          path(bam1),
          path(bam1_text),
          val(cond2),
          val(meta2),
          path(bam2),
          path(bam2_text),
          path(rmats_temp)

    val rmats_read_len
    val rmats_splice_diff_cutoff
    val rmats_novel_splice_site
    val rmats_min_intron_len
    val rmats_max_exon_len
    val rmats_paired_stats

    output:
    path "${output_dir}/${rmats_output_dir}/*",
        emit: rmats_post

    tuple val(contrast),
        path("${output_dir}/${rmats_output_dir}/RMATS_gene_summary.tsv"),
        emit: gene_summary
    
    tuple val(contrast),
        path("${output_dir}/${rmats_output_dir}/*.MATS.JC.annotated.tsv"),
        emit: jc_annotated

    path "${output_dir}/${logfile}",
        emit: log

    tuple val("${task.process}"),
        val('rmats'),
        eval('rmats.py --version 2>&1 | sed -e "s/v//g"'),
        topic: versions,
        emit: versions_rmats

    when:
    task.ext.when == null || task.ext.when

    script:

    def contrast_name = contrast?.id ?: (
        cond2 ? "${cond1}_vs_${cond2}" : "${cond1}"
    )

    output_dir = contrast_name
    rmats_output_dir = rmats_paired_stats
        ? 'rmats_post_paired'
        : 'rmats_post'

    /*
     * Only meta1 is required because both groups should have the same
     * read type and strandedness.
     */
    def meta = meta1 instanceof List ? meta1[0] : meta1
    def args = task.ext.args ?: ''

    /*
     * Determine whether the input is single-end or paired-end.
     */
    def read_type = meta.single_end ? 'single' : 'paired'

    /*
     * Convert EASI/nf-core strandedness terminology into rMATS libType.
     */
    def strandedness = 'fr-unstranded'

    if (meta.strandedness == 'forward') {
        strandedness = 'fr-secondstrand'
    }
    else if (meta.strandedness == 'reverse') {
        strandedness = 'fr-firststrand'
    }

    /*
     * Optional novel splice-site handling.
     */
    def novel_splice_sites = rmats_novel_splice_site
        ? '--novelSS'
        : ''

    def min_intron_len = ''
    def max_exon_len = ''

    if (rmats_novel_splice_site) {
        min_intron_len = rmats_min_intron_len
            ? "--mil ${rmats_min_intron_len}"
            : '--mil 50'

        max_exon_len = rmats_max_exon_len
            ? "--mel ${rmats_max_exon_len}"
            : '--mel 500'
    }

    /*
     * cond1 is the target/treatment condition and is passed as b1.
     * cond2 is the reference/control condition and is passed as b2.
     *
     * Therefore:
     * IncLevelDifference = PSI_target - PSI_reference
     */
    def b1 = bam1_text
        ? "--b1 ${bam1_text}"
        : ''

    def b2 = ''
    def stat_flag = ''
    def paired_stats = ''

    if (cond2 && bam2_text) {
        b2 = "--b2 ${bam2_text}"

        paired_stats = rmats_paired_stats
            ? '--paired-stats'
            : ''
    }
    else {
        stat_flag = '--statoff'
    }

    logfile = rmats_paired_stats
        ? 'rmats_post_paired.log'
        : 'rmats_post.log'

    """
    mkdir -p ${output_dir}/${rmats_output_dir}

    rmats.py \\
        ${args} \\
        --gtf ${gtf} \\
        ${b1} \\
        ${b2} \\
        --od ${output_dir}/${rmats_output_dir} \\
        --tmp ./ \\
        -t ${read_type} \\
        --libType ${strandedness} \\
        --readLength ${rmats_read_len} \\
        --variable-read-length \\
        --nthread ${task.cpus} \\
        --tstat ${task.cpus} \\
        --cstat ${rmats_splice_diff_cutoff} \\
        --task post \\
        ${paired_stats} \\
        ${novel_splice_sites} \\
        ${min_intron_len} \\
        ${max_exon_len} \\
        ${stat_flag} \\
        --allow-clipping \\
        1> ${output_dir}/${logfile}

    python - <<'PYTHON'
    import csv
    from pathlib import Path
    import math
    from collections import defaultdict


    output_dir = Path("${output_dir}/${rmats_output_dir}")

    contrast = "${contrast_name}"
    target = "${cond1}"
    reference = "${cond2}"


    column_renames = {
        "IJC_SAMPLE_1": "IJC_target",
        "SJC_SAMPLE_1": "SJC_target",
        "IJC_SAMPLE_2": "IJC_reference",
        "SJC_SAMPLE_2": "SJC_reference",
        "IncLevel1": "PSI_target",
        "IncLevel2": "PSI_reference",
        "IncLevelDifference": "dPSI_target_minus_reference",
    }


    for input_file in sorted(output_dir.glob("*.MATS.*.txt")):

        output_file = input_file.with_name(
            input_file.name.removesuffix(".txt")
            + ".annotated.tsv"
        )

        with input_file.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as source:

            reader = csv.reader(
                source,
                delimiter="\\t",
            )

            try:
                original_header = next(reader)
            except StopIteration:
                continue

            # rMATS combined files can contain more than one column named ID.
            # Rename them deterministically so the annotated table has unique
            # column names.
            id_number = 0
            annotated_header = []

            for column in original_header:
                if column == "ID":
                    id_number += 1

                    if id_number == 1:
                        annotated_header.append("event_id")
                    else:
                        annotated_header.append(
                            f"event_id_{id_number}"
                        )
                else:
                    annotated_header.append(
                        column_renames.get(column, column)
                    )

            with output_file.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as destination:

                writer = csv.writer(
                    destination,
                    delimiter="\\t",
                    lineterminator="\\n",
                )

                writer.writerow(
                    [
                        "contrast",
                        "target_condition",
                        "reference_condition",
                        *annotated_header,
                    ]
                )

                for row in reader:
                    writer.writerow(
                        [
                            contrast,
                            target,
                            reference,
                            *row,
                        ]
                    )
                    
    event_type_map = {
        "SE": "SE",
        "A5SS": "A5SS",
        "A3SS": "A3SS",
        "MXE": "MXE",
        "RI": "RI",
    }

    gene_events = defaultdict(list)

    # Use only JCEC results so that each rMATS event is counted once.
    # Reading both JC and JCEC would duplicate the same events.
    for annotated_file in sorted(
        output_dir.glob("*.MATS.JCEC.annotated.tsv")
    ):
        event_type = annotated_file.name.split(".")[0]

        if event_type not in event_type_map:
            continue

        with annotated_file.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as source:

            reader = csv.DictReader(
                source,
                delimiter="\\t",
            )

            for row in reader:
                gene_id = (
                    row.get("GeneID", "")
                    .strip()
                    .strip('"')
                )

                gene_name = (
                    row.get("geneSymbol", "")
                    .strip()
                    .strip('"')
                )

                if not gene_id and not gene_name:
                    continue

                try:
                    effect = float(
                        row.get(
                            "dPSI_target_minus_reference",
                            "nan",
                        )
                    )
                except (TypeError, ValueError):
                    effect = math.nan

                try:
                    fdr = float(
                        row.get("FDR", "nan")
                    )
                except (TypeError, ValueError):
                    fdr = math.nan

                gene_events[
                    (gene_id, gene_name)
                ].append(
                    {
                        "event_id": row.get(
                            "event_id",
                            "",
                        ),
                        "event_type": event_type_map[
                            event_type
                        ],
                        "effect": effect,
                        "significant": (
                            not math.isnan(fdr)
                            and fdr < 0.05
                        ),
                    }
                )

    gene_summary_file = (
        output_dir/"RMATS_gene_summary.tsv"
    )

    with gene_summary_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as destination:

        writer = csv.DictWriter(
            destination,
            fieldnames=[
                "gene_id",
                "gene_name",
                "method",
                "detected",
                "significant",
                "n_events",
                "n_significant_events",
                "max_abs_effect",
                "effect_direction",
                "event_types",
            ],
            delimiter="\\t",
            lineterminator="\\n",
        )

        writer.writeheader()

        rows = []

        for (
            gene_id,
            gene_name,
        ), events in gene_events.items():

            significant_events = [
                event
                for event in events
                if event["significant"]
            ]

            direction_events = (
                significant_events
                if significant_events
                else events
            )

            effects = [
                event["effect"]
                for event in direction_events
                if not math.isnan(event["effect"])
                and event["effect"] != 0
            ]

            has_positive = any(
                effect > 0
                for effect in effects
            )

            has_negative = any(
                effect < 0
                for effect in effects
            )

            if has_positive and has_negative:
                effect_direction = "mixed"
            elif has_positive:
                effect_direction = "increased"
            elif has_negative:
                effect_direction = "decreased"
            else:
                effect_direction = "none"

            absolute_effects = [
                abs(event["effect"])
                for event in events
                if not math.isnan(event["effect"])
            ]

            rows.append(
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "method": "rMATS",
                    "detected": True,
                    "significant": bool(
                        significant_events
                    ),
                    "n_events": len(events),
                    "n_significant_events": len(
                        significant_events
                    ),
                    "max_abs_effect": (
                        round(
                            max(absolute_effects),
                            3,
                        )
                        if absolute_effects
                        else ""
                    ),
                    "effect_direction":
                        effect_direction,
                    "event_types": "|".join(
                        sorted(
                            {
                                event["event_type"]
                                for event in events
                            }
                        )
                    ),
                }
            )

        rows.sort(
            key=lambda row: (
                not row["significant"],
                -row["n_significant_events"],
                -(
                    row["max_abs_effect"]
                    if isinstance(
                        row["max_abs_effect"],
                        float,
                    )
                    else 0
                ),
                row["gene_name"],
            )
        )

        writer.writerows(rows)
    PYTHON
    """

}