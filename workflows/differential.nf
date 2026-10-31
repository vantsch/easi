/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EASI_PREPROCESSING        } from '../subworkflows/local/preprocessing/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { AGAT_CONVERTSPGXF2GXF     } from '../modules/nf-core/agat/convertspgxf2gxf/main'
include { FEATURECOUNTS_MERGE       } from '../modules/local/featurecounts_merge/main'
include { SUBREAD_FEATURECOUNTS     } from '../modules/nf-core/subread/featurecounts/main'
include { DESEQ2_DIFFERENTIAL       } from '../modules/local/deseq2/differential/main'
include { MAJIQ                     } from '../subworkflows/local/majiq/main'
include { RMATS                     } from '../subworkflows/local/rmats/main'
include { REGTOOLS_JUNCTIONSEXTRACT } from '../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTER        } from '../modules/local/leafcutter/cluster/main'
include { LEAFCUTTER_DS             } from '../modules/local/leafcutter/ds/main'
include { LEAFCUTTER_DOWNSTREAM     } from '../modules/local/leafcutter/downstream/main'
include { ISA                       } from '../subworkflows/local/isa/main'
include { GENE_LEVEL_INTEGRATION    } from '../modules/local/integration/gene_level/main'
include { SPLICE_REGION_INTEGRATION } from '../modules/local/integration/splice_region/main'
include { SASHIMI_PLOT              } from '../modules/local/visualization/sashimi/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EASI_DIFFERENTIAL {

    take:
    ch_samplesheet
    ch_existing_bams
    ch_existing_salmon
    ch_all_samples
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE ANALYSIS SELECTION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (
        !params.run_expression.toBoolean() &&
        !params.run_splicing.toBoolean()
    ) {
        error """
        At least one analysis branch must be enabled:

          --run_expression true
          --run_splicing true
        """.stripIndent()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DEFINE CONTRAST(S)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (!params.contrasts && (!params.contrast_reference || !params.contrast_target)) {
        error """
        A differential comparison requires either:

        --contrasts  /path/to/contrasts.csv

        or both:

        --contrast_reference  <reference_condition>
        --contrast_target     <target_condition>

        Example (single):

        --contrast_reference WT --contrast_target KO1

        Example (multi):

        contrasts.csv
            contrast_name,contrast_reference,contrast_target
            KO1_vs_WT,WT,KO1
            KO2_vs_WT,WT,KO2
        """.stripIndent()
    }

    def ch_contrasts

    if (params.contrasts) {
        ch_contrasts = Channel
            .fromPath(params.contrasts, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def name = row.contrast_name ?: "${row.contrast_target}_vs_${row.contrast_reference}"

                tuple(
                    [id: name],
                    row.contrast_reference,
                    row.contrast_target
                )
            }
    } else {
        def contrast_name = params.contrast_name ?:
            "${params.contrast_target}_vs_${params.contrast_reference}"

        ch_contrasts = Channel.of(
            tuple(
                [id: contrast_name],
                params.contrast_reference,
                params.contrast_target
            )
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SHARED PREPROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    EASI_PREPROCESSING(
        ch_samplesheet,
        ch_existing_bams,
        ch_existing_salmon,
        ch_all_samples
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREPROCESSING OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_bam_bai = EASI_PREPROCESSING.out.bam_bai

    ch_bams = EASI_PREPROCESSING.out.bam

    ch_salmon_results = EASI_PREPROCESSING.out.salmon_results

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        REFERENCE ANNOTATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_gtf = Channel.value([
        [:],
        file(params.gtf, checkIfExists: true)
    ])

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREPROCESSING QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_multiqc_files = ch_multiqc_files.mix(
        EASI_PREPROCESSING.out.qc_files
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EXPRESSION ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.run_expression.toBoolean()) {

        /*
        --------------------------------------------------------------------------------
            featureCounts
        --------------------------------------------------------------------------------
        */

        ch_bam_for_featurecounts = ch_bams
            .combine(ch_gtf)
            .map { meta, bam, gtf_meta, gtf ->
                tuple(
                    meta,
                    [bam],
                    gtf
                )
            }

        SUBREAD_FEATURECOUNTS(
            ch_bam_for_featurecounts
        )

        ch_multiqc_files = ch_multiqc_files.mix(
            SUBREAD_FEATURECOUNTS.out.summary
                .map { _meta, file -> file }
        )

        /*
        --------------------------------------------------------------------------------
            Merge featureCounts outputs
        --------------------------------------------------------------------------------
        */

        ch_merge_script = Channel.value(
            file(
                "${projectDir}/modules/local/featurecounts_merge/bin/merge_featurecounts.py",
                checkIfExists: true
            )
        )

        FEATURECOUNTS_MERGE(
            SUBREAD_FEATURECOUNTS.out.counts
                .map { meta, file -> file }
                .collect(),
            ch_merge_script
        )

        /*
        --------------------------------------------------------------------------------
            DESeq2
        --------------------------------------------------------------------------------
        */

        ch_deseq2_contrast = ch_contrasts
            .map { meta, ref, tgt ->
                tuple(
                    meta,
                    params.contrast_variable,
                    ref,
                    tgt,
                    params.deseq2_formula,
                    params.deseq2_comparison
                )
            }

        ch_deseq2_input = FEATURECOUNTS_MERGE.out.matrix
            .map { counts ->
                tuple(
                    [id: 'counts'],
                    file(params.input, checkIfExists: true),
                    counts
                )
            }

        ch_empty_file = Channel.value([
            [id: 'none'],
            []
        ])

        DESEQ2_DIFFERENTIAL(
            ch_deseq2_contrast,
            ch_deseq2_input.first(),
            ch_empty_file,
            ch_empty_file,
            file(
                "${projectDir}/modules/local/deseq2/differential/templates/deseq2_report.Rmd",
                checkIfExists: true
            )
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPLICING ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.run_splicing.toBoolean()) {

        /*
        --------------------------------------------------------------------------------
            Define splicing contrasts
        --------------------------------------------------------------------------------
        */

        ch_splicing_contrast = ch_contrasts
            .map { meta, ref, tgt ->
                [
                    contrast : meta.id,
                    treatment: tgt,
                    control  : ref
                ]
            }

        /*
        --------------------------------------------------------------------------------
            MAJIQ
        --------------------------------------------------------------------------------
        */

        if (params.run_majiq.toBoolean()) {

            /*
             * Convert GTF annotation to GFF3
             */

            AGAT_CONVERTSPGXF2GXF(
                ch_gtf
            )

            ch_gff3 = AGAT_CONVERTSPGXF2GXF.out.output_gff

            /*
             * Alternative splicing analysis with MAJIQ
             */

            MAJIQ(
                ch_gff3,
                ch_bams,
                ch_splicing_contrast
            )
        }

        /*
        --------------------------------------------------------------------------------
            rMATS
        --------------------------------------------------------------------------------
        */

        /*
         * BAM metadata already contains the resolved
         * strandedness for each sample.
         *
         * This includes both newly generated and reused BAMs.
         */

        ch_rmats_bam = ch_bams

        if (params.run_rmats.toBoolean()) {

            RMATS(
                ch_rmats_bam,
                ch_gtf,
                ch_splicing_contrast,
                params.rmats_read_len,
                params.rmats_splice_diff_cutoff,
                params.rmats_novel_splice_site,
                params.rmats_min_intron_len,
                params.rmats_max_exon_len,
                params.rmats_paired_stats
            )
        }

        /*
        --------------------------------------------------------------------------------
            LeafCutter
        --------------------------------------------------------------------------------
        */

        if (params.run_leafcutter.toBoolean()) {

            /*
             * Associate each BAM with its own RegTools
             * strand specificity using the unique sample ID.
             */

            ch_regtools_input = ch_bam_bai
                .map { meta, bam, bai ->
                    tuple(meta.id, meta, bam, bai)
                }
                .join(
                    EASI_PREPROCESSING.out.regtools_strandness
                        .map { meta, strand ->
                            tuple(meta.id, strand)
                        },
                    failOnMismatch: true,
                    failOnDuplicate: true
                )
                .map { sample_id, meta, bam, bai, strand ->
                    tuple(meta, bam, bai, strand)
                }

            /*
             * Extract splice junctions with RegTools
             */

            REGTOOLS_JUNCTIONSEXTRACT(
                ch_regtools_input
            )

            /*
            * Build LeafCutter input with deterministic sample ordering.
            *
            * Target samples are listed first, followed by reference samples.
            * Within each group, samples are sorted by their unique IDs.
            *
            * This ensures reproducible group orientation and input ordering.
            */

            ch_leafcutter_cluster_input = ch_contrasts
                .combine(
                    REGTOOLS_JUNCTIONSEXTRACT.out.junc
                        .collect(flat: false)
                        .map { it -> [it] }
                )
                .map { c_meta, ref, tgt, collected ->

                    // Select samples belonging to this contrast.
                    def selected = collected
                        .findAll { meta, junc ->
                            meta.condition == ref || meta.condition == tgt
                        }

                    // Deterministic ordering:
                    // 1. Target group
                    // 2. Reference group
                    // 3. Sample ID within each group
                    def ordered = selected.sort { a, b ->

                        def rank_a = a[0].condition == tgt ? 0 : 1
                        def rank_b = b[0].condition == tgt ? 0 : 1

                        def group_comparison = rank_a <=> rank_b

                        if (group_comparison != 0) {
                            return group_comparison
                        }

                        return a[0].id.toString() <=> b[0].id.toString()
                    }

                    def junc_files = ordered.collect { meta, junc ->
                        junc
                    }

                    def groups = ordered.collect { meta, junc ->
                        [meta.id, meta.condition]
                    }

                    tuple(
                        c_meta,
                        junc_files,
                        groups
                    )
                }

            /*
             * Cluster splice junctions
             */

            LEAFCUTTER_CLUSTER(
                ch_leafcutter_cluster_input
            )

            /*
             * LeafCutter differential-splicing analysis
             */

            ch_leafcutter_ds_input = LEAFCUTTER_CLUSTER.out.numerators
                .join(
                    LEAFCUTTER_CLUSTER.out.groups_file
                )
                .map { meta, counts, groups ->
                    tuple(
                        meta,
                        counts,
                        groups
                    )
                }

            LEAFCUTTER_DS(
                ch_leafcutter_ds_input
            )

            /*
             * Annotate and summarise LeafCutter results
             */

            ch_leafcutter_downstream = LEAFCUTTER_DS.out.cluster_significance
                .join(
                    LEAFCUTTER_DS.out.effect_sizes,
                    by: 0
                )
                .join(
                    LEAFCUTTER_CLUSTER.out.numerators,
                    by: 0
                )
                .join(
                    LEAFCUTTER_CLUSTER.out.groups_file,
                    by: 0
                )
                .combine(ch_gtf)
                .map {
                    meta,
                    cluster_significance,
                    effect_sizes,
                    counts,
                    groups,
                    gtf_meta,
                    gtf ->

                    tuple(
                        meta,
                        cluster_significance,
                        effect_sizes,
                        counts,
                        groups,
                        gtf
                    )
                }

            LEAFCUTTER_DOWNSTREAM(
                ch_leafcutter_downstream
            )
        }

        /*
        --------------------------------------------------------------------------------
            IsoformSwitchAnalyzeR
        --------------------------------------------------------------------------------
        */

        if (params.run_isa.toBoolean()) {

            ch_isa_contrasts = ch_contrasts
                .map { meta, ref, tgt ->
                    [
                        id: meta.id,
                        control: ref,
                        treatment: tgt
                    ]
                }

            ISA(
                ch_salmon_results,
                file(params.input, checkIfExists: true),
                ch_isa_contrasts
            )
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            GENE-LEVEL INTEGRATION
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        if (
            params.run_expression.toBoolean() &&
            params.run_splicing.toBoolean() &&
            params.run_majiq.toBoolean() &&
            params.run_rmats.toBoolean() &&
            params.run_leafcutter.toBoolean() &&
            params.run_isa.toBoolean()
        ) {

            /*
            --------------------------------------------------------------------------------
                Combine DESeq2 results and standardized splicing gene summaries
            --------------------------------------------------------------------------------
            */

            ch_deseq2_results = DESEQ2_DIFFERENTIAL.out.results

            ch_majiq = MAJIQ.out.gene_summary

            /*
             * Restore canonical contrast names for rMATS gene results.
             */

            ch_rmats = RMATS.out.gene_summary
                .map { meta, result_file ->
                    tuple(meta.id, result_file)
                }
                .join(
                    ch_splicing_contrast.map { c ->
                        tuple(
                            "${c.treatment}_vs_${c.control}",
                            [id: c.contrast]
                        )
                    }
                )
                .map { old_id, result_file, canonical_meta ->
                    tuple(canonical_meta, result_file)
                }

            ch_leaf = LEAFCUTTER_DOWNSTREAM.out.gene_summary

            /*
             * Normalize ISA metadata for integration joins.
             */

            ch_isa = ISA.out.gene_summary
                .map { meta, result_file ->
                    tuple([id: meta.id], result_file)
                }

            /*
             * Join expression and splicing gene summaries.
             */

            ch_gene_level_integration = ch_deseq2_results
                .join(
                    ch_majiq,
                    by: 0
                )
                .join(
                    ch_rmats,
                    by: 0
                )
                .join(
                    ch_leaf,
                    by: 0
                )
                .join(
                    ch_isa,
                    by: 0
                )

            ch_gene_ranking_config = Channel.value(
                file(
                    "${projectDir}/conf/gene_ranking.yml",
                    checkIfExists: true
                )
            )

            GENE_LEVEL_INTEGRATION(
                ch_gene_level_integration,
                ch_gene_ranking_config
            )
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SPLICE-REGION INTEGRATION AND SASHIMI VISUALISATION
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        if (
            params.run_majiq.toBoolean() &&
            params.run_rmats.toBoolean() &&
            params.run_leafcutter.toBoolean()
        ) {

            /*
             * Extract the MAJIQ junction table from the delta-PSI
             * modulize output.
             */

            ch_majiq_junctions = MAJIQ.out.deltapsi_tsv_files
                .map { meta, files ->

                    def junction_file = files.find { file ->
                        file.name == 'junctions.tsv'
                    }

                    if (!junction_file) {
                        error """
                        MAJIQ delta-PSI modulize output does not contain
                        junctions.tsv for contrast '${meta.id}'.
                        """.stripIndent()
                    }

                    tuple(meta, junction_file)
                }

            /*
             * Join the three detailed result tables by contrast metadata.
             */

            /*
             * Restore canonical IDs for rMATS junction tables.
             */

            ch_rmats_jc_canonical = RMATS.out.jc_annotated
                .map { meta, result_files ->
                    tuple(meta.id, result_files)
                }
                .join(
                    ch_splicing_contrast.map { c ->
                        tuple(
                            "${c.treatment}_vs_${c.control}",
                            [id: c.contrast]
                        )
                    }
                )
                .map { old_id, result_files, canonical_meta ->
                    tuple(canonical_meta, result_files)
                }

            ch_splice_region_integration = ch_majiq_junctions
                .join(
                    ch_rmats_jc_canonical,
                    by: 0
                )
                .join(
                    LEAFCUTTER_DOWNSTREAM.out.introns_annotated,
                    by: 0
                )

            SPLICE_REGION_INTEGRATION(
                ch_splice_region_integration
            )

            /*
             * Combine ranked events with all sample BAMs and indexes.
             *
             * Result:
             *   meta,
             *   ranked_events,
             *   sample_metadata,
             *   BAM list,
             *   BAI list
             */

            ch_sashimi_input = SPLICE_REGION_INTEGRATION.out.ranked_splicing_events
                .join(
                    ch_contrasts.map { meta, ref, tgt ->
                        tuple(meta, ref, tgt)
                    },
                    by: 0
                )
                .combine(
                    ch_bam_bai
                        .collect(flat: false)
                        .map { it -> [it] }
                )
                .map { meta, ranked_events, ref, tgt, all_bams ->

                    def ordered_rows = all_bams
                        .findAll { row ->
                            row[0].condition == ref || row[0].condition == tgt
                        }
                        .sort { left, right ->
                            left[0].id.toString() <=> right[0].id.toString()
                        }

                    tuple(
                        meta,
                        ranked_events,
                        ordered_rows.collect { it[0] },
                        ordered_rows.collect { it[1] },
                        ordered_rows.collect { it[2] }
                    )
                }

            ch_sashimi_gtf = ch_gtf
                .map { _gtf_meta, gtf ->
                    gtf
                }

            SASHIMI_PLOT(
                ch_sashimi_input,
                ch_sashimi_gtf
            )
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SOFTWARE VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def topic_versions = channel.topic('versions')
        .distinct()
        .branch { entry ->
            versions_file : entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [
                process[process.lastIndexOf(':') + 1..-1],
                "  ${tool}: ${version}"
            ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()

            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(
        ch_versions.mix(topic_versions.versions_file)
    )
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_' + 'easi_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true
        )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MULTIQC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    def ch_summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: 'nextflow_schema.json'
    )

    def ch_workflow_summary = channel.value(
        paramsSummaryMultiqc(ch_summary_params)
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(
            name: 'workflow_summary_mqc.yaml'
        )
    )

    def ch_multiqc_custom_methods_description =
        multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file(
                "${projectDir}/assets/methods_description_template.yml",
                checkIfExists: true
            )

    def ch_methods_description = channel.value(
        ch_multiqc_custom_methods_description.text
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC(
        ch_multiqc_files
            .flatten()
            .collect()
            .map { files ->
                [
                    [id: 'easi'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file(
                            "${projectDir}/assets/multiqc_config.yml",
                            checkIfExists: true
                        ),
                    multiqc_logo
                        ? file(multiqc_logo, checkIfExists: true)
                        : [],
                    [],
                    []
                ]
            }
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:

    multiqc_report = MULTIQC.out.report
        .map { _meta, report -> [report] }
        .toList()

    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/