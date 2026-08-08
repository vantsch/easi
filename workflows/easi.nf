/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW      } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE                } from '../modules/nf-core/trimgalore/main'
include { STAR_ALIGN                } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX            } from '../modules/nf-core/samtools/index/main'
include { BEDOPS_GTF2BED            } from '../modules/nf-core/bedops/gtf2bed/main'
include { AGAT_CONVERTSPGXF2GXF     } from '../modules/nf-core/agat/convertspgxf2gxf/main'
include { RSEQC_INFEREXPERIMENT     } from '../modules/nf-core/rseqc/inferexperiment/main'
include { RSEQC_STRANDNESS          } from '../modules/local/rseqc/strandness/main'
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
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_easi_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EASI {

    take:
    ch_samplesheet
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
        DEFINE CONTRAST
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (!params.contrast_reference || !params.contrast_target) {
        error """
        A differential comparison requires both:

          --contrast_reference
          --contrast_target

        Example:

          --contrast_reference HBR
          --contrast_target UHR
        """.stripIndent()
    }

    def contrast_name = params.contrast_name ?:
        "${params.contrast_target}_vs_${params.contrast_reference}"

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SHARED PREPROCESSING AND ALIGNMENT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW(ch_samplesheet)

    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_RAW.out.zip.map { _meta, file -> file }
    )

    //
    // MODULE: Trim reads
    //
    TRIMGALORE(ch_samplesheet)

    ch_multiqc_files = ch_multiqc_files.mix(
        TRIMGALORE.out.log.map { _meta, file -> file }
    )

    //
    // MODULE: Run FastQC on trimmed reads
    //
    FASTQC_TRIMMED(TRIMGALORE.out.reads)

    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_TRIMMED.out.zip.map { _meta, file -> file }
    )

    //
    // MODULE: Align trimmed reads with STAR
    //
    ch_star_index = Channel.value([
        [:],
        file(params.star_index, checkIfExists: true)
    ])

    ch_gtf = Channel.value([
        [:],
        file(params.gtf, checkIfExists: true)
    ])

    STAR_ALIGN(
        TRIMGALORE.out.reads,
        ch_star_index,
        ch_gtf,
        false
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        STAR_ALIGN.out.log_final.map { _meta, file -> file }
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        STAR_ALIGN.out.log_out.map { _meta, file -> file }
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        STAR_ALIGN.out.log_progress.map { _meta, file -> file }
    )

    //
    // MODULE: Index STAR BAM files with samtools
    //
    SAMTOOLS_INDEX(
        STAR_ALIGN.out.bam_sorted_aligned
    )

    ch_bam_bai = STAR_ALIGN.out.bam_sorted_aligned
        .join(SAMTOOLS_INDEX.out.index)
        .map { meta, bam, bai ->
            tuple(meta, bam, bai)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SHARED QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // MODULE: Convert GTF annotation to BED12 for RSeQC
    //
    BEDOPS_GTF2BED(ch_gtf)

    ch_bed12 = BEDOPS_GTF2BED.out.bed
        .map { meta, bed ->
            bed
        }

    //
    // MODULE: Infer RNA-seq library strandedness
    //
    RSEQC_INFEREXPERIMENT(
        ch_bam_bai,
        ch_bed12
    )

    RSEQC_STRANDNESS(
        RSEQC_INFEREXPERIMENT.out.txt
            .map { meta, report -> report }
            .collect()
    )

    ch_regtools_strandness = RSEQC_STRANDNESS.out.strandness
        .map { strandness ->

            def normalized = strandness.trim()

            def regtools_codes = [
                'REVERSE'   : 'RF',
                'FORWARD'   : 'FR',
                'UNSTRANDED': 'XS'
            ]

            if (!regtools_codes.containsKey(normalized)) {
                error "Unsupported inferred strandness: ${normalized}"
            }

            regtools_codes[normalized]
        }

    ch_salmon_lib_type = RSEQC_STRANDNESS.out.strandness
        .map { strandness ->

            def normalized = strandness.trim()

            def salmon_codes = [
                'REVERSE'   : 'ISR',
                'FORWARD'   : 'ISF',
                'UNSTRANDED': 'IU'
            ]

            if (!salmon_codes.containsKey(normalized)) {
                error "Unsupported inferred strandness for Salmon: ${normalized}"
            }

            salmon_codes[normalized]
        }

    ch_multiqc_files = ch_multiqc_files.mix(
        RSEQC_INFEREXPERIMENT.out.txt.map { meta, file -> file }
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EXPRESSION ANALYSIS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.run_expression.toBoolean()) {

        //
        // MODULE: Count reads with featureCounts
        //
        ch_bam_for_featurecounts = STAR_ALIGN.out.bam_sorted_aligned
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
            SUBREAD_FEATURECOUNTS.out.summary.map { _meta, file -> file }
        )

        //
        // MODULE: Merge featureCounts outputs
        //
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

        //
        // MODULE: Differential expression with DESeq2
        //
        ch_deseq2_contrast = Channel.value([
            [id: contrast_name],
            params.contrast_variable,
            params.contrast_reference,
            params.contrast_target,
            params.deseq2_formula,
            params.deseq2_comparison
        ])

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
            ch_deseq2_input,
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

        ch_splicing_contrast = Channel.value([
            contrast : contrast_name,
            treatment: params.contrast_target,
            control  : params.contrast_reference
        ])

        /*
        --------------------------------------------------------------------------------
            MAJIQ
        --------------------------------------------------------------------------------
        */

        if (params.run_majiq.toBoolean()) {

            //
            // MODULE: Convert GTF annotation to GFF3
            //
            AGAT_CONVERTSPGXF2GXF(ch_gtf)

            ch_gff3 = AGAT_CONVERTSPGXF2GXF.out.output_gff

            //
            // SUBWORKFLOW: Alternative splicing analysis with MAJIQ
            //
            MAJIQ(
                ch_gff3,
                STAR_ALIGN.out.bam_sorted_aligned,
                ch_splicing_contrast
            )
        }

        /*
        --------------------------------------------------------------------------------
            rMATS
        --------------------------------------------------------------------------------
        */

        ch_rmats_strandness = RSEQC_STRANDNESS.out.strandness
            .map { strandness ->

                def normalized = strandness.trim().toUpperCase()

                def rmats_strandness = [
                    'REVERSE'   : 'reverse',
                    'FORWARD'   : 'forward',
                    'UNSTRANDED': 'unstranded'
                ]

                if (!rmats_strandness.containsKey(normalized)) {
                    error "Unsupported inferred strandness for rMATS: ${normalized}"
                }

                rmats_strandness[normalized]
            }

        ch_rmats_bam = STAR_ALIGN.out.bam_sorted_aligned
            .combine(ch_rmats_strandness)
            .map { meta, bam, strandedness ->
                tuple(
                    meta + [strandedness: strandedness],
                    bam
                )
            }

        if (params.run_rmats.toBoolean()) {

            //
            // SUBWORKFLOW: Alternative splicing analysis with rMATS
            //
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

            //
            // MODULE: Extract splice junctions with RegTools
            //
            REGTOOLS_JUNCTIONSEXTRACT(
                ch_bam_bai,
                ch_regtools_strandness
            )

            //
            // MODULE: Construct LeafCutter clustering input
            //
            ch_leafcutter_cluster_input = REGTOOLS_JUNCTIONSEXTRACT.out.junc
                .collect(flat: false)
                .map { collected ->

                    def junc_files = collected.collect { meta, junc ->
                        junc
                    }

                    def groups = collected.collect { meta, junc ->
                        [
                            meta.id,
                            meta.condition
                        ]
                    }

                    tuple(
                        [id: contrast_name],
                        junc_files,
                        groups
                    )
                }

            //
            // MODULE: Cluster splice junctions with LeafCutter
            //
            LEAFCUTTER_CLUSTER(
                ch_leafcutter_cluster_input
            )

            //
            // MODULE: LeafCutter differential-splicing analysis
            //
            ch_leafcutter_ds_input = LEAFCUTTER_CLUSTER.out.numerators
                .join(LEAFCUTTER_CLUSTER.out.groups_file)
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

            //
            // MODULE: Annotate and summarise LeafCutter results
            //
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

            ISA(
                TRIMGALORE.out.reads,
                ch_salmon_lib_type,
                file(params.input, checkIfExists: true),
                contrast_name
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
            ch_rmats = RMATS.out.gene_summary
            ch_leaf = LEAFCUTTER_DOWNSTREAM.out.gene_summary
            ch_isa = ISA.out.gene_summary


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
                .map { files ->

                    def junction_file = files.find { file ->
                        file.name == 'junctions.tsv'
                    }

                    if (!junction_file) {
                        error """
                        MAJIQ delta-PSI modulize output does not contain
                        junctions.tsv for contrast '${contrast_name}'.
                        """.stripIndent()
                    }

                    tuple(
                        [id: contrast_name],
                        junction_file
                    )
                }

            /*
             * Join the three detailed result tables by contrast metadata.
             */
            ch_splice_region_integration = ch_majiq_junctions
                .join(
                    RMATS.out.jc_annotated,
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
             * Collect the sample metadata, BAM files and BAM indexes into
             * aligned lists for the contrast-level sashimi process.
             *
             * collect(flat: false) is important because each
             * (meta, bam, bai) tuple must remain intact.
             */
            ch_sashimi_samples = ch_bam_bai
                .collect(
                    flat: false
                )
                .map { rows ->

                    def ordered_rows = rows.sort { left, right ->
                        left[0].id.toString()
                            <=> right[0].id.toString()
                    }

                    def sample_metadata = ordered_rows.collect { row ->
                        row[0]
                    }

                    def bam_files = ordered_rows.collect { row ->
                        row[1]
                    }

                    def bai_files = ordered_rows.collect { row ->
                        row[2]
                    }

                    tuple(
                        sample_metadata,
                        bam_files,
                        bai_files
                    )
                }

            /*
             * Combine the ranked event table with the complete set of
             * sample BAMs and indexes.
             *
             * Result:
             *   meta,
             *   ranked_events,
             *   sample_metadata,
             *   BAM list,
             *   BAI list
             */
            ch_sashimi_input =
                SPLICE_REGION_INTEGRATION.out.ranked_splicing_events
                    .combine(
                        ch_sashimi_samples
                    )
                    .map {
                        meta,
                        ranked_events,
                        sample_metadata,
                        bam_files,
                        bai_files ->

                        tuple(
                            meta,
                            ranked_events,
                            sample_metadata,
                            bam_files,
                            bai_files
                        )
                    }

            /*
             * SASHIMI_PLOT expects the GTF as a plain path rather than
             * the (meta, path) tuple used by other modules.
             */
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
        methodsDescriptionText(ch_multiqc_custom_methods_description)
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