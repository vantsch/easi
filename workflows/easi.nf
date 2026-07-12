/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_RAW     } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/fastqc/main'
include { MULTIQC                  } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE               } from '../modules/nf-core/trimgalore/main'
include { STAR_ALIGN               } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX           } from '../modules/nf-core/samtools/index/main'
include { BEDOPS_GTF2BED           } from '../modules/nf-core/bedops/gtf2bed/main'
include { AGAT_CONVERTSPGXF2GXF    } from '../modules/nf-core/agat/convertspgxf2gxf/main'
include { RSEQC_INFEREXPERIMENT    } from '../modules/nf-core/rseqc/inferexperiment/main'
include { FEATURECOUNTS_MERGE      } from '../modules/local/featurecounts_merge/main'
include { SUBREAD_FEATURECOUNTS    } from '../modules/nf-core/subread/featurecounts/main'
include { DESEQ2_DIFFERENTIAL      } from '../modules/local/deseq2/differential/main'
include { MAJIQ_BUILDGFF3          } from '../modules/local/majiq/buildgff3/main'
include { MAJIQ_BUILDSJ            } from '../modules/local/majiq/buildsj/main'
include { MAJIQ_BUILDUPDATE        } from '../modules/local/majiq/buildupdate/main'
include { MAJIQ_PSICOVERAGE        } from '../modules/local/majiq/psicoverage/main'
include { MAJIQ_QUANTIFY           } from '../modules/local/majiq/quantify/main'
include { MAJIQ_DELTAPSI           } from '../modules/local/majiq/deltapsi/main'
include { MAJIQ_HETEROGEN          } from '../modules/local/majiq/heterogen/main'
include { MAJIQ_SGCOVERAGE         } from '../modules/local/majiq/sgcoverage/main'
include { MAJIQ_MODULIZE as QUANTIFY_MODULIZE }     from '../modules/local/majiq/modulize/main'
include { MAJIQ_MODULIZE as DELTAPSI_MODULIZE }     from '../modules/local/majiq/modulize/main'
include { MAJIQ_MODULIZE as HETEROGEN_MODULIZE }    from '../modules/local/majiq/modulize/main'
include { paramsSummaryMap         } from 'plugin/nf-schema'
include { paramsSummaryMultiqc     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText   } from '../subworkflows/local/utils_nfcore_easi_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EASI {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // DEFINE CONTRAST
    //
    if (!params.contrast_reference || !params.contrast_target) {
        error """
        A differential comparison requires both:
          --contrast_reference
          --contrast_target

        Example:
          --contrast_reference HBR --contrast_target UHR
        """.stripIndent()
    }

    def contrast_name = params.contrast_name ?:
        "${params.contrast_target}_vs_${params.contrast_reference}"

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.map { _meta, file -> file })

    //
    // MODULE: Trim reads
    //
    TRIMGALORE(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.map { _meta, file -> file })

    //
    // MODULE: Run FastQC on trimmed reads
    //
    FASTQC_TRIMMED(TRIMGALORE.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.map { _meta, file -> file })

    //
    // MODULE: Align trimmed reads with STAR
    //
    ch_star_index = Channel.value([[:], file(params.star_index, checkIfExists: true)])
    ch_gtf        = Channel.value([[:], file(params.gtf, checkIfExists: true)])

    STAR_ALIGN(
        TRIMGALORE.out.reads,
        ch_star_index,
        ch_gtf,
        false
    )

    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.map { _meta, file -> file })
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_out.map { _meta, file -> file })
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_progress.map { _meta, file -> file })

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

    //
    // MODULE: Convert GTF annotation to BED12 for RSeQC
    //
    BEDOPS_GTF2BED(ch_gtf)

    ch_bed12 = BEDOPS_GTF2BED.out.bed
        .map { meta, bed -> bed }

    //
    // MODULE: Convert GTF annotation to GFF3 for MAJIQ
    //
    AGAT_CONVERTSPGXF2GXF(ch_gtf)

    ch_gff3 = AGAT_CONVERTSPGXF2GXF.out.output_gff

    //
    // MODULE: Infer RNA-seq library strandedness
    //
    RSEQC_INFEREXPERIMENT(
        ch_bam_bai,
        ch_bed12
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        RSEQC_INFEREXPERIMENT.out.txt.map { meta, file -> file }
    )
    
    //
    // MODULE: Count reads with featureCounts
    //
    ch_bam_for_featurecounts = STAR_ALIGN.out.bam_sorted_aligned
        .combine(ch_gtf)
        .map { meta, bam, gtf_meta, gtf ->
        [ meta, [bam], gtf ]
        }

    SUBREAD_FEATURECOUNTS(ch_bam_for_featurecounts)

    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.map { _meta, file -> file })

    //
    // MODULE: Merge featureCounts outputs
    //
    ch_merge_script = Channel.value(
        file("${projectDir}/modules/local/featurecounts_merge/bin/merge_featurecounts.py", checkIfExists: true)
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
            [
                [id: 'counts'],
                file(params.input, checkIfExists: true),
                counts
            ]
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
        file("${projectDir}/modules/local/deseq2/differential/templates/deseq2_report.Rmd")
    )

    //
    // MODULE: Build the MAJIQ reference splice graph
    //
    MAJIQ_BUILDGFF3(
        ch_gff3
    )
    ch_majiq_splicegraph = MAJIQ_BUILDGFF3.out.splicegraph

    //
    // MODULE: Build per-sample MAJIQ splice-junction files
    //
    ch_bam_for_majiq = STAR_ALIGN.out.bam_sorted_aligned
        .map { meta, bam ->
            tuple(meta, [bam])
        }

    ch_majiq_buildsj_input = ch_bam_for_majiq
        .combine(MAJIQ_BUILDGFF3.out.splicegraph)

    MAJIQ_BUILDSJ(
        ch_majiq_buildsj_input
    )
    
    //
    // MODULE: MAJIQ_BUILDUPDATE
    //
    ch_majiq_sj_condition_map = MAJIQ_BUILDSJ.out.sj
        .collect()
        .map { collected ->
            collected.collate(2)
        }
        .map { pairs ->
            pairs.collect { pair ->
                [pair[1].name.toString(), pair[0].condition]
            }
        }

    ch_majiq_sj_files = MAJIQ_BUILDSJ.out.sj
        .collect()
        .map { collected ->
            collected.collate(2)
        }
        .map { pairs ->
            pairs.collect { pair ->
                pair[1]
            }
        }

        MAJIQ_BUILDUPDATE(
            ch_majiq_sj_files,
            ch_majiq_splicegraph,
            ch_majiq_sj_condition_map
        )

    //
    // MODULE: MAJIQ PSICOVERAGE
    //
    ch_majiq_finished_splicegraph = MAJIQ_BUILDUPDATE.out.splicegraph

    ch_majiq_condition_samples_sj = MAJIQ_BUILDSJ.out.sj
        .map { meta, sj_file ->
            tuple(meta.condition, sj_file)
        }
        .groupTuple()
        .combine(ch_majiq_finished_splicegraph)

    MAJIQ_PSICOVERAGE(
        ch_majiq_condition_samples_sj
    )
    //
    // MODULE: MAJIQ QUANTIFY
    //
    MAJIQ_QUANTIFY(
        MAJIQ_PSICOVERAGE.out.psi_coverage
            .combine(MAJIQ_BUILDUPDATE.out.splicegraph)
    )

    //
    // MODULE: MAJIQ DELTAPSI AND HETEROGEN
    //
    ch_majiq_contrasts = Channel.value([
        contrast : contrast_name,
        treatment: params.contrast_target,
        control  : params.contrast_reference
    ])
    
    ch_majiq_contrast_comparison = ch_majiq_contrasts
        .map { it -> [it.treatment, it] }
        .combine(MAJIQ_PSICOVERAGE.out.psi_coverage, by: 0)
        .map { it -> it[1] + [psicov1: it[2]] }
        .map { it -> [it.control, it] }
        .combine(MAJIQ_PSICOVERAGE.out.psi_coverage, by: 0)
        .map { it -> it[1] + [psicov2: it[2]] }

    ch_majiq_contrast_input = ch_majiq_contrast_comparison
        .combine(ch_majiq_finished_splicegraph)
        .map { it ->
            [
                it[0].contrast,
                it[0].treatment,
                it[0].control,
                it[0].psicov1,
                it[0].psicov2,
                it[2]
            ]
        }

    MAJIQ_DELTAPSI(
        ch_majiq_contrast_input
    )

    MAJIQ_HETEROGEN(
        ch_majiq_contrast_input
    )

    //
    // MODULE: MAJIQ SGCOVERAGE
    //

    MAJIQ_SGCOVERAGE(
        ch_majiq_condition_samples_sj
    )
    
    //
    // MODULE: MAJIQ MODULIZE - QUANTIFY
    //
    ch_majiq_modulize_input_quantify = MAJIQ_SGCOVERAGE.out.sgc_files
        .collect { it[1] }
        .combine(
            MAJIQ_PSICOVERAGE.out.psi_coverage
                .collect { it[1] }
        )
        .toList()
        .combine(ch_majiq_finished_splicegraph)

    QUANTIFY_MODULIZE(
        ch_majiq_modulize_input_quantify
    )   

    //
    // MODULE: MAJIQ MODULIZE - DELTAPSI
    //
    ch_majiq_modulize_input_deltapsi = MAJIQ_SGCOVERAGE.out.sgc_files
        .collect { it[1] }
        .combine(
            MAJIQ_DELTAPSI.out.dpsicov.collect()
        )
        .toList()
        .combine(ch_majiq_finished_splicegraph)

    DELTAPSI_MODULIZE(
        ch_majiq_modulize_input_deltapsi
    )

    //
    // MODULE: MAJIQ MODULIZE - HETEROGEN
    //
    ch_majiq_modulize_input_heterogen = MAJIQ_SGCOVERAGE.out.sgc_files
        .collect { it[1] }
        .combine(
            MAJIQ_HETEROGEN.out.hetcov.collect()
        )
        .toList()
        .combine(ch_majiq_finished_splicegraph)

    HETEROGEN_MODULIZE(
        ch_majiq_modulize_input_heterogen
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'easi_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'easi'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
