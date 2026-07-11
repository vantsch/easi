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
include { FEATURECOUNTS_MERGE      } from '../modules/local/featurecounts_merge/main'
include { DESEQ2_DIFFERENTIAL      } from '../modules/local/deseq2/differential/main'
include { SUBREAD_FEATURECOUNTS    } from '../modules/nf-core/subread/featurecounts/main'
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
        [id: 'UHR_vs_HBR'],
        'condition',
        'HBR',
        'UHR',
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
