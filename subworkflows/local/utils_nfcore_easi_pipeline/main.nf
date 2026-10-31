//
// Subworkflow with functionality specific to the nf-core/easi pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version
    validate_params
    monochrome_logs
    nextflow_cli_args
    outdir
    input
    help
    help_full
    show_hidden

    main:

    ch_versions = channel.empty()

    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    def before_text = ""
    def after_text = ""

    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/easi ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""

    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/easi/blob/master/CITATIONS.md
"""

    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    validateInputParameters()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PARSE FASTQ AND REUSABLE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        The field order follows assets/schema_input.json.

        Fields marked as metadata are combined into meta.

        Expected output:

        meta, fastq_1, fastq_2, bam, bai, salmon_quant

        Metadata includes:
            id
            condition
            replicate
            donor (optional)
            batch (optional)
            strandedness (optional, defaults to auto)

        All rows are grouped by sample ID.

        FASTQ samples retain the original format:
            tuple(meta, reads)

        Reusable inputs:
            tuple(meta, bam, bai)
            tuple(meta, salmon_quant_directory)

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    channel
        .fromList(
            samplesheetToList(
                input,
                "${projectDir}/assets/schema_input.json"
            )
        )
        .map { meta, fastq_1, fastq_2, bam, bai, salmon_quant ->

            /*
                Identify available FASTQ files.
            */

            def has_fastq_1 = hasSamplesheetValue(fastq_1)
            def has_fastq_2 = hasSamplesheetValue(fastq_2)

            def reads = []

            if (has_fastq_1) {

                reads.add(fastq_1)

                if (has_fastq_2) {
                    reads.add(fastq_2)
                }
            }

            /*
                Preserve sample metadata.

                Add a default strandedness value when the
                optional samplesheet column is absent.
            */

            def sample_meta = new LinkedHashMap(meta)

            sample_meta.strandedness = (
                meta.strandedness ?: 'auto'
            ).toString().trim().toLowerCase()

            /*
                Determine single-end/paired-end status
                when FASTQs are available.
            */

            /*
            Sequencing layout.

            When FASTQs are available, determine the layout
            directly from the presence of read 2.

            When FASTQs are absent, use the declared
            single_end metadata.
            */

            if (has_fastq_1) {

                def inferred_single_end = !has_fastq_2

                if (hasSamplesheetValue(meta.single_end) &&
                    meta.single_end != inferred_single_end) {
                    error(
                        "Sample ${meta.id}: declared single_end=${meta.single_end} " +
                        "disagrees with the supplied FASTQ files."
                    )
                }

                sample_meta.single_end = inferred_single_end

            } else {

                if (!(meta.single_end instanceof Boolean)) {
                    error(
                        "Sample ${meta.id}: single_end must be true or false " +
                        "when FASTQ files are not provided."
                    )
                }

                sample_meta.single_end = meta.single_end
            }

            /*
                Construct sample record.
            */

            def sample_record = [
                meta        : sample_meta,
                reads       : reads,
                bam         : hasSamplesheetValue(bam) ? bam : null,
                bai         : hasSamplesheetValue(bai) ? bai : null,
                salmon_quant: hasSamplesheetValue(salmon_quant) ? salmon_quant : null
            ]

            return tuple(
                meta.id,
                sample_record
            )
        }
        .groupTuple()
        .map { sample_id, records ->
            validateInputSamplesheet(sample_id, records)
        }
        .set { ch_all_samples }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FASTQ PREPROCESSING INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Preprocess reads when either:

        - No reusable BAM is available
        - No reusable Salmon quantification is available

        Fully processed samples do not enter this channel.

        Output format:
            tuple(meta, reads)

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_samplesheet = ch_all_samples
        .filter { sample ->
            sample.bam == null || sample.salmon_quant == null
        }
        .map { sample ->
            tuple(
                sample.meta,
                sample.reads
            )
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREVIOUSLY ALIGNED BAM INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_existing_bams = ch_all_samples
        .filter { sample ->
            sample.bam != null
        }
        .map { sample ->
            tuple(
                sample.meta,
                sample.bam,
                sample.bai
            )
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREVIOUS SALMON QUANTIFICATIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_existing_salmon = ch_all_samples
        .filter { sample ->
            sample.salmon_quant != null
        }
        .map { sample ->
            tuple(
                sample.meta,
                sample.salmon_quant
            )
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:

    samplesheet     = ch_samplesheet
    existing_bams   = ch_existing_bams
    existing_salmon = ch_existing_salmon
    all_samples     = ch_all_samples
    versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    monochrome_logs
    multiqc_report

    main:

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )

    def multiqc_reports = multiqc_report.toList()

    workflow.onComplete {

        if (email || email_on_fail) {

            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {

        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateInputParameters() {
    genomeExistsError()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK WHETHER AN OPTIONAL SAMPLESHEET VALUE EXISTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nf-schema represents missing optional file columns
    as empty lists.

    Treat null, empty collections and empty strings
    as absent.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def hasSamplesheetValue(value) {

    if (value == null) {
        return false
    }

    if (value instanceof Collection) {
        return !value.isEmpty()
    }

    return !value.toString().trim().isEmpty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE AND COMBINE SAMPLE RECORDS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Supports:
      - One FASTQ pair per sample
      - Multiple sequencing lanes per sample
      - Existing BAM/BAI files
      - Existing Salmon quantifications
      - Optional donor and batch metadata
      - Optional per-sample strandedness

    Multiple FASTQ rows are supported.

    For now, reusable BAM/Salmon inputs must occupy
    one row per sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateInputSamplesheet(sample_id, records) {

    /*
        Check that sample records exist.
    */

    if (!records || records.isEmpty()) {
        error "Sample ${sample_id}: no input records found."
    }

    def first = records[0]

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE STRANDEDNESS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def allowed_strandedness = [
        'auto',
        'forward',
        'reverse',
        'unstranded'
    ]

    records.each { record ->

        def strandedness = record.meta.strandedness

        if (!(strandedness in allowed_strandedness)) {

            error(
                "Sample ${sample_id}: invalid strandedness '${strandedness}'. " +
                "Allowed values: ${allowed_strandedness.join(', ')}"
            )
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE METADATA CONSISTENCY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Metadata must agree between sequencing lanes
        belonging to the same sample.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def metadata_fields = [
        'id',
        'condition',
        'replicate',
        'donor',
        'batch',
        'strandedness',
        'single_end'
    ]

    records.each { record ->

        metadata_fields.each { field ->

            def expected = first.meta[field]
            def observed = record.meta[field]

            if (expected != observed) {

                error(
                    "Sample ${sample_id}: inconsistent '${field}' " +
                    "across samplesheet rows."
                )
            }
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE REUSABLE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def has_reusable_inputs = records.any { record ->

        record.bam != null ||
        record.bai != null ||
        record.salmon_quant != null
    }

    if (records.size() > 1 && has_reusable_inputs) {

        error(
            "Sample ${sample_id}: reusable BAM/Salmon inputs " +
            "must be provided in a single samplesheet row."
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE SINGLE-END / PAIRED-END CONSISTENCY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def combined_reads = []

    def endedness = records
        .findAll { record ->
            !record.reads.isEmpty()
        }
        .collect { record ->
            record.meta.single_end
        }
        .unique()

    if (endedness.size() > 1) {

        error(
            "Sample ${sample_id}: multiple FASTQ runs must " +
            "all be single-end or all be paired-end."
        )
    }

    /*
        Combine FASTQ files from multiple sequencing lanes.
    */

    records.each { record ->
        combined_reads.addAll(record.reads)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE BAM AND SALMON INPUT COMBINATIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def bam = first.bam
    def bai = first.bai
    def salmon_quant = first.salmon_quant

    if ((bam == null) != (bai == null)) {

        error(
            "Sample ${sample_id}: BAM and BAI must be " +
            "provided together."
        )
    }

    if (combined_reads.isEmpty() &&
        (bam == null || salmon_quant == null)) {

        error(
            "Sample ${sample_id}: provide FASTQs or both " +
            "reusable BAM/BAI and Salmon quantification."
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FINAL SAMPLE METADATA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def combined_meta = new LinkedHashMap(first.meta)

    if (!endedness.isEmpty()) {
        combined_meta.single_end = endedness[0]
    }

    /*
        Ensure strandedness is present in the final metadata.
    */

    combined_meta.strandedness = (
        combined_meta.strandedness ?: 'auto'
    ).toString().trim().toLowerCase()

    /*
        Return the combined sample record.
    */

    return [
        meta        : combined_meta,
        reads       : combined_reads,
        bam         : bam,
        bai         : bai,
        salmon_quant: salmon_quant
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GET ATTRIBUTE FROM GENOME CONFIGURATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getGenomeAttribute(attribute) {

    if (params.genomes &&
        params.genome &&
        params.genomes.containsKey(params.genome)) {

        if (params.genomes[params.genome].containsKey(attribute)) {

            return params.genomes[params.genome][attribute]
        }
    }

    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE GENOME SELECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def genomeExistsError() {

    if (params.genomes &&
        params.genome &&
        !params.genomes.containsKey(params.genome)) {

        def error_string =
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(', ')}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

        error(error_string)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC METHODS INFORMATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def toolCitationText() {

    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016)",
        "."
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {

    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results from multiple tools. Bioinformatics, 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
    ].join(' ').trim()

    return reference_text
}