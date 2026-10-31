
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EASI SHARED PREPROCESSING

    Shared sample-level preprocessing for EASI.

    Handles:
      - FASTQ quality control
      - Read trimming
      - STAR alignment
      - BAM indexing
      - Previously aligned BAM reuse
      - Per-sample RSeQC strandedness inference
      - Preservation of declared library strandedness
      - Single-end and paired-end Salmon library types
      - Salmon quantification
      - Previously generated Salmon quantification reuse

    Downstream analysis modes consume the same BAM and Salmon outputs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW     } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE               } from '../../../modules/nf-core/trimgalore/main'
include { STAR_ALIGN               } from '../../../modules/nf-core/star/align/main'
include { SAMTOOLS_INDEX           } from '../../../modules/nf-core/samtools/index/main'
include { BEDOPS_GTF2BED           } from '../../../modules/nf-core/bedops/gtf2bed/main'
include { RSEQC_INFEREXPERIMENT    } from '../../../modules/nf-core/rseqc/inferexperiment/main'
include { RSEQC_STRANDNESS         } from '../../../modules/local/rseqc/strandness/main'
include { SALMON_INDEX             } from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT             } from '../../../modules/nf-core/salmon/quant/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SHARED PREPROCESSING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EASI_PREPROCESSING {

    take:
    ch_samplesheet
    ch_existing_bams
    ch_existing_salmon
    ch_all_samples

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        REFERENCE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_star_index = Channel.value([
        [:],
        file(params.star_index, checkIfExists: true)
    ])

    ch_gtf = Channel.value([
        [:],
        file(params.gtf, checkIfExists: true)
    ])

    ch_bed_gtf = ch_gtf

    ch_isa_gtf = Channel.value(
        file(
            params.isa_gtf ?: params.gtf,
            checkIfExists: true
        )
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IDENTIFY SAMPLES REQUIRING COMPUTATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_star_required = ch_all_samples
        .filter { sample ->
            sample.bam == null
        }
        .map { sample ->
            tuple(sample.meta.id, true)
        }

    ch_salmon_required = ch_all_samples
        .filter { sample ->
            sample.salmon_quant == null
        }
        .map { sample ->
            tuple(sample.meta.id, true)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RAW READ QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    FASTQC_RAW(
        ch_samplesheet
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        READ TRIMMING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    TRIMGALORE(
        ch_samplesheet
    )

    ch_trimmed_reads = TRIMGALORE.out.reads

    FASTQC_TRIMMED(
        ch_trimmed_reads
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SELECT READS FOR STAR
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_reads_for_star = ch_trimmed_reads
        .map { meta, reads ->
            tuple(meta.id, meta, reads)
        }
        .join(ch_star_required, failOnDuplicate: true)
        .map { sample_id, meta, reads, required ->
            tuple(meta, reads)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STAR ALIGNMENT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    STAR_ALIGN(
        ch_reads_for_star,
        ch_star_index,
        ch_gtf,
        false
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INDEX NEW BAM FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    SAMTOOLS_INDEX(
        STAR_ALIGN.out.bam_sorted_aligned
    )

    ch_new_bam_bai = STAR_ALIGN.out.bam_sorted_aligned
        .map { meta, bam ->
            tuple(meta.id, meta, bam)
        }
        .join(
            SAMTOOLS_INDEX.out.index.map { meta, bai ->
                tuple(meta.id, bai)
            },
            failOnDuplicate: true
        )
        .map { sample_id, meta, bam, bai ->
            tuple(meta, bam, bai)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COMBINE NEW AND EXISTING BAM FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_all_bam_bai = ch_new_bam_bai
        .mix(ch_existing_bams)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SELECT SAMPLES REQUIRING STRANDEDNESS INFERENCE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Only samples with strandedness = auto undergo RSeQC inference.

        Explicitly declared forward, reverse and unstranded values
        are preserved.

        Both newly generated and previously aligned BAMs are supported.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_auto_bams = ch_all_bam_bai
        .filter { meta, bam, bai ->
            (meta.strandedness ?: 'auto').toString().toLowerCase() == 'auto'
        }

    ch_declared_bams = ch_all_bam_bai
        .filter { meta, bam, bai ->
            (meta.strandedness ?: 'auto').toString().toLowerCase() != 'auto'
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GENERATE BED REFERENCE FOR RSEQC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    BEDOPS_GTF2BED(
        ch_bed_gtf
    )

    ch_bed12 = BEDOPS_GTF2BED.out.bed
        .map { meta, bed ->
            bed
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN RSEQC FOR AUTO SAMPLES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    RSEQC_INFEREXPERIMENT(
        ch_auto_bams,
        ch_bed12
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CLASSIFY STRANDEDNESS FOR EACH SAMPLE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        RSEQC_STRANDNESS processes one report per sample.

        Expected output:
            tuple(meta, classification)

        classification:
            FORWARD
            REVERSE
            UNSTRANDED
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    RSEQC_STRANDNESS(
        RSEQC_INFEREXPERIMENT.out.txt
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        NORMALIZE INFERRED STRANDEDNESS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_inferred_strandedness = RSEQC_STRANDNESS.out.strandness
        .map { meta, result ->

            def classification = result.toString().trim().toLowerCase()

            if (!(classification in ['forward', 'reverse', 'unstranded'])) {
                error """
                Invalid RSeQC strandedness classification.

                Sample: ${meta.id}
                Classification: ${classification}
                """.stripIndent()
            }

            def updated_meta = meta + [
                strandedness: classification
            ]

            tuple(updated_meta.id, updated_meta)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PRESERVE DECLARED STRANDEDNESS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_declared_strandedness = ch_declared_bams
        .map { meta, bam, bai ->

            def classification = meta.strandedness.toString().toLowerCase()

            if (!(classification in ['forward', 'reverse', 'unstranded'])) {
                error """
                Invalid declared strandedness.

                Sample: ${meta.id}
                Strandedness: ${classification}
                """.stripIndent()
            }

            def updated_meta = meta + [
                strandedness: classification
            ]

            tuple(updated_meta.id, updated_meta)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COMBINE INFERRED AND DECLARED STRANDEDNESS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Every sample receives exactly one resolved strandedness value.

        The sample ID is used for downstream joins.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_resolved_strandedness = ch_inferred_strandedness
        .mix(ch_declared_strandedness)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ATTACH RESOLVED METADATA TO BAM FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_resolved_bam_bai = ch_all_bam_bai
        .map { meta, bam, bai ->
            tuple(meta.id, bam, bai)
        }
        .join(
            ch_resolved_strandedness,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        .map { sample_id, bam, bai, meta ->
            tuple(meta, bam, bai)
        }

    ch_resolved_bams = ch_resolved_bam_bai
        .map { meta, bam, bai ->
            tuple(meta, bam)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREPARE PER-SAMPLE STRANDEDNESS FOR REGTOOLS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_regtools_strandness = ch_resolved_strandedness
        .map { sample_id, meta ->

            def regtools_codes = [
                'reverse'   : 'RF',
                'forward'   : 'FR',
                'unstranded': 'XS'
            ]

            def code = regtools_codes[meta.strandedness]

            if (!code) {
                error "Unsupported RegTools strandedness for ${sample_id}: ${meta.strandedness}"
            }

            tuple(meta, code)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SELECT READS FOR SALMON
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_reads_for_salmon = ch_trimmed_reads
        .map { meta, reads ->
            tuple(meta.id, reads)
        }
        .join(
            ch_salmon_required,
            failOnDuplicate: true
        )
        .map { sample_id, reads, required ->
            tuple(sample_id, reads)
        }
        .join(
            ch_resolved_strandedness,
            failOnDuplicate: true
        )
        .map { sample_id, reads, meta ->
            tuple(meta, reads)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONDITIONAL SALMON INDEXING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Only create the Salmon index if at least one sample
        requires new quantification.

        This avoids unnecessary indexing when all samples
        already have reusable Salmon results.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_salmon_needed = ch_reads_for_salmon
        .toList()
        .filter { reads ->
            !reads.isEmpty()
        }

    ch_genome_fasta = ch_salmon_needed
        .map {
            file(params.fasta, checkIfExists: true)
        }

    ch_transcript_fasta = ch_salmon_needed
        .map {
            file(params.transcript_fasta, checkIfExists: true)
        }

    SALMON_INDEX(
        ch_genome_fasta,
        ch_transcript_fasta
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SALMON QUANTIFICATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        The Salmon module derives the library type from:

            meta.single_end
            meta.strandedness

        Single-end:
            unstranded -> U
            forward    -> SF
            reverse    -> SR

        Paired-end:
            unstranded -> IU
            forward    -> ISF
            reverse    -> ISR

        An empty lib_type argument allows the module to use
        per-sample metadata.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_transcript_fasta_quant = Channel.value(
        file(params.transcript_fasta, checkIfExists: true)
    )

    SALMON_QUANT(
        ch_reads_for_salmon,
        SALMON_INDEX.out.index,
        ch_isa_gtf,
        ch_transcript_fasta_quant,
        false,
        ''
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COMBINE NEW AND EXISTING SALMON RESULTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_all_salmon_results = SALMON_QUANT.out.results
        .mix(ch_existing_salmon)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT QUALITY CONTROL OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_qc_files = FASTQC_RAW.out.zip
        .map { meta, file -> file }
        .mix(
            TRIMGALORE.out.log
                .map { meta, file -> file }
        )
        .mix(
            FASTQC_TRIMMED.out.zip
                .map { meta, file -> file }
        )
        .mix(
            STAR_ALIGN.out.log_final
                .map { meta, file -> file }
        )
        .mix(
            STAR_ALIGN.out.log_out
                .map { meta, file -> file }
        )
        .mix(
            STAR_ALIGN.out.log_progress
                .map { meta, file -> file }
        )
        .mix(
            RSEQC_INFEREXPERIMENT.out.txt
                .map { meta, file -> file }
        )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:

    // BAMs with resolved per-sample strandedness
    bam_bai = ch_resolved_bam_bai

    bam = ch_resolved_bams

    // Salmon quantifications
    salmon_results = ch_all_salmon_results

    // Per-sample library metadata
    resolved_strandedness = ch_resolved_strandedness

    // Per-sample RegTools codes
    regtools_strandness = ch_regtools_strandness

    // Trimmed reads
    trimmed_reads = ch_trimmed_reads

    // Quality control
    qc_files = ch_qc_files

    // Newly aligned BAMs
    new_bam_bai = ch_new_bam_bai
}
