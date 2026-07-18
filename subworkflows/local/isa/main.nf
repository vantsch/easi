/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SALMON_INDEX          }    from '../../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT          }    from '../../../modules/nf-core/salmon/quant/main'
include { ISA_PART1             }    from '../../../modules/local/isa/part1/main'
include { PFAM_SCAN             }    from '../../../modules/local/pfam/scan/main'
include { PFAM_FORMAT_ISA       }    from '../../../modules/local/pfam/format_isa/main'
include { ISA_PART2             }    from '../../../modules/local/isa/part2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ISA SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISA {

    take:
    reads
    salmon_lib_type
    samplesheet
    contrast_name

    main:

    ch_genome_fasta = Channel.value(
        file(params.fasta, checkIfExists: true)
    )

    ch_transcript_fasta = Channel.value(
        file(params.transcript_fasta, checkIfExists: true)
    )

    ch_isa_gtf = Channel.value(
        file(
            params.isa_gtf ?: params.gtf,
            checkIfExists: true
        )
    )

    SALMON_INDEX(
        ch_genome_fasta,
        ch_transcript_fasta
    )

    SALMON_QUANT(
        reads,
        SALMON_INDEX.out.index,
        ch_isa_gtf,
        ch_transcript_fasta,
        false,
        salmon_lib_type
    )

    ch_isa_part1_input = SALMON_QUANT.out.results
        .collect(flat: false)
        .map { collected ->

            def result_dirs = collected.collect { item ->
                item[1]
            }

            tuple(
                [id: contrast_name],
                result_dirs
            )
        }

    ch_isa_contrast_file = Channel
        .of(
            "treatment,control\n" +
            "${params.contrast_target}," +
            "${params.contrast_reference}\n"
        )
        .collectFile(
            name: "${contrast_name}_isa_contrast.csv",
            newLine: false
        )

    ISA_PART1(
        ch_isa_part1_input,
        ch_isa_gtf,
        ch_transcript_fasta,
        samplesheet,
        ch_isa_contrast_file,
        params.isa_alpha,
        params.isa_dif_cutoff
    )

    if (!params.pfam_hmm) {
        error "ISA requires --pfam_hmm"
    }

    if (!params.pfam_clans) {
        error "ISA requires --pfam_clans"
    }

    def pfam_hmm = file(
        params.pfam_hmm,
        checkIfExists: true
    )

    ch_pfam_db = Channel.value(
        tuple(
            pfam_hmm,
            file("${params.pfam_hmm}.h3f", checkIfExists: true),
            file("${params.pfam_hmm}.h3i", checkIfExists: true),
            file("${params.pfam_hmm}.h3m", checkIfExists: true),
            file("${params.pfam_hmm}.h3p", checkIfExists: true)
        )
    )

    PFAM_SCAN(
        ISA_PART1.out.aa_fasta,
        ch_pfam_db
    )

    ch_pfam_clans = Channel.value(
        file(params.pfam_clans, checkIfExists: true)
    )

    PFAM_FORMAT_ISA(
        PFAM_SCAN.out.domtblout,
        ch_pfam_clans
    )

    ch_isa_part2_input = ISA_PART1.out.switchlist
        .join(
            PFAM_FORMAT_ISA.out.pfam_isa,
            by: 0
        )

    ISA_PART2(
        ch_isa_part2_input
    )

    emit:
    salmon_results             = SALMON_QUANT.out.results
    switchlist_part1           = ISA_PART1.out.switchlist
    isa_part1_switches         = ISA_PART1.out.switches
    isa_part1_isoform_features = ISA_PART1.out.isoform_features
    isa_part1_orf_analysis     = ISA_PART1.out.orf_analysis
    aa_fasta                   = ISA_PART1.out.aa_fasta
    nt_fasta                   = ISA_PART1.out.nt_fasta
    sequence_files             = ISA_PART1.out.sequences
    sequence_manifest          = ISA_PART1.out.sequence_manifest
    pfam_domtblout             = PFAM_SCAN.out.domtblout
    pfam_tblout                = PFAM_SCAN.out.tblout
    pfam_isa                   = PFAM_FORMAT_ISA.out.pfam_isa
    switchlist_part2           = ISA_PART2.out.switchlist
    domain_analysis            = ISA_PART2.out.domain_analysis
    isa_part2_isoform_features = ISA_PART2.out.isoform_features
    switch_plots               = ISA_PART2.out.switch_plots
    consequence_summary        = ISA_PART2.out.consequence_summary
}