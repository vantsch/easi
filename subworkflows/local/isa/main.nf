
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISA_PART1             } from '../../../modules/local/isa/part1/main'
include { PFAM_SCAN             } from '../../../modules/local/pfam/scan/main'
include { PFAM_FORMAT_ISA       } from '../../../modules/local/pfam/format_isa/main'
include { ISA_PART2             } from '../../../modules/local/isa/part2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ISA SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISA {

    take:
    salmon_results
    samplesheet
    contrast_name

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        REFERENCE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_transcript_fasta = Channel.value(
        file(params.transcript_fasta, checkIfExists: true)
    )

    ch_isa_gtf = Channel.value(
        file(
            params.isa_gtf ?: params.gtf,
            checkIfExists: true
        )
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PREPARE SALMON RESULTS FOR EACH CONTRAST
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Salmon quantification is performed by EASI_PREPROCESSING.

        This workflow accepts:
          - Newly generated Salmon results
          - Previously generated Salmon results

        Expected channel structure:
            tuple(meta, salmon_quantification_directory)

        Results are filtered separately for each contrast.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_isa_part1_input = salmon_results
        .collect(flat: false)
        .combine(contrast_name)
        .map { combined ->

            def contrast = combined[-1]
            def collected = combined[0..-2]

            def selected = collected.findAll { item ->
                item[0].condition in [
                    contrast.control,
                    contrast.treatment
                ]
            }

            def controls = selected.findAll { item ->
                item[0].condition == contrast.control
            }

            def treatments = selected.findAll { item ->
                item[0].condition == contrast.treatment
            }

            if (!controls || !treatments) {
                error "ISA contrast ${contrast.id}: missing samples."
            }

            tuple(
                contrast,
                selected.collect { item -> item[1] }
            )
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ISA PART 1
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ISA_PART1(
        ch_isa_part1_input,
        ch_isa_gtf,
        ch_transcript_fasta,
        samplesheet,
        params.isa_alpha,
        params.isa_dif_cutoff
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PFAM REQUIREMENTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (!params.pfam_hmm) {
        error "ISA requires --pfam_hmm"
    }

    if (!params.pfam_clans) {
        error "ISA requires --pfam_clans"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PFAM DATABASE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

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

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PFAM SCAN
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

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

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ISA PART 2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_isa_part2_input = ISA_PART1.out.switchlist
        .join(
            PFAM_FORMAT_ISA.out.pfam_isa,
            by: 0
        )

    ISA_PART2(
        ch_isa_part2_input,
        ch_isa_gtf
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:

    salmon_results             = salmon_results
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
    gene_summary               = ISA_PART2.out.gene_summary
}
