/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_BAMLIST            } from '../../../modules/local/rmats/create_bamlist/main'
include { RMATS_RUN                 } from '../../../modules/local/rmats/run/main'
include { RMATS_SUMMARIZE           } from '../../../modules/local/rmats/summarize/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN rMATS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RMATS {

    take:
    ch_bam
    ch_gtf
    ch_contrast

    rmats_read_len
    rmats_splice_diff_cutoff
    rmats_novel_splice_site
    rmats_min_intron_len
    rmats_max_exon_len
    rmats_paired_stats

    main:

    /*
     * Collect STAR BAMs and separate them according to the
     * control and treatment conditions in ch_contrast.
     */
    ch_rmats_groups = ch_bam
        .collect(flat: false)
        .combine(ch_contrast)
        .map { combined ->

            def contrast = combined[-1]
            def bams     = combined[0..-2]

            /*
            * rMATS calculates:
            *
            * IncLevelDifference = IncLevel1 - IncLevel2
            *
            * Therefore:
            * group1 / b1 = treatment/target
            * group2 / b2 = control/reference
            */
            def cond1 = contrast.treatment
            def cond2 = contrast.control

            def group1 = bams.findAll { it[0].condition == cond1 }
            def group2 = bams.findAll { it[0].condition == cond2 }

            if (!group1) {
                error "No BAM files found for rMATS treatment condition: ${cond1}"
            }

            if (!group2) {
                error "No BAM files found for rMATS control condition: ${cond2}"
            }

            tuple(
                [id: "${cond1}_vs_${cond2}"],
                cond1,
                group1.collect { it[0] },
                group1.collect { it[1] },
                cond2,
                group2.collect { it[0] },
                group2.collect { it[1] }
            )
        }


    /*
     * Create comma-separated b1.txt and b2.txt files.
     */
    CREATE_BAMLIST(
        ch_rmats_groups.map {
            contrast,
            cond1,
            meta1,
            bam1,
            cond2,
            meta2,
            bam2 ->

            tuple(
                contrast,
                cond1,
                bam1,
                cond2,
                bam2
            )
        }
    )


    /*
     * Reattach BAMs and metadata to the generated BAM-list files.
     */
    ch_bam_lists = CREATE_BAMLIST.out.bam_list1
        .join(
            CREATE_BAMLIST.out.bam_list2,
            by: 0
        )


    ch_rmats_input = ch_rmats_groups
        .join(
            ch_bam_lists,
            by: 0
        )
        .map {
            contrast,
            cond1,
            meta1,
            bam1,
            cond2,
            meta2,
            bam2,
            bam1_text,
            bam2_text ->

            tuple(
                contrast,
                cond1,
                meta1,
                bam1,
                bam1_text,
                cond2,
                meta2,
                bam2,
                bam2_text
            )
        }


    /*
    * Extract the GTF path from the nf-core-style metadata tuple.
    */
    ch_rmats_gtf = ch_gtf.map { _meta, gtf ->
        gtf
    }

    /*
    * Run the rMATS preparation stage.
    */
    RMATS_RUN(
        ch_rmats_gtf,
        ch_rmats_input,
        rmats_read_len,
        rmats_splice_diff_cutoff,
        rmats_novel_splice_site,
        rmats_min_intron_len,
        rmats_max_exon_len
    )


    /*
     * Attach the temporary rMATS output to the comparison input.
     */
    ch_rmats_summary_input = ch_rmats_input
        .join(
            RMATS_RUN.out.rmats_temp,
            by: 0
        )
        .map {
            contrast,
            cond1,
            meta1,
            bam1,
            bam1_text,
            cond2,
            meta2,
            bam2,
            bam2_text,
            rmats_temp ->

            tuple(
                contrast,
                cond1,
                meta1,
                bam1,
                bam1_text,
                cond2,
                meta2,
                bam2,
                bam2_text,
                rmats_temp
            )
        }


    /*
     * Run rMATS statistics and generate final event tables.
     */
    RMATS_SUMMARIZE(
        ch_rmats_gtf,
        ch_rmats_summary_input,
        rmats_read_len,
        rmats_splice_diff_cutoff,
        rmats_novel_splice_site,
        rmats_min_intron_len,
        rmats_max_exon_len,
        rmats_paired_stats
    )


    emit:
    results  = RMATS_SUMMARIZE.out.rmats_post
    gene_summary = RMATS_SUMMARIZE.out.gene_summary
    post_log = RMATS_SUMMARIZE.out.log
    prep_log = RMATS_RUN.out.log
}