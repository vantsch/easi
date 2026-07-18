/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAJIQ_BUILDGFF3                         } from '../../../modules/local/majiq/buildgff3/main'
include { MAJIQ_BUILDSJ                           } from '../../../modules/local/majiq/buildsj/main'
include { MAJIQ_BUILDUPDATE                       } from '../../../modules/local/majiq/buildupdate/main'
include { MAJIQ_PSICOVERAGE                       } from '../../../modules/local/majiq/psicoverage/main'
include { MAJIQ_QUANTIFY                          } from '../../../modules/local/majiq/quantify/main'
include { MAJIQ_DELTAPSI                          } from '../../../modules/local/majiq/deltapsi/main'
include { MAJIQ_HETEROGEN                         } from '../../../modules/local/majiq/heterogen/main'
include { MAJIQ_SGCOVERAGE                        } from '../../../modules/local/majiq/sgcoverage/main'
include { MAJIQ_MODULIZE as QUANTIFY_MODULIZE     } from '../../../modules/local/majiq/modulize/main'
include { MAJIQ_MODULIZE as DELTAPSI_MODULIZE     } from '../../../modules/local/majiq/modulize/main'
include { MAJIQ_MODULIZE as HETEROGEN_MODULIZE    } from '../../../modules/local/majiq/modulize/main'
include { MAJIQ_SUMMARIZE                         } from '../../../modules/local/majiq/downstream/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAJIQ SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAJIQ {

    take:
    ch_gff3
    ch_bam
    ch_contrast

    main:

    /*
     * Build the MAJIQ reference splice graph
     */
    MAJIQ_BUILDGFF3(
        ch_gff3
    )

    ch_majiq_splicegraph = MAJIQ_BUILDGFF3.out.splicegraph

    /*
     * Build per-sample MAJIQ splice-junction files
     */
    ch_bam_for_majiq = ch_bam
        .map { meta, bam ->
            tuple(meta, [bam])
        }

    ch_majiq_buildsj_input = ch_bam_for_majiq
        .combine(ch_majiq_splicegraph)

    MAJIQ_BUILDSJ(
        ch_majiq_buildsj_input
    )

    /*
     * Update the reference splice graph using detected junctions
     */
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

    ch_majiq_finished_splicegraph = MAJIQ_BUILDUPDATE.out.splicegraph

    /*
     * Calculate per-condition PSI coverage
     */
    ch_majiq_condition_samples_sj = MAJIQ_BUILDSJ.out.sj
        .map { meta, sj_file ->
            tuple(meta.condition, sj_file)
        }
        .groupTuple()
        .combine(ch_majiq_finished_splicegraph)

    MAJIQ_PSICOVERAGE(
        ch_majiq_condition_samples_sj
    )

    /*
     * Quantify PSI
     */
    MAJIQ_QUANTIFY(
        MAJIQ_PSICOVERAGE.out.psi_coverage
            .combine(ch_majiq_finished_splicegraph)
    )

    /*
     * Construct comparison input for deltapsi and heterogen
     */
    ch_majiq_contrast_comparison = ch_contrast
        .map { contrast ->
            [contrast.treatment, contrast]
        }
        .combine(MAJIQ_PSICOVERAGE.out.psi_coverage, by: 0)
        .map { values ->
            values[1] + [psicov1: values[2]]
        }
        .map { contrast ->
            [contrast.control, contrast]
        }
        .combine(MAJIQ_PSICOVERAGE.out.psi_coverage, by: 0)
        .map { values ->
            values[1] + [psicov2: values[2]]
        }

    ch_majiq_contrast_input = ch_majiq_contrast_comparison
        .combine(ch_majiq_finished_splicegraph)
        .map { values ->
            [
                values[0].contrast,
                values[0].treatment,
                values[0].control,
                values[0].psicov1,
                values[0].psicov2,
                values[2]
            ]
        }

    MAJIQ_DELTAPSI(
        ch_majiq_contrast_input
    )

    MAJIQ_HETEROGEN(
        ch_majiq_contrast_input
    )

    /*
     * Generate splice-graph coverage files
     */
    MAJIQ_SGCOVERAGE(
        ch_majiq_condition_samples_sj
    )

    /*
     * Modulize quantified PSI results
     */
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

    /*
     * Modulize deltapsi results
     */
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

    /*
     * Modulize heterogen results
     */
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

    /*
     * Summaize MAJIQ results into a single Excel workbook
     */
    MAJIQ_SUMMARIZE(DELTAPSI_MODULIZE.out.tsv_files,
        file("${projectDir}/modules/local/majiq/downstream/bin/report_events.py")
    )

    emit:
    splicegraph          = ch_majiq_finished_splicegraph
    sj                   = MAJIQ_BUILDSJ.out.sj
    psi_coverage         = MAJIQ_PSICOVERAGE.out.psi_coverage
    deltapsi             = MAJIQ_DELTAPSI.out.dpsicov
    heterogen            = MAJIQ_HETEROGEN.out.hetcov
    sgc                  = MAJIQ_SGCOVERAGE.out.sgc_files
    quantify_modulized   = QUANTIFY_MODULIZE.out.modulize
    deltapsi_modulized   = DELTAPSI_MODULIZE.out.modulize
    heterogen_modulized  = HETEROGEN_MODULIZE.out.modulize
    deltapsi_report      = MAJIQ_SUMMARIZE.out.workbook
}