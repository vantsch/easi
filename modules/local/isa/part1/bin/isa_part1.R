#!/usr/bin/env Rscript

# IsoformSwitchAnalyzeR Part 1
#
# Imports Salmon transcript quantification, identifies isoform switches,
# predicts ORFs, and exports nucleotide and amino-acid sequences.
#
# Pfam and SignalP should be run on the exported amino-acid FASTA before
# ISA Part 2 is executed.

# -------------------------------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) < 6) {
    stop(
        paste(
            "Expected arguments:",
            "<gtf>",
            "<transcript_fasta>",
            "<samplesheet>",
            "<contrastsheet>",
            "<alpha>",
            "<dIF_cutoff>"
        )
    )
}

# Salmon result directories are staged as direct subdirectories of the
# Nextflow process working directory.
salmon_output <- getwd()

gtf                  <- argv[1]
transcript_sequences <- argv[2]
samplesheet          <- argv[3]
contrastsheet        <- argv[4]
alpha                <- as.numeric(argv[5])
dIFcutoff            <- as.numeric(argv[6])


# -------------------------------------------------------------------------
# Validate inputs
# -------------------------------------------------------------------------

if (!file.exists(gtf)) {
    stop("GTF file does not exist: ", gtf)
}

if (!file.exists(transcript_sequences)) {
    stop(
        "Transcript FASTA does not exist: ",
        transcript_sequences
    )
}

if (!file.exists(samplesheet)) {
    stop("Samplesheet does not exist: ", samplesheet)
}

if (!file.exists(contrastsheet)) {
    stop("Contrastsheet does not exist: ", contrastsheet)
}

if (is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
}

if (is.na(dIFcutoff) || dIFcutoff < 0 || dIFcutoff > 1) {
    stop("dIFcutoff must be between 0 and 1.")
}


# -------------------------------------------------------------------------
# Load package
# -------------------------------------------------------------------------

suppressPackageStartupMessages(
    library(IsoformSwitchAnalyzeR)
)


# -------------------------------------------------------------------------
# ISA Part 1
# -------------------------------------------------------------------------

isoformSwitchAnalysisPart1 <- function(
    switchAnalyzeRlist,
    pathToGTF = NULL,
    pathToOutput = NULL,
    alpha = NULL,
    dIFcutoff = NULL
) {
    isConditional <-
        switchAnalyzeRlist$sourceId != "preDefinedSwitches"

    # ---------------------------------------------------------------------
    # Filter low-information transcripts and genes
    # ---------------------------------------------------------------------

    if (isConditional) {
        switchAnalyzeRlist <- preFilter(
            switchAnalyzeRlist = switchAnalyzeRlist,
            removeSingleIsoformGenes = TRUE,
            quiet = FALSE,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )
    }

    # ---------------------------------------------------------------------
    # Differential isoform usage testing
    # ---------------------------------------------------------------------

    if (isConditional) {
        replicate_counts <-
            switchAnalyzeRlist$conditions$nrReplicates

        message(
            "Replicate counts detected by ISA: ",
            paste(replicate_counts, collapse = ", ")
        )

        if (any(replicate_counts > 5)) {
            message("Running isoformSwitchTestSatuRn()...")

            switchAnalyzeRlist <- isoformSwitchTestSatuRn(
                switchAnalyzeRlist = switchAnalyzeRlist,
                reduceToSwitchingGenes = TRUE,
                quiet = FALSE,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

        } else {
            message("Running isoformSwitchTestDEXSeq()...")

            switchAnalyzeRlist <- isoformSwitchTestDEXSeq(
                switchAnalyzeRlist = switchAnalyzeRlist,
                reduceToSwitchingGenes = FALSE,
                quiet = FALSE,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )
        }

        if (
            is.null(switchAnalyzeRlist$isoformSwitchAnalysis) ||
            nrow(switchAnalyzeRlist$isoformSwitchAnalysis) == 0
        ) {
            stop(
                "No isoform switches were identified ",
                "with the current cutoffs."
            )
        }
    }

    # ---------------------------------------------------------------------
    # Import annotated ORFs and predict missing ORFs
    # ---------------------------------------------------------------------

    if (is.null(switchAnalyzeRlist$orfAnalysis)) {
        message("Adding annotated ORFs from the GTF...")

        switchAnalyzeRlist <- addORFfromGTF(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToGTF = pathToGTF,
            quiet = FALSE
        )

        missing_orfs <- (
            !is.null(switchAnalyzeRlist$orfAnalysis) &&
            "orf_origin" %in%
                colnames(switchAnalyzeRlist$orfAnalysis) &&
            any(
                switchAnalyzeRlist$orfAnalysis$orf_origin ==
                    "not_annotated_yet",
                na.rm = TRUE
            )
        )

        if (missing_orfs) {
            message(
                "Predicting ORFs for isoforms without ",
                "annotated ORFs..."
            )

            switchAnalyzeRlist <- analyzeNovelIsoformORF(
                switchAnalyzeRlist = switchAnalyzeRlist,
                analysisAllIsoformsWithoutORF = TRUE
            )
        }
    }

    # ---------------------------------------------------------------------
    # Export nucleotide and amino-acid sequences
    # ---------------------------------------------------------------------

    message(
        "Extracting nucleotide and amino-acid sequences ",
        "for Pfam and SignalP..."
    )

    switchAnalyzeRlist <- extractSequence(
        switchAnalyzeRlist = switchAnalyzeRlist,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        addToSwitchAnalyzeRlist = TRUE,
        writeToFile = TRUE,
        pathToOutput = pathToOutput,
        quiet = FALSE,
        alpha = alpha,
        dIFcutoff = dIFcutoff
    )

    return(switchAnalyzeRlist)
}


# -------------------------------------------------------------------------
# Import Salmon expression
# -------------------------------------------------------------------------

message("Importing Salmon quantification results...")

salmonQuant <- importIsoformExpression(
    parentDir = salmon_output,
    showProgress = TRUE,
    quiet = FALSE
)


# -------------------------------------------------------------------------
# Build design matrix
# -------------------------------------------------------------------------

design <- read.csv(
    samplesheet,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_design_columns <- c(
    "sample",
    "condition"
)

missing_design_columns <- setdiff(
    required_design_columns,
    colnames(design)
)

if (length(missing_design_columns) > 0) {
    stop(
        "Samplesheet is missing required columns: ",
        paste(missing_design_columns, collapse = ", ")
    )
}

design <- design[
    ,
    required_design_columns,
    drop = FALSE
]

colnames(design) <- c(
    "sampleID",
    "condition"
)

design$sampleID <- as.character(design$sampleID)
design$condition <- as.character(design$condition)


# -------------------------------------------------------------------------
# Build comparison table
# -------------------------------------------------------------------------

contrasts <- read.csv(
    contrastsheet,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_contrast_columns <- c(
    "treatment",
    "control"
)

missing_contrast_columns <- setdiff(
    required_contrast_columns,
    colnames(contrasts)
)

if (length(missing_contrast_columns) > 0) {
    stop(
        "Contrastsheet is missing required columns: ",
        paste(missing_contrast_columns, collapse = ", ")
    )
}

contrasts <- contrasts[
    ,
    required_contrast_columns,
    drop = FALSE
]

if (nrow(contrasts) != 1) {
    stop(
        "ISA Part 1 currently expects exactly one comparison, but found ",
        nrow(contrasts),
        "."
    )
}

# EASI convention:
#
# treatment = target
# control   = reference
#
# IsoformSwitchAnalyzeR calculates:
#
# dIF = condition_2 - condition_1
#
# Therefore:
#
# condition_1 = reference/control
# condition_2 = target/treatment

comparisons <- data.frame(
    condition_1 = as.character(contrasts$control),
    condition_2 = as.character(contrasts$treatment),
    stringsAsFactors = FALSE
)

reference_condition <- comparisons$condition_1[1]
target_condition    <- comparisons$condition_2[1]

if (
    is.na(reference_condition) ||
    is.na(target_condition) ||
    reference_condition == "" ||
    target_condition == ""
) {
    stop(
        "The treatment and control conditions must not be empty."
    )
}

if (reference_condition == target_condition) {
    stop(
        "The treatment and control conditions must be different."
    )
}

available_conditions <- unique(design$condition)

missing_conditions <- setdiff(
    c(reference_condition, target_condition),
    available_conditions
)

if (length(missing_conditions) > 0) {
    stop(
        "Contrast conditions missing from the samplesheet: ",
        paste(missing_conditions, collapse = ", ")
    )
}

unexpected_conditions <- setdiff(
    available_conditions,
    c(reference_condition, target_condition)
)

if (length(unexpected_conditions) > 0) {
    stop(
        "The samplesheet contains conditions not included in the ",
        "requested ISA contrast: ",
        paste(unexpected_conditions, collapse = ", ")
    )
}

design$condition <- factor(
    design$condition,
    levels = c(
        reference_condition,
        target_condition
    )
)

message("ISA design matrix:")
print(design)

message("ISA comparisons:")
print(comparisons)

message(
    "ISA contrast orientation: ",
    target_condition,
    " vs ",
    reference_condition,
    "; expected dIF = ",
    target_condition,
    " - ",
    reference_condition
)


# -------------------------------------------------------------------------
# Import expression and annotation into ISA
# -------------------------------------------------------------------------

message("Building SwitchAnalyzeRlist...")

SwitchList <- importRdata(
    isoformCountMatrix = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix = design,
    isoformExonAnnoation = gtf,
    isoformNtFasta = transcript_sequences,
    comparisonsToMake = comparisons,

    # Salmon and the generated transcript FASTA contain unversioned
    # ENST identifiers, while the GTF importer may construct identifiers
    # such as ENST00000479818.5.
    ignoreAfterPeriod = TRUE,

    showProgress = TRUE,
    quiet = FALSE
)

message("SwitchAnalyzeRlist imported successfully.")


# -------------------------------------------------------------------------
# Validate imported contrast orientation
# -------------------------------------------------------------------------

if (is.null(SwitchList$isoformFeatures)) {
    stop(
        "SwitchAnalyzeRlist does not contain isoformFeatures after import."
    )
}

required_import_columns <- c(
    "condition_1",
    "condition_2"
)

missing_import_columns <- setdiff(
    required_import_columns,
    colnames(SwitchList$isoformFeatures)
)

if (length(missing_import_columns) > 0) {
    stop(
        "Cannot validate ISA contrast orientation after importRdata(). ",
        "Missing columns: ",
        paste(missing_import_columns, collapse = ", ")
    )
}

observed_comparisons <- unique(
    SwitchList$isoformFeatures[
        ,
        required_import_columns,
        drop = FALSE
    ]
)

observed_comparisons <- observed_comparisons[
    complete.cases(observed_comparisons),
    ,
    drop = FALSE
]

if (
    nrow(observed_comparisons) != 1 ||
    observed_comparisons$condition_1[1] != reference_condition ||
    observed_comparisons$condition_2[1] != target_condition
) {
    observed_text <- if (nrow(observed_comparisons) == 0) {
        "no complete comparison was found"
    } else {
        paste(
            apply(
                observed_comparisons,
                1,
                paste,
                collapse = " -> "
            ),
            collapse = "; "
        )
    }

    stop(
        "ISA contrast orientation mismatch after importRdata(). ",
        "Expected ",
        reference_condition,
        " -> ",
        target_condition,
        "; observed ",
        observed_text,
        "."
    )
}

message(
    "Imported ISA comparison validated: condition_1 = ",
    reference_condition,
    ", condition_2 = ",
    target_condition
)


# -------------------------------------------------------------------------
# Run ISA Part 1
# -------------------------------------------------------------------------

dir.create(
    "results",
    recursive = TRUE,
    showWarnings = FALSE
)

SwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist = SwitchList,
    pathToGTF = gtf,
    pathToOutput = "results/",
    alpha = alpha,
    dIFcutoff = dIFcutoff
)


# -------------------------------------------------------------------------
# Validate final dIF direction
# -------------------------------------------------------------------------

if (is.null(SwitchList$isoformFeatures)) {
    stop(
        "SwitchAnalyzeRlist does not contain isoformFeatures after ",
        "Part 1 analysis."
    )
}

features <- SwitchList$isoformFeatures

required_direction_columns <- c(
    "condition_1",
    "condition_2",
    "IF1",
    "IF2",
    "dIF"
)

missing_direction_columns <- setdiff(
    required_direction_columns,
    colnames(features)
)

if (length(missing_direction_columns) > 0) {
    stop(
        "Cannot validate ISA dIF direction. Missing columns: ",
        paste(missing_direction_columns, collapse = ", ")
    )
}

final_comparisons <- unique(
    features[
        ,
        c("condition_1", "condition_2"),
        drop = FALSE
    ]
)

final_comparisons <- final_comparisons[
    complete.cases(final_comparisons),
    ,
    drop = FALSE
]

if (
    nrow(final_comparisons) != 1 ||
    final_comparisons$condition_1[1] != reference_condition ||
    final_comparisons$condition_2[1] != target_condition
) {
    stop(
        "ISA contrast orientation changed during Part 1 analysis."
    )
}

direction_check <- with(
    features,
    is.na(dIF) |
        is.na(IF1) |
        is.na(IF2) |
        abs(dIF - (IF2 - IF1)) < 1e-8
)

if (!all(direction_check)) {
    failed_rows <- sum(!direction_check)

    stop(
        "ISA dIF validation failed for ",
        failed_rows,
        " isoform rows: dIF is not equal to IF2 - IF1."
    )
}

message(
    "ISA contrast validated successfully: dIF = ",
    target_condition,
    " - ",
    reference_condition
)


# -------------------------------------------------------------------------
# Save intermediate ISA object
# -------------------------------------------------------------------------

saveRDS(
    SwitchList,
    file = "switchlist_part1.rds"
)


# -------------------------------------------------------------------------
# Save preliminary tables
# -------------------------------------------------------------------------

write.table(
    SwitchList$isoformSwitchAnalysis,
    file = "isa_part1_switches.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
)

write.table(
    SwitchList$isoformFeatures,
    file = "isa_part1_isoform_features.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
)

if (!is.null(SwitchList$orfAnalysis)) {
    write.table(
        SwitchList$orfAnalysis,
        file = "isa_part1_orf_analysis.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
    )
}


# -------------------------------------------------------------------------
# List sequence files produced by extractSequence()
# -------------------------------------------------------------------------

sequence_files <- list.files(
    path = "results",
    recursive = TRUE,
    full.names = TRUE
)

writeLines(
    sequence_files,
    con = "isa_part1_sequence_files.txt"
)

message("Sequence files generated by ISA:")

if (length(sequence_files) > 0) {
    message(
        paste(
            sequence_files,
            collapse = "\n"
        )
    )
} else {
    warning(
        "No sequence files were detected under results/. ",
        "Check the output from extractSequence()."
    )
}


# -------------------------------------------------------------------------
# Completion and session information
# -------------------------------------------------------------------------

message("ISA Part 1 completed successfully.")

citation("IsoformSwitchAnalyzeR")
sessionInfo()