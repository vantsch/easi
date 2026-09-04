#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(
        paste0(
            "Usage: isa_part2.R ",
            "<switchlist_part1.rds> <pfam_isa.tsv> ",
            "<isa_annotation.gtf> ",
            "[alpha] [dIF_cutoff]"
        ),
        call. = FALSE
    )
}

switchlist_file <- args[1]
pfam_file       <- args[2]
annotation_gtf  <- args[3]

alpha <- if (length(args) >= 4) {
    as.numeric(args[4])
} else {
    0.05
}

dIF_cutoff <- if (length(args) >= 5) {
    as.numeric(args[5])
} else {
    0.10
}

if (
    length(alpha) != 1 ||
    is.na(alpha) ||
    alpha <= 0 ||
    alpha > 1
) {
    stop(
        "alpha must be greater than 0 and no greater than 1.",
        call. = FALSE
    )
}

if (
    length(dIF_cutoff) != 1 ||
    is.na(dIF_cutoff) ||
    dIF_cutoff < 0
) {
    stop(
        "dIF cutoff must be zero or greater.",
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
})

read_transcript_gene_mapping <- function(
    gtf_file,
    chunk_size = 100000L
) {

    message(
        "Reading transcript-to-gene mapping from: ",
        gtf_file
    )

    if (
        length(gtf_file) != 1 ||
        is.na(gtf_file) ||
        gtf_file == "" ||
        !file.exists(gtf_file)
    ) {
        stop(
            paste0(
                "GTF file does not exist: ",
                gtf_file
            ),
            call. = FALSE
        )
    }

    extract_attribute <- function(
        attributes,
        attribute_name
    ) {

        pattern <- paste0(
            ".*(?:^|;[[:space:]]*)",
            attribute_name,
            "[[:space:]]+\"([^\"]+)\".*"
        )

        found <- grepl(
            pattern,
            attributes,
            perl = TRUE
        )

        result <- rep(
            NA_character_,
            length(attributes)
        )

        result[found] <- sub(
            pattern,
            "\\1",
            attributes[found],
            perl = TRUE
        )

        result
    }

    connection <- if (
        grepl(
            "\\.gz$",
            gtf_file,
            ignore.case = TRUE
        )
    ) {
        gzfile(
            gtf_file,
            open = "rt"
        )
    } else {
        file(
            gtf_file,
            open = "rt"
        )
    }

    on.exit(
        close(connection),
        add = TRUE
    )

    mapping_chunks <- list()
    chunk_index <- 0L

    repeat {

        lines <- readLines(
            connection,
            n = chunk_size,
            warn = FALSE
        )

        if (
            length(lines) == 0
        ) {
            break
        }

        # Remove comments and empty lines.
        lines <- lines[
            nzchar(lines) &
            !startsWith(lines, "#")
        ]

        if (
            length(lines) == 0
        ) {
            next
        }

        fields <- strsplit(
            lines,
            split = "\t",
            fixed = TRUE
        )

        valid_gtf_rows <- lengths(fields) >= 9L

        if (
            !any(valid_gtf_rows)
        ) {
            next
        }

        fields <- fields[
            valid_gtf_rows
        ]

        feature <- vapply(
            fields,
            function(x) {
                x[[3]]
            },
            character(1)
        )

        transcript_rows <- feature == "transcript"

        if (
            !any(transcript_rows)
        ) {
            next
        }

        transcript_fields <- fields[
            transcript_rows
        ]

        attributes <- vapply(
            transcript_fields,
            function(x) {
                x[[9]]
            },
            character(1)
        )

        transcript_id <- extract_attribute(
            attributes,
            "transcript_id"
        )

        gene_id <- extract_attribute(
            attributes,
            "gene_id"
        )

        gene_name <- extract_attribute(
            attributes,
            "gene_name"
        )

        valid_mapping <-
            !is.na(transcript_id) &
            transcript_id != "" &
            !is.na(gene_id) &
            gene_id != ""

        if (
            !any(valid_mapping)
        ) {
            next
        }

        chunk_index <- chunk_index + 1L

        mapping_chunks[[chunk_index]] <- data.frame(
            isoform_id = transcript_id[
                valid_mapping
            ],
            gene_id = gene_id[
                valid_mapping
            ],
            gene_name = gene_name[
                valid_mapping
            ],
            stringsAsFactors = FALSE
        )
    }

    if (
        length(mapping_chunks) == 0
    ) {
        stop(
            paste0(
                "No transcript mappings were extracted from ",
                gtf_file,
                ". Confirm that the annotation contains ",
                "transcript features with transcript_id and gene_id."
            ),
            call. = FALSE
        )
    }

    transcript_mapping <- do.call(
        rbind,
        mapping_chunks
    )

    rm(mapping_chunks)

    transcript_mapping$isoform_id <- trimws(
        transcript_mapping$isoform_id
    )

    transcript_mapping$gene_id <- trimws(
        transcript_mapping$gene_id
    )

    transcript_mapping$gene_name <- trimws(
        transcript_mapping$gene_name
    )

    # Remove Ensembl version suffixes.
    transcript_mapping$isoform_id <- sub(
        "\\.[0-9]+$",
        "",
        transcript_mapping$isoform_id
    )

    transcript_mapping$gene_id <- sub(
        "\\.[0-9]+$",
        "",
        transcript_mapping$gene_id
    )

    valid_rows <-
        !is.na(transcript_mapping$isoform_id) &
        transcript_mapping$isoform_id != "" &
        !is.na(transcript_mapping$gene_id) &
        transcript_mapping$gene_id != ""

    transcript_mapping <- transcript_mapping[
        valid_rows,
        ,
        drop = FALSE
    ]

    # Detect genuine cases where one transcript maps to
    # more than one gene.
    transcript_gene_pairs <- unique(
        transcript_mapping[
            c(
                "isoform_id",
                "gene_id"
            )
        ]
    )

    duplicated_transcript <-
        duplicated(
            transcript_gene_pairs$isoform_id
        ) |
        duplicated(
            transcript_gene_pairs$isoform_id,
            fromLast = TRUE
        )

    if (
        any(duplicated_transcript)
    ) {

        conflicting_transcripts <- unique(
            transcript_gene_pairs$isoform_id[
                duplicated_transcript
            ]
        )

        preview <- paste(
            head(
                conflicting_transcripts,
                10
            ),
            collapse = ", "
        )

        stop(
            paste0(
                length(conflicting_transcripts),
                " transcript IDs map to more than one gene ID. ",
                "Examples: ",
                preview
            ),
            call. = FALSE
        )
    }

    # Keep one annotation row per transcript.
    transcript_mapping <- transcript_mapping[
        !duplicated(
            transcript_mapping$isoform_id
        ),
        ,
        drop = FALSE
    ]

    rownames(transcript_mapping) <- NULL

    message(
        "Loaded mappings for ",
        format(
            nrow(transcript_mapping),
            big.mark = ","
        ),
        " transcripts and ",
        format(
            length(
                unique(
                    transcript_mapping$gene_id
                )
            ),
            big.mark = ","
        ),
        " genes."
    )

    transcript_mapping
}

message("Reading ISA Part 1 object: ", switchlist_file)

SwitchList <- readRDS(
    switchlist_file
)

message("Importing Pfam annotations: ", pfam_file)

SwitchList <- analyzePFAM(
    switchAnalyzeRlist   = SwitchList,
    pathToPFAMresultFile = pfam_file,
    showProgress         = FALSE,
    quiet                = TRUE
)

message("Analyzing alternative splicing")

SwitchList <- analyzeAlternativeSplicing(
    switchAnalyzeRlist = SwitchList,
    quiet              = TRUE
)

message("Analyzing switch consequences")

SwitchList <- analyzeSwitchConsequences(
    switchAnalyzeRlist = SwitchList,
    consequencesToAnalyze = c(
        "intron_retention",
        "ORF_seq_similarity",
        "NMD_status",
        "domains_identified"
    ),
    quiet = TRUE
)

############################################################
## Standardised gene-level summary                         ##
############################################################

message("Generating standardised ISA gene-level summary")

isoform_features <- SwitchList$isoformFeatures

############################################################
## Restore Ensembl gene IDs from the annotation GTF        ##
############################################################

transcript_mapping <- read_transcript_gene_mapping(
    annotation_gtf
)

isoform_features$isoform_id <- sub(
    "\\.[0-9]+$",
    "",
    as.character(
        isoform_features$isoform_id
    )
)

mapping_index <- match(
    isoform_features$isoform_id,
    transcript_mapping$isoform_id
)

mapped_gene_id <- transcript_mapping$gene_id[
    mapping_index
]

mapped_gene_name <- transcript_mapping$gene_name[
    mapping_index
]

mapped <- (
    !is.na(mapped_gene_id) &
    mapped_gene_id != ""
)

message(
    "Mapped ",
    sum(mapped),
    " of ",
    nrow(isoform_features),
    " ISA isoform rows to Ensembl gene IDs."
)

if (!all(mapped)) {

    unmapped_isoforms <- unique(
        isoform_features$isoform_id[
            !mapped
        ]
    )

    preview <- paste(
        head(
            unmapped_isoforms,
            10
        ),
        collapse = ", "
    )

    stop(
        paste0(
            length(unmapped_isoforms),
            " unique ISA transcript IDs could not be found in ",
            "the supplied annotation GTF. Examples: ",
            preview
        ),
        call. = FALSE
    )
}

# Retain the original ISA gene field for traceability.
isoform_features$original_gene_id <- as.character(
    isoform_features$gene_id
)

# Standard EASI identifiers:
#   gene_id   = Ensembl gene ID
#   gene_name = gene symbol
isoform_features$gene_id <- mapped_gene_id

has_mapped_gene_name <- (
    !is.na(mapped_gene_name) &
    mapped_gene_name != ""
)

isoform_features$gene_name[
    has_mapped_gene_name
] <- mapped_gene_name[
    has_mapped_gene_name
]

SwitchList$isoformFeatures <- isoform_features

required_columns <- c(
    "gene_id",
    "gene_name",
    "isoform_id",
    "dIF",
    "isoform_switch_q_value",
    "gene_switch_q_value"
)

missing_columns <- setdiff(
    required_columns,
    colnames(isoform_features)
)

if (length(missing_columns) > 0) {
    stop(
        paste0(
            "SwitchList$isoformFeatures is missing required columns: ",
            paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
    )
}

isoform_features$gene_id <- as.character(
    isoform_features$gene_id
)

isoform_features$gene_name <- as.character(
    isoform_features$gene_name
)

isoform_features$isoform_id <- as.character(
    isoform_features$isoform_id
)

isoform_features$dIF <- suppressWarnings(
    as.numeric(
        isoform_features$dIF
    )
)

isoform_features$isoform_switch_q_value <- suppressWarnings(
    as.numeric(
        isoform_features$isoform_switch_q_value
    )
)

isoform_features$gene_switch_q_value <- suppressWarnings(
    as.numeric(
        isoform_features$gene_switch_q_value
    )
)

valid_gene_rows <-
    !is.na(isoform_features$gene_id) &
    isoform_features$gene_id != ""

isoform_features <- isoform_features[
    valid_gene_rows,
    ,
    drop = FALSE
]

# A gene-level ISA test is considered detected when ISA produced
# a gene-level switch q-value.
isoform_features$detected_event <-
    !is.na(isoform_features$gene_switch_q_value)

# An isoform contributes to a significant gene-level switch when:
#   1. the gene-level switch q-value passes alpha; and
#   2. the isoform-level absolute dIF passes the requested cutoff.
isoform_features$significant_event <-
    isoform_features$detected_event &
    isoform_features$gene_switch_q_value < alpha &
    !is.na(isoform_features$dIF) &
    abs(isoform_features$dIF) >= dIF_cutoff

split_by_gene <- split(
    isoform_features,
    isoform_features$gene_id,
    drop = TRUE
)

gene_summary_rows <- lapply(
    split_by_gene,
    function(gene_data) {

        available_names <- unique(
            gene_data$gene_name[
                !is.na(gene_data$gene_name) &
                gene_data$gene_name != ""
            ]
        )

        gene_name <- if (length(available_names) > 0) {
            available_names[[1]]
        } else {
            NA_character_
        }

        detected <- any(
            gene_data$detected_event
        )

        significant <- any(
            gene_data$significant_event
        )

        detected_effects <- gene_data$dIF[
            gene_data$detected_event &
                !is.na(gene_data$dIF)
        ]

        significant_isoforms <- gene_data[
            gene_data$significant_event &
                !is.na(gene_data$dIF),
            ,
            drop = FALSE
        ]

        n_significant_isoforms <- length(
            unique(
                significant_isoforms$isoform_id
            )
        )

        has_positive <- any(
            significant_isoforms$dIF > 0,
            na.rm = TRUE
        )

        has_negative <- any(
            significant_isoforms$dIF < 0,
            na.rm = TRUE
        )

        effect_direction <- if (!significant) {

            "none"

        } else if (n_significant_isoforms == 1) {

            "single_isoform_change"

        } else if (
            n_significant_isoforms == 2 &&
            has_positive &&
            has_negative
        ) {

            "isoform_switch"

        } else {

            "complex_switch"
        }

        max_abs_effect <- if (
            length(detected_effects) == 0
        ) {
            NA_real_
        } else {
            round(
                max(
                    abs(detected_effects),
                    na.rm = TRUE
                ),
                3
            )
        }

        data.frame(
            gene_id = as.character(
                gene_data$gene_id[[1]]
            ),
            gene_name = gene_name,
            method = "ISA",
            detected = detected,
            significant = significant,

            # ISA produces one gene-level switch test per gene
            # for each comparison.
            n_events = as.integer(
                detected
            ),
            n_significant_events = as.integer(
                significant
            ),

            max_abs_effect = max_abs_effect,
            effect_direction = effect_direction,
            event_types = "isoform_switch",
            stringsAsFactors = FALSE
        )
    }
)

if (length(gene_summary_rows) > 0) {

    gene_summary <- do.call(
        rbind,
        gene_summary_rows
    )

    gene_summary <- gene_summary[
        order(
            -as.integer(gene_summary$significant),
            -gene_summary$n_significant_events,
            -gene_summary$max_abs_effect,
            gene_summary$gene_name,
            gene_summary$gene_id,
            na.last = TRUE
        ),
        ,
        drop = FALSE
    ]

    rownames(gene_summary) <- NULL

} else {

    gene_summary <- data.frame(
        gene_id = character(0),
        gene_name = character(0),
        method = character(0),
        detected = logical(0),
        significant = logical(0),
        n_events = integer(0),
        n_significant_events = integer(0),
        max_abs_effect = numeric(0),
        effect_direction = character(0),
        event_types = character(0),
        stringsAsFactors = FALSE
    )
}

############################################################
## Switch plots                                            ##
############################################################

message("Generating switch plots")

dir.create(
    "switch_plots",
    recursive = TRUE,
    showWarnings = FALSE
)

switchPlotTopSwitches(
    switchAnalyzeRlist = SwitchList,
    n = Inf,
    pathToOutput = "switch_plots",
    filterForConsequences = TRUE,
    splitFunctionalConsequences = FALSE,
    IFcutoff = 0,
    additionalArguments = list(
        plotTopology = FALSE
    ),
    quiet = TRUE
)

############################################################
## Consequence summary                                     ##
############################################################

message("Generating consequence summary")

pdf(
    file = "common_switch_consequences.pdf",
    width = 10,
    height = 7
)

print(
    extractConsequenceSummary(
        switchAnalyzeRlist = SwitchList,
        plot = TRUE,
        returnResult = FALSE
    )
)

dev.off()

############################################################
## Write outputs                                           ##
############################################################

message("Writing ISA Part 2 outputs")

saveRDS(
    SwitchList,
    "switchlist_part2.rds"
)

write.table(
    SwitchList$domainAnalysis,
    file = "isa_part2_domain_analysis.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    SwitchList$isoformFeatures,
    file = "isa_part2_isoform_features.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    gene_summary,
    file = "isa_part2.gene_summary.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

message("ISA Part 2 completed")