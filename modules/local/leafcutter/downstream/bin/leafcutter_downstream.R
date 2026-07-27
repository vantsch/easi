#!/usr/bin/env Rscript

############################################################
## LeafCutter downstream annotation and table generation   ##
############################################################

options(stringsAsFactors = FALSE)

############################################################
## General functions                                       ##
############################################################

`%||%` <- function(x, fallback) {
    if (
        is.null(x) ||
        length(x) == 0 ||
        is.na(x) ||
        identical(x, "")
    ) {
        fallback
    } else {
        x
    }
}


is_valid_string <- function(x) {
    !is.null(x) &&
        length(x) == 1 &&
        !is.na(x) &&
        nzchar(trimws(x))
}


print_usage <- function() {
    cat(
        paste0(
            "\nLeafCutter downstream annotation\n\n",
            "Required arguments:\n",
            "  --cluster_file FILE       LeafCutter *_cluster_significance.txt\n",
            "  --effect_file FILE        LeafCutter *_effect_sizes.txt\n",
            "  --counts_file FILE        LeafCutter *_perind_numers.counts.gz\n",
            "  --groups_file FILE        LeafCutter groups_file.txt\n",
            "  --gtf_file FILE           GTF or GTF.GZ annotation\n",
            "  --output_prefix STRING    Output prefix\n\n",
            "Optional arguments:\n",
            "  --organism STRING         Organism label [default: human]\n",
            "  --alpha NUMBER            Adjusted p-value threshold [default: 0.05]\n",
            "  --min_abs_effect NUMBER   Minimum absolute delta PSI [default: 0.10]\n",
            "  --help                    Show this message\n\n"
        )
    )
}


parse_cli_args <- function(args) {

    if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
        print_usage()
        quit(status = 0)
    }

    parsed <- list()
    index <- 1L

    while (index <= length(args)) {

        argument <- args[[index]]

        if (!startsWith(argument, "--")) {
            stop(
                "Unexpected positional argument: ",
                argument,
                call. = FALSE
            )
        }

        key <- sub("^--", "", argument)

        if (index == length(args)) {
            stop(
                "No value supplied for --",
                key,
                call. = FALSE
            )
        }

        value <- args[[index + 1L]]

        if (startsWith(value, "--")) {
            stop(
                "No value supplied for --",
                key,
                call. = FALSE
            )
        }

        parsed[[key]] <- value
        index <- index + 2L
    }

    parsed
}


validate_file <- function(path, label) {

    if (!is_valid_string(path)) {
        stop(
            "Missing required argument: ",
            label,
            call. = FALSE
        )
    }

    if (!file.exists(path)) {
        stop(
            label,
            " does not exist: ",
            path,
            call. = FALSE
        )
    }
}


require_columns <- function(data, columns, table_name) {

    missing_columns <- setdiff(columns, colnames(data))

    if (length(missing_columns) > 0) {
        stop(
            table_name,
            " is missing required columns: ",
            paste(missing_columns, collapse = ", "),
            call. = FALSE
        )
    }
}


write_tsv <- function(data, filename) {

    data.table::fwrite(
        data,
        file = filename,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )
}


collapse_unique <- function(values) {

    values <- as.character(values)
    values <- values[!is.na(values) & values != ""]

    if (length(values) == 0) {
        return(NA_character_)
    }

    paste(unique(values), collapse = ";")
}


collapse_semicolon_values <- function(values) {

    values <- as.character(values)
    values <- values[!is.na(values) & values != ""]

    if (length(values) == 0) {
        return(NA_character_)
    }

    split_values <- unlist(
        strsplit(
            values,
            split = ";",
            fixed = TRUE
        ),
        use.names = FALSE
    )

    split_values <- split_values[
        !is.na(split_values) &
            split_values != ""
    ]

    if (length(split_values) == 0) {
        return(NA_character_)
    }

    paste(unique(split_values), collapse = ";")
}


safe_min <- function(x) {

    x <- suppressWarnings(as.integer(x))

    if (length(x) == 0 || all(is.na(x))) {
        return(NA_integer_)
    }

    min(x, na.rm = TRUE)
}


safe_max <- function(x) {

    x <- suppressWarnings(as.integer(x))

    if (length(x) == 0 || all(is.na(x))) {
        return(NA_integer_)
    }

    max(x, na.rm = TRUE)
}


read_leafcutter_table <- function(file) {

    connection <- if (grepl("\\.gz$", file, ignore.case = TRUE)) {
        gzfile(file, open = "rt")
    } else {
        file
    }

    on.exit({
        if (inherits(connection, "connection")) {
            close(connection)
        }
    }, add = TRUE)

    read.delim(
        connection,
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = ""
    )
}

############################################################
## GTF parsing                                              ##
############################################################

extract_gtf_attribute <- function(attributes, key) {

    pattern <- paste0(
        "(?:^|;[[:space:]]*)",
        key,
        "[[:space:]]+\"([^\"]+)\""
    )

    matches <- regexec(
        pattern,
        attributes,
        perl = TRUE
    )

    values <- regmatches(attributes, matches)

    vapply(
        values,
        function(value) {
            if (length(value) >= 2) {
                value[[2]]
            } else {
                NA_character_
            }
        },
        character(1)
    )
}


read_gtf_genes <- function(gtf_file) {

    message("Reading GTF annotation: ", gtf_file)

    if (grepl("\\.gz$", gtf_file, ignore.case = TRUE)) {
        gtf_command <- paste(
            "gzip -dc",
            shQuote(gtf_file)
        )
    } else {
        gtf_command <- paste(
            "cat",
            shQuote(gtf_file)
        )
    }

    gtf <- data.table::fread(
        cmd = gtf_command,
        sep = "\t",
        header = FALSE,
        quote = "",
        fill = TRUE,
        data.table = FALSE,
        showProgress = FALSE,
        select = seq_len(9)
    )

    if (ncol(gtf) < 9) {
        stop(
            "The GTF must contain nine tab-delimited columns.",
            call. = FALSE
        )
    }

    colnames(gtf) <- c(
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute"
    )

    # fread does not apply a comment character when reading through cmd.
    gtf <- gtf[
        !startsWith(gtf[["seqname"]], "#"),
        ,
        drop = FALSE
    ]

    gtf[["start"]] <- as.integer(gtf[["start"]])
    gtf[["end"]] <- as.integer(gtf[["end"]])

    gtf[["gene_id"]] <- extract_gtf_attribute(
        gtf[["attribute"]],
        "gene_id"
    )

    gtf[["gene_name"]] <- extract_gtf_attribute(
        gtf[["attribute"]],
        "gene_name"
    )

    gtf[["gene_biotype"]] <- extract_gtf_attribute(
        gtf[["attribute"]],
        "gene_biotype"
    )

    missing_biotype <- is.na(gtf[["gene_biotype"]])

    if (any(missing_biotype)) {
        gtf[["gene_biotype"]][missing_biotype] <-
            extract_gtf_attribute(
                gtf[["attribute"]][missing_biotype],
                "gene_type"
            )
    }

    gene_rows <- gtf[
        gtf[["feature"]] == "gene" &
            !is.na(gtf[["gene_id"]]),
        ,
        drop = FALSE
    ]

    if (nrow(gene_rows) > 0) {

        genes <- gene_rows[
            ,
            c(
                "seqname",
                "start",
                "end",
                "strand",
                "gene_id",
                "gene_name",
                "gene_biotype"
            ),
            drop = FALSE
        ]

    } else {

        warning(
            "No explicit gene features found in the GTF; ",
            "constructing gene ranges from exon features."
        )

        exon_rows <- gtf[
            gtf[["feature"]] == "exon" &
                !is.na(gtf[["gene_id"]]),
            ,
            drop = FALSE
        ]

        if (nrow(exon_rows) == 0) {
            stop(
                "No usable gene or exon features were found in the GTF.",
                call. = FALSE
            )
        }

        genes <- exon_rows |>
            dplyr::group_by(
                seqname,
                strand,
                gene_id,
                gene_name,
                gene_biotype
            ) |>
            dplyr::summarise(
                start = min(start, na.rm = TRUE),
                end = max(end, na.rm = TRUE),
                .groups = "drop"
            ) |>
            as.data.frame()

        genes <- genes[
            ,
            c(
                "seqname",
                "start",
                "end",
                "strand",
                "gene_id",
                "gene_name",
                "gene_biotype"
            ),
            drop = FALSE
        ]
    }

    genes <- genes[
        !duplicated(
            genes[
                ,
                c(
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    "gene_id"
                ),
                drop = FALSE
            ]
        ),
        ,
        drop = FALSE
    ]

    message("Loaded ", nrow(genes), " gene ranges.")

    genes
}

############################################################
## LeafCutter identifier parsing                            ##
############################################################

parse_leafcutter_intron_ids <- function(intron_ids) {

    intron_ids <- as.character(intron_ids)

    pattern <- "^(.+):([0-9]+):([0-9]+):(.+)$"

    matches <- regexec(
        pattern,
        intron_ids,
        perl = TRUE
    )

    parsed <- regmatches(intron_ids, matches)

    parsed_df <- data.frame(
        intron_id = intron_ids,
        seqname = NA_character_,
        intron_start = NA_integer_,
        intron_end = NA_integer_,
        cluster_suffix = NA_character_,
        cluster = NA_character_,
        strand = NA_character_,
        stringsAsFactors = FALSE
    )

    valid <- lengths(parsed) == 5

    if (any(valid)) {

        parsed_df[["seqname"]][valid] <- vapply(
            parsed[valid],
            `[[`,
            character(1),
            2
        )

        parsed_df[["intron_start"]][valid] <- as.integer(
            vapply(
                parsed[valid],
                `[[`,
                character(1),
                3
            )
        )

        parsed_df[["intron_end"]][valid] <- as.integer(
            vapply(
                parsed[valid],
                `[[`,
                character(1),
                4
            )
        )

        parsed_df[["cluster_suffix"]][valid] <- vapply(
            parsed[valid],
            `[[`,
            character(1),
            5
        )

        # Cluster-level LeafCutter results use identifiers such as:
        # chr22:clu_7_-
        parsed_df[["cluster"]][valid] <- paste(
            parsed_df[["seqname"]][valid],
            parsed_df[["cluster_suffix"]][valid],
            sep = ":"
        )

        parsed_df[["strand"]][valid] <- sub(
            "^.*_([+-])$",
            "\\1",
            parsed_df[["cluster_suffix"]][valid]
        )

        invalid_strand <-
            !parsed_df[["strand"]][valid] %in% c("+", "-")

        valid_indices <- which(valid)

        parsed_df[["strand"]][
            valid_indices[invalid_strand]
        ] <- NA_character_
    }

    if (any(!valid)) {
        warning(
            sum(!valid),
            " LeafCutter intron identifiers could not be parsed."
        )
    }

    parsed_df
}


find_intron_id_column <- function(data, table_name) {

    preferred_names <- c(
        "intron",
        "intron_id",
        "intronID",
        "event",
        "feature"
    )

    preferred_match <- preferred_names[
        preferred_names %in% colnames(data)
    ]

    if (length(preferred_match) > 0) {
        return(preferred_match[[1]])
    }

    candidate_scores <- vapply(
        data,
        function(column) {

            if (!is.character(column) && !is.factor(column)) {
                return(0)
            }

            values <- as.character(column)

            if (length(values) == 0) {
                return(0)
            }

            mean(
                grepl(
                    "^.+:[0-9]+:[0-9]+:.+$",
                    values
                ),
                na.rm = TRUE
            )
        },
        numeric(1)
    )

    if (
        length(candidate_scores) == 0 ||
        max(candidate_scores, na.rm = TRUE) == 0
    ) {
        stop(
            "Could not identify an intron-ID column in ",
            table_name,
            ". Available columns: ",
            paste(colnames(data), collapse = ", "),
            call. = FALSE
        )
    }

    names(which.max(candidate_scores))
}


harmonise_seqnames <- function(query_seqnames, reference_seqnames) {

    query_seqnames <- as.character(query_seqnames)
    reference_seqnames <- as.character(reference_seqnames)

    query_nonmissing <- query_seqnames[
        !is.na(query_seqnames)
    ]

    reference_nonmissing <- reference_seqnames[
        !is.na(reference_seqnames)
    ]

    if (
        length(query_nonmissing) == 0 ||
        length(reference_nonmissing) == 0
    ) {
        return(query_seqnames)
    }

    query_has_chr <- mean(
        startsWith(query_nonmissing, "chr")
    ) > 0.5

    reference_has_chr <- mean(
        startsWith(reference_nonmissing, "chr")
    ) > 0.5

    if (query_has_chr && !reference_has_chr) {
        query_seqnames <- sub(
            "^chr",
            "",
            query_seqnames
        )
    }

    if (!query_has_chr && reference_has_chr) {
        query_seqnames <- ifelse(
            is.na(query_seqnames),
            NA_character_,
            paste0("chr", query_seqnames)
        )
    }

    query_seqnames
}

############################################################
## Genomic annotation                                      ##
############################################################

annotate_introns_with_genes <- function(intron_table, genes) {

    annotation <- intron_table

    annotation[["gene_id"]] <- NA_character_
    annotation[["gene_name"]] <- NA_character_
    annotation[["gene_biotype"]] <- NA_character_
    annotation[["n_overlapping_genes"]] <- 0L

    annotation[["seqname_gtf"]] <- harmonise_seqnames(
        annotation[["seqname"]],
        genes[["seqname"]]
    )

    valid_introns <-
        !is.na(annotation[["seqname_gtf"]]) &
        !is.na(annotation[["intron_start"]]) &
        !is.na(annotation[["intron_end"]])

    if (!any(valid_introns)) {
        return(annotation)
    }

    intron_strand <- annotation[["strand"]][valid_introns]
    intron_strand[is.na(intron_strand)] <- "*"

    intron_gr <- GenomicRanges::GRanges(
        seqnames = annotation[["seqname_gtf"]][valid_introns],
        ranges = IRanges::IRanges(
            start = annotation[["intron_start"]][valid_introns],
            end = annotation[["intron_end"]][valid_introns]
        ),
        strand = intron_strand
    )

    gene_strand <- genes[["strand"]]
    gene_strand[
        is.na(gene_strand) |
            !gene_strand %in% c("+", "-", "*")
    ] <- "*"

    gene_gr <- GenomicRanges::GRanges(
        seqnames = genes[["seqname"]],
        ranges = IRanges::IRanges(
            start = genes[["start"]],
            end = genes[["end"]]
        ),
        strand = gene_strand
    )

    overlaps <- suppressWarnings(
        GenomicRanges::findOverlaps(
            intron_gr,
            gene_gr,
            ignore.strand = FALSE
        )
    )

    if (length(overlaps) == 0) {
        warning("No LeafCutter introns overlapped GTF gene ranges.")
        return(annotation)
    }

    valid_indices <- which(valid_introns)

    overlap_df <- data.frame(
        intron_row = valid_indices[
            S4Vectors::queryHits(overlaps)
        ],
        gene_row = S4Vectors::subjectHits(overlaps),
        stringsAsFactors = FALSE
    )

    overlap_df[["gene_id"]] <-
        genes[["gene_id"]][overlap_df[["gene_row"]]]

    overlap_df[["gene_name"]] <-
        genes[["gene_name"]][overlap_df[["gene_row"]]]

    overlap_df[["gene_biotype"]] <-
        genes[["gene_biotype"]][overlap_df[["gene_row"]]]

    collapsed <- overlap_df |>
        dplyr::group_by(intron_row) |>
        dplyr::summarise(
            gene_id = collapse_unique(gene_id),
            gene_name = collapse_unique(gene_name),
            gene_biotype = collapse_unique(gene_biotype),
            n_overlapping_genes = dplyr::n_distinct(
                gene_id[
                    !is.na(gene_id) &
                        gene_id != ""
                ]
            ),
            .groups = "drop"
        ) |>
        as.data.frame()

    target_rows <- collapsed[["intron_row"]]

    annotation[["gene_id"]][target_rows] <-
        collapsed[["gene_id"]]

    annotation[["gene_name"]][target_rows] <-
        collapsed[["gene_name"]]

    annotation[["gene_biotype"]][target_rows] <-
        collapsed[["gene_biotype"]]

    annotation[["n_overlapping_genes"]][target_rows] <-
        collapsed[["n_overlapping_genes"]]

    annotation
}


append_intron_annotation <- function(
    input_table,
    intron_column,
    genes
) {

    intron_ids <- as.character(
        input_table[[intron_column]]
    )

    parsed <- parse_leafcutter_intron_ids(intron_ids)

    annotation <- annotate_introns_with_genes(
        parsed,
        genes
    )

    annotation_columns <- setdiff(
        colnames(annotation),
        "intron_id"
    )

    output <- cbind(
        input_table,
        annotation[
            ,
            annotation_columns,
            drop = FALSE
        ]
    )

    if (!"intron_id" %in% colnames(output)) {
        output[["intron_id"]] <- intron_ids
    }

    preferred_columns <- c(
        "intron_id",
        "seqname",
        "intron_start",
        "intron_end",
        "strand",
        "cluster",
        "cluster_suffix",
        "gene_id",
        "gene_name",
        "gene_biotype",
        "n_overlapping_genes"
    )

    preferred_columns <- preferred_columns[
        preferred_columns %in% colnames(output)
    ]

    output[
        ,
        c(
            preferred_columns,
            setdiff(
                colnames(output),
                preferred_columns
            )
        ),
        drop = FALSE
    ]
}

############################################################
## Parse options and validate inputs                        ##
############################################################

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))

opt <- list(
    cluster_file = cli[["cluster_file"]],
    effect_file = cli[["effect_file"]],
    counts_file = cli[["counts_file"]],
    groups_file = cli[["groups_file"]],
    gtf_file = cli[["gtf_file"]],
    output_prefix = cli[["output_prefix"]],
    organism = cli[["organism"]] %||% "human",
    alpha = as.numeric(cli[["alpha"]] %||% "0.05"),
    min_abs_effect = as.numeric(
        cli[["min_abs_effect"]] %||% "0.10"
    )
)

required_files <- c(
    cluster_file = "cluster_file",
    effect_file = "effect_file",
    counts_file = "counts_file",
    groups_file = "groups_file",
    gtf_file = "gtf_file"
)

for (option_name in names(required_files)) {
    validate_file(
        opt[[option_name]],
        required_files[[option_name]]
    )
}

if (!is_valid_string(opt[["output_prefix"]])) {
    stop(
        "Missing required argument: output_prefix",
        call. = FALSE
    )
}

if (
    is.na(opt[["alpha"]]) ||
    opt[["alpha"]] <= 0 ||
    opt[["alpha"]] > 1
) {
    stop(
        "--alpha must be greater than 0 and no greater than 1.",
        call. = FALSE
    )
}

if (
    is.na(opt[["min_abs_effect"]]) ||
    opt[["min_abs_effect"]] < 0
) {
    stop(
        "--min_abs_effect must be zero or greater.",
        call. = FALSE
    )
}

############################################################
## Pathway enrichment helper functions                     ##
############################################################

split_gene_values <- function(values) {
  
  values <- as.character(values)
  values <- values[
    !is.na(values) &
      values != ""
  ]
  
  if (length(values) == 0) {
    return(character(0))
  }
  
  genes <- unlist(
    strsplit(
      values,
      split = ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )
  
  genes <- trimws(genes)
  
  unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )
}


strip_ensembl_version <- function(gene_ids) {
  
  gene_ids <- as.character(gene_ids)
  
  sub(
    "\\.[0-9]+$",
    "",
    gene_ids
  )
}


write_empty_enrichment_table <- function(filename) {
  
  empty_table <- data.frame(
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalue = numeric(0),
    geneID = character(0),
    Count = integer(0),
    stringsAsFactors = FALSE
  )
  
  write_tsv(
    empty_table,
    filename
  )
}


save_enrichment_plots <- function(
    enrichment_result,
    output_prefix,
    title,
    show_category = 20
) {
  
  enrichment_df <- as.data.frame(enrichment_result)
  
  if (nrow(enrichment_df) == 0) {
    return(invisible(NULL))
  }
  
  dot_plot <- enrichplot::dotplot(
    enrichment_result,
    showCategory = show_category
  ) +
    ggplot2::ggtitle(title)
  
  ggplot2::ggsave(
    filename = paste0(
      output_prefix,
      ".dotplot.pdf"
    ),
    plot = dot_plot,
    width = 10,
    height = 7
  )
  
  ggplot2::ggsave(
    filename = paste0(
      output_prefix,
      ".dotplot.png"
    ),
    plot = dot_plot,
    width = 10,
    height = 7,
    dpi = 200
  )
  
  bar_plot <- enrichplot::barplot(
    enrichment_result,
    showCategory = show_category
  ) +
    ggplot2::ggtitle(title)
  
  ggplot2::ggsave(
    filename = paste0(
      output_prefix,
      ".barplot.pdf"
    ),
    plot = bar_plot,
    width = 9,
    height = 7
  )
  
  ggplot2::ggsave(
    filename = paste0(
      output_prefix,
      ".barplot.png"
    ),
    plot = bar_plot,
    width = 9,
    height = 7,
    dpi = 200
  )
  
  invisible(NULL)
}


write_enrichment_result <- function(
    enrichment_result,
    filename
) {
  
  enrichment_df <- as.data.frame(enrichment_result)
  
  if (nrow(enrichment_df) > 0) {
    
    enrichment_df <- enrichment_df[
      order(
        enrichment_df[["p.adjust"]],
        enrichment_df[["pvalue"]],
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    write_tsv(
      enrichment_df,
      filename
    )
    
  } else {
    
    write_empty_enrichment_table(filename)
  }
  
  invisible(enrichment_df)
}

############################################################
## Load libraries                                           ##
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
})

############################################################
## Read inputs                                              ##
############################################################

message("Reading LeafCutter cluster results.")
cluster_results <- read_leafcutter_table(
    opt[["cluster_file"]]
)

message("Reading LeafCutter intron effect sizes.")
effect_results <- read_leafcutter_table(
    opt[["effect_file"]]
)

message("Reading LeafCutter numerator counts.")
counts_results <- read_leafcutter_table(
    opt[["counts_file"]]
)

message("Reading LeafCutter sample groups.")
groups <- data.table::fread(
    opt[["groups_file"]],
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE
)

if (ncol(groups) < 2) {
    stop(
        "The LeafCutter groups file must contain at least two columns.",
        call. = FALSE
    )
}

group_column_names <- c(
    "sample",
    "group",
    if (ncol(groups) > 2) {
        paste0("covariate_", seq_len(ncol(groups) - 2))
    } else {
        character(0)
    }
)

colnames(groups) <- group_column_names

require_columns(
    cluster_results,
    c(
        "status",
        "loglr",
        "df",
        "p",
        "cluster",
        "p.adjust"
    ),
    "LeafCutter cluster-significance table"
)

############################################################
## Read GTF and annotate introns                             ##
############################################################

genes <- read_gtf_genes(opt[["gtf_file"]])

effect_intron_column <- find_intron_id_column(
    effect_results,
    "LeafCutter effect-size table"
)

message(
    "Using effect-size column '",
    effect_intron_column,
    "' as the intron identifier."
)

effect_results_annotated <- append_intron_annotation(
    effect_results,
    effect_intron_column,
    genes
)

counts_intron_column <- find_intron_id_column(
    counts_results,
    "LeafCutter numerator-count table"
)

message(
    "Using count-table column '",
    counts_intron_column,
    "' as the intron identifier."
)

counts_results_annotated <- append_intron_annotation(
    counts_results,
    counts_intron_column,
    genes
)

############################################################
## Cluster-level annotation and significance                ##
############################################################

cluster_annotation <- effect_results_annotated |>
    dplyr::filter(
        !is.na(cluster),
        cluster != ""
    ) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
        seqname = collapse_unique(seqname),
        cluster_start = safe_min(intron_start),
        cluster_end = safe_max(intron_end),
        strand = collapse_unique(strand),
        gene_id = collapse_semicolon_values(gene_id),
        gene_name = collapse_semicolon_values(gene_name),
        gene_biotype = collapse_semicolon_values(
            gene_biotype
        ),
        n_introns = dplyr::n_distinct(intron_id),
        n_annotated_introns = sum(
            !is.na(gene_id) &
                gene_id != ""
        ),
        max_abs_deltapsi = if (
            "deltapsi" %in% colnames(
                dplyr::cur_data_all()
            ) &&
            !all(is.na(deltapsi))
        ) {
            max(abs(as.numeric(deltapsi)), na.rm = TRUE)
        } else {
            NA_real_
        },
        .groups = "drop"
    ) |>
    as.data.frame()

cluster_results_annotated <- merge(
    cluster_results,
    cluster_annotation,
    by = "cluster",
    all.x = TRUE,
    sort = FALSE
)

cluster_results_annotated[["tested"]] <-
    cluster_results_annotated[["status"]] == "Success" &
    !is.na(cluster_results_annotated[["p.adjust"]])

cluster_results_annotated[["significant"]] <-
    cluster_results_annotated[["tested"]] &
    cluster_results_annotated[["p.adjust"]] <
        opt[["alpha"]]

cluster_results_annotated <- cluster_results_annotated[
    order(
        cluster_results_annotated[["p.adjust"]],
        na.last = TRUE
    ),
    ,
    drop = FALSE
]

############################################################
## Add cluster statistics to intron-level tables            ##
############################################################

cluster_statistics <- cluster_results_annotated[
    ,
    c(
        "cluster",
        "status",
        "loglr",
        "df",
        "p",
        "p.adjust",
        "tested",
        "significant"
    ),
    drop = FALSE
]

colnames(cluster_statistics)[
    colnames(cluster_statistics) != "cluster"
] <- paste0(
    "cluster_",
    colnames(cluster_statistics)[
        colnames(cluster_statistics) != "cluster"
    ]
)

effect_results_annotated <- merge(
    effect_results_annotated,
    cluster_statistics,
    by = "cluster",
    all.x = TRUE,
    sort = FALSE
)

counts_results_annotated <- merge(
    counts_results_annotated,
    cluster_statistics,
    by = "cluster",
    all.x = TRUE,
    sort = FALSE
)

if ("deltapsi" %in% colnames(effect_results_annotated)) {

    effect_results_annotated[["deltapsi"]] <-
        as.numeric(effect_results_annotated[["deltapsi"]])

    effect_results_annotated[["passes_effect_threshold"]] <-
        !is.na(effect_results_annotated[["deltapsi"]]) &
        abs(effect_results_annotated[["deltapsi"]]) >=
            opt[["min_abs_effect"]]

} else {

    warning(
        "No 'deltapsi' column was found in the effect-size table. ",
        "The effect-size threshold will not be applied."
    )

    effect_results_annotated[["passes_effect_threshold"]] <-
        NA
}

effect_results_annotated[["significant_intron"]] <-
    effect_results_annotated[["cluster_significant"]] %in% TRUE &
    effect_results_annotated[["passes_effect_threshold"]] %in% TRUE


############################################################
## Standardised gene-level summary                         ##
############################################################

message("Generating standardised LeafCutter gene-level summary.")


# Generate a reliable gene-ID-to-name lookup directly from the GTF.
gene_name_lookup <- genes[
    !is.na(genes[["gene_id"]]) &
        genes[["gene_id"]] != "",
    c(
        "gene_id",
        "gene_name"
    ),
    drop = FALSE
]

gene_name_lookup <- gene_name_lookup[
    !duplicated(gene_name_lookup[["gene_id"]]),
    ,
    drop = FALSE
]


# Expand semicolon-separated gene annotations so that each
# gene-cluster association occupies one row.
expanded_gene_cluster_rows <- lapply(
    seq_len(nrow(effect_results_annotated)),
    function(row_index) {

        gene_ids <- split_gene_values(
            effect_results_annotated[["gene_id"]][row_index]
        )

        if (length(gene_ids) == 0) {
            return(NULL)
        }

        data.frame(
            gene_id = gene_ids,
            cluster = rep(
                effect_results_annotated[["cluster"]][row_index],
                length(gene_ids)
            ),
            deltapsi = rep(
                suppressWarnings(
                    as.numeric(
                        effect_results_annotated[["deltapsi"]][row_index]
                    )
                ),
                length(gene_ids)
            ),
            cluster_tested = rep(
                effect_results_annotated[["cluster_tested"]][row_index]
                    %in% TRUE,
                length(gene_ids)
            ),
            significant_intron = rep(
                effect_results_annotated[["significant_intron"]][row_index]
                    %in% TRUE,
                length(gene_ids)
            ),
            stringsAsFactors = FALSE
        )
    }
)


expanded_gene_cluster_rows <- Filter(
    Negate(is.null),
    expanded_gene_cluster_rows
)


if (length(expanded_gene_cluster_rows) > 0) {

    gene_cluster_introns <- do.call(
        rbind,
        expanded_gene_cluster_rows
    )

    gene_cluster_introns <- merge(
        gene_cluster_introns,
        gene_name_lookup,
        by = "gene_id",
        all.x = TRUE,
        sort = FALSE
    )


    # Reduce each gene-cluster association to one event.
    #
    # The representative signed effect is the intron ΔPSI with
    # the largest absolute magnitude in that cluster.
    gene_cluster_events <- gene_cluster_introns |>
        dplyr::group_by(
            gene_id,
            gene_name,
            cluster
        ) |>
        dplyr::summarise(
            detected = any(
                cluster_tested %in% TRUE
            ),
            significant = any(
                significant_intron %in% TRUE
            ),
            representative_effect = if (
                length(deltapsi) == 0 ||
                all(is.na(deltapsi))
            ) {
                NA_real_
            } else {
                deltapsi[
                    which.max(abs(deltapsi))
                ]
            },
            .groups = "drop"
        ) |>
        as.data.frame()


    gene_summary <- gene_cluster_events |>
        dplyr::group_by(
            gene_id,
            gene_name
        ) |>
        dplyr::summarise(
            method = "LeafCutter",
            detected = any(
                detected %in% TRUE
            ),
            significant = any(
                significant %in% TRUE
            ),
            n_events = dplyr::n_distinct(
                cluster[
                    detected %in% TRUE
                ]
            ),
            n_significant_events = dplyr::n_distinct(
                cluster[
                    significant %in% TRUE
                ]
            ),
            max_abs_effect = if (
                all(is.na(representative_effect))
            ) {
                NA_real_
            } else {
                max(
                    abs(representative_effect),
                    na.rm = TRUE
                )
            },
            effect_direction = {

                direction_effects <- representative_effect[
                    significant %in% TRUE &
                        !is.na(representative_effect)
                ]

                # For non-significant genes, retain the direction
                # of detected events rather than returning none.
                if (length(direction_effects) == 0) {
                    direction_effects <- representative_effect[
                        detected %in% TRUE &
                            !is.na(representative_effect)
                    ]
                }

                has_positive <- any(
                    direction_effects > 0
                )

                has_negative <- any(
                    direction_effects < 0
                )

                if (has_positive && has_negative) {
                    "mixed"
                } else if (has_positive) {
                    "increased"
                } else if (has_negative) {
                    "decreased"
                } else {
                    "none"
                }
            },
            event_types = "intron_cluster",
            .groups = "drop"
        ) |>
        as.data.frame()


    gene_summary <- gene_summary[
        order(
            -as.integer(gene_summary[["significant"]]),
            -gene_summary[["n_significant_events"]],
            -gene_summary[["max_abs_effect"]],
            gene_summary[["gene_name"]],
            gene_summary[["gene_id"]],
            na.last = TRUE
        ),
        ,
        drop = FALSE
    ]


    gene_summary <- gene_summary[
        ,
        c(
            "gene_id",
            "gene_name",
            "method",
            "detected",
            "significant",
            "n_events",
            "n_significant_events",
            "max_abs_effect",
            "effect_direction",
            "event_types"
        ),
        drop = FALSE
    ]

} else {

    warning(
        "No annotated LeafCutter gene-cluster associations ",
        "were available for the gene-level summary."
    )

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
## Summary table                                            ##
############################################################

summary_table <- data.frame(
    metric = c(
        "organism",
        "alpha",
        "minimum_absolute_delta_psi",
        "samples",
        "groups",
        "clusters_total",
        "clusters_tested",
        "clusters_significant",
        "introns_with_effect_sizes",
        "significant_introns",
        "annotated_introns",
        "unannotated_introns"
    ),
    value = c(
        opt[["organism"]],
        opt[["alpha"]],
        opt[["min_abs_effect"]],
        nrow(groups),
        length(unique(groups[["group"]])),
        nrow(cluster_results_annotated),
        sum(
            cluster_results_annotated[["tested"]],
            na.rm = TRUE
        ),
        sum(
            cluster_results_annotated[["significant"]],
            na.rm = TRUE
        ),
        nrow(effect_results_annotated),
        sum(
            effect_results_annotated[["significant_intron"]],
            na.rm = TRUE
        ),
        sum(
            !is.na(effect_results_annotated[["gene_id"]]) &
                effect_results_annotated[["gene_id"]] != ""
        ),
        sum(
            is.na(effect_results_annotated[["gene_id"]]) |
                effect_results_annotated[["gene_id"]] == ""
        )
    ),
    stringsAsFactors = FALSE
)

############################################################
## Pathway over-representation analysis                    ##
############################################################

message("Preparing genes for pathway enrichment analysis.")

if (opt[["organism"]] == "human") {
  
  suppressPackageStartupMessages(
    library(org.Hs.eg.db)
  )
  
  annotation_database <- org.Hs.eg.db
  msig_species <- "Homo sapiens"
  kegg_organism <- "hsa"
  
} else if (opt[["organism"]] == "mouse") {
  
  suppressPackageStartupMessages(
    library(org.Mm.eg.db)
  )
  
  annotation_database <- org.Mm.eg.db
  msig_species <- "Mus musculus"
  kegg_organism <- "mmu"
  
} else {
  
  stop(
    "Unsupported organism for enrichment analysis: ",
    opt[["organism"]],
    ". Supported values are 'human' and 'mouse'.",
    call. = FALSE
  )
}


# All genes associated with successfully tested LeafCutter clusters.
tested_cluster_gene_ids <- split_gene_values(
  cluster_results_annotated[["gene_id"]][
    cluster_results_annotated[["tested"]] %in% TRUE
  ]
)

# Genes associated with statistically significant clusters.
significant_cluster_gene_ids <- split_gene_values(
  cluster_results_annotated[["gene_id"]][
    cluster_results_annotated[["significant"]] %in% TRUE
  ]
)

# More stringent set requiring a cluster-level delta-PSI effect.
effect_filtered_gene_ids <- split_gene_values(
  cluster_results_annotated[["gene_id"]][
    cluster_results_annotated[["significant"]] %in% TRUE &
      !is.na(
        cluster_results_annotated[["max_abs_deltapsi"]]
      ) &
      cluster_results_annotated[["max_abs_deltapsi"]] >=
      opt[["min_abs_effect"]]
  ]
)

tested_cluster_gene_ids <- strip_ensembl_version(
  tested_cluster_gene_ids
)

significant_cluster_gene_ids <- strip_ensembl_version(
  significant_cluster_gene_ids
)

effect_filtered_gene_ids <- strip_ensembl_version(
  effect_filtered_gene_ids
)


message(
  "Enrichment universe contains ",
  length(tested_cluster_gene_ids),
  " tested genes."
)

message(
  "Significant foreground contains ",
  length(significant_cluster_gene_ids),
  " genes."
)

message(
  "Effect-filtered foreground contains ",
  length(effect_filtered_gene_ids),
  " genes."
)


# Map tested Ensembl IDs to symbols and Entrez IDs.
gene_annotation <- AnnotationDbi::select(
  annotation_database,
  keys = tested_cluster_gene_ids,
  keytype = "ENSEMBL",
  columns = c(
    "SYMBOL",
    "ENTREZID"
  )
)

gene_annotation <- gene_annotation[
  !is.na(gene_annotation[["ENSEMBL"]]),
  ,
  drop = FALSE
]

gene_annotation <- gene_annotation[
  !duplicated(
    gene_annotation[
      ,
      c(
        "ENSEMBL",
        "SYMBOL",
        "ENTREZID"
      ),
      drop = FALSE
    ]
  ),
  ,
  drop = FALSE
]


tested_symbols <- unique(
  gene_annotation[["SYMBOL"]][
    !is.na(gene_annotation[["SYMBOL"]]) &
      gene_annotation[["SYMBOL"]] != ""
  ]
)

tested_entrez <- unique(
  gene_annotation[["ENTREZID"]][
    !is.na(gene_annotation[["ENTREZID"]]) &
      gene_annotation[["ENTREZID"]] != ""
  ]
)


significant_symbols <- unique(
  gene_annotation[["SYMBOL"]][
    gene_annotation[["ENSEMBL"]] %in%
      significant_cluster_gene_ids &
      !is.na(gene_annotation[["SYMBOL"]]) &
      gene_annotation[["SYMBOL"]] != ""
  ]
)

significant_entrez <- unique(
  gene_annotation[["ENTREZID"]][
    gene_annotation[["ENSEMBL"]] %in%
      significant_cluster_gene_ids &
      !is.na(gene_annotation[["ENTREZID"]]) &
      gene_annotation[["ENTREZID"]] != ""
  ]
)


effect_filtered_symbols <- unique(
  gene_annotation[["SYMBOL"]][
    gene_annotation[["ENSEMBL"]] %in%
      effect_filtered_gene_ids &
      !is.na(gene_annotation[["SYMBOL"]]) &
      gene_annotation[["SYMBOL"]] != ""
  ]
)

effect_filtered_entrez <- unique(
  gene_annotation[["ENTREZID"]][
    gene_annotation[["ENSEMBL"]] %in%
      effect_filtered_gene_ids &
      !is.na(gene_annotation[["ENTREZID"]]) &
      gene_annotation[["ENTREZID"]] != ""
  ]
)


# Write the gene sets used for enrichment.
enrichment_gene_sets <- data.frame(
  gene_set = c(
    rep(
      "tested_universe",
      length(tested_cluster_gene_ids)
    ),
    rep(
      "significant",
      length(significant_cluster_gene_ids)
    ),
    rep(
      "significant_effect_filtered",
      length(effect_filtered_gene_ids)
    )
  ),
  gene_id = c(
    tested_cluster_gene_ids,
    significant_cluster_gene_ids,
    effect_filtered_gene_ids
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  enrichment_gene_sets,
  paste0(
    opt[["output_prefix"]],
    ".leafcutter.enrichment_gene_sets.tsv"
  )
) 

############################################################
## Hallmark ORA                                            ##
############################################################

message("Running Hallmark over-representation analysis.")

hallmark_sets <- msigdbr::msigdbr(
  species = msig_species,
  category = "H"
)

hallmark_term2gene <- unique(
  hallmark_sets[
    ,
    c(
      "gs_name",
      "gene_symbol"
    ),
    drop = FALSE
  ]
)

colnames(hallmark_term2gene) <- c(
  "term",
  "gene"
)


run_hallmark_ora <- function(
    foreground,
    universe,
    result_label
) {
  
  output_base <- paste0(
    opt[["output_prefix"]],
    ".leafcutter.",
    result_label,
    ".hallmark_ora"
  )
  
  if (length(foreground) < 3) {
    
    warning(
      "Fewer than three genes available for ",
      result_label,
      " Hallmark ORA."
    )
    
    write_empty_enrichment_table(
      paste0(output_base, ".tsv")
    )
    
    return(NULL)
  }
  
  result <- clusterProfiler::enricher(
    gene = foreground,
    universe = universe,
    TERM2GENE = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500
  )
  
  write_enrichment_result(
    result,
    paste0(output_base, ".tsv")
  )
  
  save_enrichment_plots(
    result,
    output_base,
    paste0(
      "Hallmark ORA: ",
      result_label
    )
  )
  
  result
}


hallmark_significant <- run_hallmark_ora(
  foreground = significant_symbols,
  universe = tested_symbols,
  result_label = "significant"
)

hallmark_effect_filtered <- run_hallmark_ora(
  foreground = effect_filtered_symbols,
  universe = tested_symbols,
  result_label = "significant_effect_filtered"
)

############################################################
## Gene Ontology ORA                                       ##
############################################################

message("Running Gene Ontology over-representation analysis.")

run_go_ora <- function(
    foreground,
    universe,
    ontology,
    result_label
) {
  
  ontology_lower <- tolower(ontology)
  
  output_base <- paste0(
    opt[["output_prefix"]],
    ".leafcutter.",
    result_label,
    ".go_",
    ontology_lower,
    "_ora"
  )
  
  if (length(foreground) < 3) {
    
    warning(
      "Fewer than three genes available for ",
      result_label,
      " GO ",
      ontology,
      " ORA."
    )
    
    write_empty_enrichment_table(
      paste0(output_base, ".tsv")
    )
    
    return(NULL)
  }
  
  result <- clusterProfiler::enrichGO(
    gene = foreground,
    universe = universe,
    OrgDb = annotation_database,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 1000,
    readable = TRUE
  )
  
  write_enrichment_result(
    result,
    paste0(output_base, ".tsv")
  )
  
  save_enrichment_plots(
    result,
    output_base,
    paste0(
      "GO ",
      ontology,
      " ORA: ",
      result_label
    )
  )
  
  result
}


for (ontology in c("BP", "MF", "CC")) {
  
  run_go_ora(
    foreground = significant_entrez,
    universe = tested_entrez,
    ontology = ontology,
    result_label = "significant"
  )
  
  run_go_ora(
    foreground = effect_filtered_entrez,
    universe = tested_entrez,
    ontology = ontology,
    result_label = "significant_effect_filtered"
  )
}

############################################################
## KEGG ORA                                                ##
############################################################

message("Running KEGG over-representation analysis.")

run_kegg_ora <- function(
    foreground,
    universe,
    result_label
) {
  
  output_base <- paste0(
    opt[["output_prefix"]],
    ".leafcutter.",
    result_label,
    ".kegg_ora"
  )
  
  if (length(foreground) < 3) {
    
    warning(
      "Fewer than three genes available for ",
      result_label,
      " KEGG ORA."
    )
    
    write_empty_enrichment_table(
      paste0(output_base, ".tsv")
    )
    
    return(NULL)
  }
  
  result <- tryCatch(
    {
      clusterProfiler::enrichKEGG(
        gene = foreground,
        universe = universe,
        organism = kegg_organism,
        keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 5,
        maxGSSize = 500,
        use_internal_data = FALSE
      )
    },
    error = function(error_condition) {
      
      warning(
        "KEGG enrichment failed: ",
        conditionMessage(error_condition)
      )
      
      NULL
    }
  )
  
  if (is.null(result)) {
    
    write_empty_enrichment_table(
      paste0(output_base, ".tsv")
    )
    
    return(NULL)
  }
  
  write_enrichment_result(
    result,
    paste0(output_base, ".tsv")
  )
  
  save_enrichment_plots(
    result,
    output_base,
    paste0(
      "KEGG ORA: ",
      result_label
    )
  )
  
  result
}


kegg_significant <- run_kegg_ora(
  foreground = significant_entrez,
  universe = tested_entrez,
  result_label = "significant"
)

kegg_effect_filtered <- run_kegg_ora(
  foreground = effect_filtered_entrez,
  universe = tested_entrez,
  result_label = "significant_effect_filtered"
)


############################################################
## Write outputs                                            ##
############################################################

prefix <- opt[["output_prefix"]]

write_tsv(
    gene_summary,
    paste0(
        prefix,
        ".leafcutter.gene_summary.tsv"
    )
)

write_tsv(
    cluster_results_annotated,
    paste0(
        prefix,
        ".leafcutter.clusters.annotated.tsv"
    )
)

write_tsv(
    cluster_results_annotated[
        cluster_results_annotated[["significant"]] %in% TRUE,
        ,
        drop = FALSE
    ],
    paste0(
        prefix,
        ".leafcutter.clusters.significant.tsv"
    )
)

write_tsv(
    effect_results_annotated,
    paste0(
        prefix,
        ".leafcutter.introns.annotated.tsv"
    )
)

write_tsv(
    effect_results_annotated[
        effect_results_annotated[["significant_intron"]] %in% TRUE,
        ,
        drop = FALSE
    ],
    paste0(
        prefix,
        ".leafcutter.introns.significant.tsv"
    )
)

write_tsv(
    counts_results_annotated,
    paste0(
        prefix,
        ".leafcutter.counts.annotated.tsv"
    )
)

write_tsv(
    groups,
    paste0(
        prefix,
        ".leafcutter.sample_groups.tsv"
    )
)

write_tsv(
    summary_table,
    paste0(
        prefix,
        ".leafcutter.summary.tsv"
    )
)

sink(
    paste0(
        prefix,
        ".leafcutter.R_sessionInfo.log"
    )
)

print(sessionInfo())
sink()

message("LeafCutter downstream annotation completed.")