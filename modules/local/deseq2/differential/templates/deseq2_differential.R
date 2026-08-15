#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
  !is.null(input) && nzchar(trimws(input))
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

  ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

  if (ext == "tsv" || ext == "txt") {
    separator <- "\\t"
  } else if (ext == "csv") {
    separator <- ","
  } else {
    stop(paste("Unknown separator for", ext))
  }

  read.delim(
    file,
    sep = separator,
    header = header,
    row.names = row.names,
    check.names = check.names
  )
}

#
#' Turn “null” or empty strings into actual NULL
#'
#' @param x Input option
#'
#' @return NULL or x
#'
nullify <- function(x) {
  if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Run DESeq2 differential expression results with shared settings
#'
#' A thin wrapper around DESeq2::results() that injects the common parameters
#' (lfcThreshold, altHypothesis, independentFiltering, alpha, pAdjustMethod, minmu)
#' from the global `opt` list and the `dds` object. Use the `...` to supply the
#' branch-specific argument (`name = …` or `contrast = …`).
#'
#' @param ... Additional arguments passed to `results()`, e.g. `name = opt\$contrast_string`
#'   or `contrast = c(variable, target, reference)`.
#' @return A `DESeqResults` object containing the differential expression test results.
#'
run_results <- function(...) {
  results(
    dds,
    lfcThreshold         = opt\$lfc_threshold,
    altHypothesis        = opt\$alt_hypothesis,
    independentFiltering = opt\$independent_filtering,
    alpha                = opt\$alpha,
    pAdjustMethod        = opt\$p_adjust_method,
    minmu                = opt\$minmu,
    ...
  )
}
#' Shrink log2 fold changes via apeglm with shared settings
#'
#' A thin wrapper around DESeq2::lfcShrink() that injects the common parameters
#' (type = "apeglm") and the `dds` object. Use `...` to supply the branch-
#' specific argument, for example a coefficient name or contrast vector.
#'
#' @param ... Additional arguments passed to `lfcShrink()`, for example
#' `coef = "condition_treated_vs_control"` or
#' `contrast = c(variable, target, reference)`.
#' @return A `DESeqResults` object with shrunken log2 fold changes.
run_shrink <- function(...) {
  lfcShrink(
    dds,
    type = "apeglm",
    ...
  )
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
  output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
  count_file = '$counts',
  sample_file = '$samplesheet',
  contrast_variable = '$contrast_variable',
  reference_level = '$reference',
  target_level = '$target',
  contrast_string    = "$comparison",             # Optional full (complex) contrast expression comparison
  formula            = "$formula",              # User-specified formula
  blocking_variables = NULL,
  control_genes_file = '$control_genes_file',
  transcript_lengths_file = '$transcript_lengths_file',
  sizefactors_from_controls = FALSE,
  gene_id_col = "gene_id",
  sample_id_col = "experiment_accession",
  organism = "human",
  subset_to_contrast_samples = FALSE,
  exclude_samples_col = NULL,
  exclude_samples_values = NULL,
  test = "Wald",
  fit_type = "parametric",
  sf_type = 'ratio',
  min_replicates_for_replace = 7,
  use_t = FALSE,
  lfc_threshold = 0,
  alt_hypothesis = 'greaterAbs',
  independent_filtering = TRUE,
  p_adjust_method = 'BH',
  alpha = 0.1,
  minmu = 0.5,
  vs_method = 'vst', # 'rlog', 'vst', or 'rlog,vst'
  shrink_lfc = TRUE,
  cores = 1,
  vs_blind = TRUE,
  vst_nsub = 1000,
  round_digits = NULL,
  seed = NULL
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}


if ( ! is.null(opt\$round_digits)){
  opt\$round_digits <- as.numeric(opt\$round_digits)
}

# If there is no option supplied, convert string "null" to NULL
keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level", "seed", "blocking_variables", "transcript_lengths_file")
opt[keys] <- lapply(opt[keys], nullify)

contrast_variable <- make.names(opt[["contrast_variable"]])

if ( ! is.null(opt\$seed)){
  set.seed(opt\$seed)
}

# Check if required parameters have been provided
required_opts <- c('contrast_variable', 'output_prefix')

if (!is_valid_string(opt[["contrast_string"]])) {
  required_opts <- c(required_opts, 'reference_level', 'target_level')
}

missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
  stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
  if (! is_valid_string(opt[[file_input]])) {
    stop(paste("Please provide", file_input), call. = FALSE)
  }

  if (! file.exists(opt[[file_input]])){
    stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
  }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DESeq2)
library(BiocParallel)
library(limma)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(AnnotationDbi)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(GO.db)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
  read_delim_flexible(
    file = opt\$count_file,
    header = TRUE,
    row.names = opt\$gene_id_col,
    check.names = FALSE
  )
sample.sheet <- read_delim_flexible(file = opt\$sample_file)

# Deal with spaces that may be in sample column
opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
  stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
  sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
  stop(paste(
    length(missing_samples),
    'specified samples missing from count table:',
    paste(missing_samples, collapse = ',')
  ))
} else{
  # Save any non-count data, will gene metadata etc we might need later
  noncount.table <-
    count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
  count.table <- count.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

blocking.vars <- c()

if (!contrast_variable %in% colnames(sample.sheet)) {
  stop(
    paste0(
      'Chosen contrast variable \"',
      contrast_variable,
      '\" not in sample sheet'
    )
  )
}

# For the standard target-vs-reference comparison, validate both levels and
# explicitly set the requested reference level in the DESeq2 design.
if (!is_valid_string(opt[["contrast_string"]])) {
  if (any(!c(opt[["reference_level"]], opt[["target_level"]]) %in% sample.sheet[[contrast_variable]])) {
    stop(
      paste(
        'Please choose reference and treatment levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
      )
    )
  }

  sample.sheet[[contrast_variable]] <- relevel(
    as.factor(sample.sheet[[contrast_variable]]),
    ref = opt[["reference_level"]]
  )
}

if (!is_valid_string(opt[["formula"]]) && is_valid_string(opt[["blocking_variables"]])) {
  blocking.vars = make.names(unlist(strsplit(opt[["blocking_variables"]], split = ';')))
  if (!all(blocking.vars %in% colnames(sample.sheet))) {
    missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
    stop(
      paste(
        'Blocking variables', missing_block,
        'do not correspond to sample sheet columns.'
      )
    )
  }
}

# Optionally, subset to only the samples involved in the contrast
if (opt[["subset_to_contrast_samples"]]){
  if (!is_valid_string(opt[["reference_level"]]) || !is_valid_string(opt[["target_level"]])) {
    stop("subset_to_contrast_samples requires reference_level and target_level.")
  }

  sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt[["target_level"]], opt[["reference_level"]])
  selected_samples <- sample.sheet[sample_selector, opt[["sample_id_col"]]]
  count.table <- count.table[, selected_samples]
  sample.sheet <- sample.sheet[selected_samples, ]
}
# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((is_valid_string(opt\$exclude_samples_col)) && (is_valid_string(opt\$exclude_samples_values))){
  exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))

  if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
    stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
  }

  print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
  sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

  selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
  count.table <- count.table[, selected_samples]
  sample.sheet <- sample.sheet[selected_samples, ]
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts

if (is_valid_string(opt\$formula)) {
  message("Using user-specified formula: ", opt\$formula)
  model  <- opt\$formula
}  else {
  model <- '~ 0'

  if (is_valid_string(opt\$blocking_variables)) {
    model <- paste(model, paste(blocking.vars, collapse = ' + '), sep=' + ')
  }

  # Make sure all the appropriate variables are factors

  for (v in c(blocking.vars, contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
  }

  # Variable of interest goes last, see
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

  model <- paste(model, contrast_variable, sep = ' + ')
}
message("Final design formula: ", model)

################################################
################################################
## Run DESeq2 processes                       ##
################################################
################################################

if (opt\$control_genes_file != ''){
  control_genes <- readLines(opt\$control_genes_file)
  if (! opt\$sizefactors_from_controls){
    count.table <- count.table[setdiff(rownames(count.table), control_genes),]
  }
}

dds <- DESeqDataSetFromMatrix(
  countData = round(count.table),
  colData = sample.sheet,
  design = as.formula(model)
)

print(resultsNames(dds))

# Build in transcript lengths. Copying what tximport does here:
# https://github.com/thelovelab/DESeq2/blob/6947d5bc629015fb8ffb2453a91b71665a164483/R/AllClasses.R#L409

if (!is.null(opt\$transcript_lengths_file)){
  lengths <-
    read_delim_flexible(
      file = opt\$transcript_lengths_file,
      header = TRUE,
      row.names = opt\$gene_id_col,
      check.names = FALSE
    )
  lengths <- lengths[rownames(count.table), colnames(count.table)]
  dimnames(lengths) <- dimnames(dds)
  assays(dds)[["avgTxLength"]] <- lengths
}

if (opt\$control_genes_file != '' && opt\$sizefactors_from_controls){
  print(paste('Estimating size factors using', length(control_genes), 'control genes'))
  dds <- estimateSizeFactors(dds, controlGenes=rownames(count.table) %in% control_genes)
}

dds <- DESeq(
  dds,
  test = opt\$test,
  fitType = opt\$fit_type,
  minReplicatesForReplace = opt\$min_replicates_for_replace,
  useT = opt\$use_t,
  sfType = opt\$sf_type,
  parallel=TRUE, BPPARAM=MulticoreParam(opt\$cores)
)

# Variance-stabilized object for PCA
vsd_obj <- vst(
  dds,
  blind = opt[["vs_blind"]],
  nsub = opt[["vst_nsub"]]
)

# PCA data
pca_df <- plotPCA(
  vsd_obj,
  intgroup = contrast_variable,
  returnData = TRUE
)

percentVar <- round(100 * attr(pca_df, "percentVar"))

# Add sample names and generic grouping column
pca_df[["sample"]] <- rownames(pca_df)
pca_df[["group"]] <- pca_df[[contrast_variable]]

# PCA plot
PCA <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = group)
) +
  geom_point(size = 4) +
  geom_label_repel(
    aes(label = sample),
    size = 5,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  labs(color = NULL) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

pdf(
  file = paste(opt[["output_prefix"]], "deseq2.PCA.labelled.pdf", sep = "."),
  width = 8,
  height = 6
)
print(PCA)
dev.off()

png(
  file = paste(opt[["output_prefix"]], "deseq2.PCA.labelled.png", sep = "."),
  width = 800,
  height = 600,
  res = 120
)
print(PCA)
dev.off()

if (!is.null(opt\$contrast_string)) {
  coef_names <- resultsNames(dds)
  if (opt\$contrast_string %in% coef_names) {
    # Direct coefficient name
    comp.results <- run_results(name = opt\$contrast_string)
    if (opt\$shrink_lfc) {
      comp.results <- run_shrink(coef = opt\$contrast_string)
    }
  } else {
    # Parse as limma-style contrast expression
    design_mat <- model.matrix(as.formula(model), data = as.data.frame(colData(dds)))
    colnames(design_mat) <- make.names(colnames(design_mat))
    numeric_contrast <- as.numeric(
      limma::makeContrasts(contrasts = opt\$contrast_string, levels = colnames(design_mat))
    )

    # Run DESeq2 results with numeric contrast
    comp.results <- run_results(contrast = numeric_contrast)
    if (opt\$shrink_lfc) {
      comp.results <- run_shrink(contrast = numeric_contrast)
    }
  }
} else {
  contrast_var_tg_ref <- c(contrast_variable,
                           opt\$target_level,
                           opt\$reference_level)

  comp.results <- run_results(contrast = contrast_var_tg_ref)
  if (opt\$shrink_lfc) {
    coef_name <- paste0(
      contrast_variable,
      "_",
      make.names(opt[["target_level"]]),
      "_vs_",
      make.names(opt[["reference_level"]])
    )

    if (!coef_name %in% resultsNames(dds)) {
      stop(
        paste0(
          "Could not identify DESeq2 coefficient '",
          coef_name,
          "' required for apeglm shrinkage. Available coefficients: ",
          paste(resultsNames(dds), collapse = ", ")
        )
      )
    }

    comp.results <- run_shrink(coef = coef_name)
  }

  opt[["contrast_string"]] <- paste(
    opt[["target_level"]],
    opt[["reference_level"]],
    sep = "_vs_"
  )
}

comp.results <- data.frame(
  gene_id = rownames(comp.results),
  comp.results,
  check.names = FALSE
)

comp.results <- comp.results[order(comp.results[["padj"]], na.last = TRUE), ]

# Annotate gene names

if (opt[["organism"]] == "human") {
  
  library(org.Hs.eg.db)
  anno_db <- org.Hs.eg.db
  
} else if (opt[["organism"]] == "mouse") {
  
  library(org.Mm.eg.db)
  anno_db <- org.Mm.eg.db
  
} else {
  
  stop("Unsupported organism: ", opt[["organism"]])
  
}

gene_ids_for_annotation <- as.character(comp.results[["gene_id"]])
gene_ids_for_annotation <- sapply(
  strsplit(gene_ids_for_annotation, ".", fixed = TRUE),
  `[`,
  1
)

anno <- AnnotationDbi::select(
  anno_db,
  keys = gene_ids_for_annotation,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

anno <- anno[!duplicated(anno[["ENSEMBL"]]), ]

comp.results[["ENSEMBL"]] <- gene_ids_for_annotation

comp.results <- merge(
  comp.results,
  anno,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

# Significant DE genes at different absolute log2FC cutoffs
res_df <- comp.results

for (lfc_cutoff in c(0.5, 1, 1.5)) {
  
  sig_df <- res_df[
    !is.na(res_df[["padj"]]) &
      res_df[["padj"]] < 0.05 &
      !is.na(res_df[["log2FoldChange"]]) &
      abs(res_df[["log2FoldChange"]]) >= lfc_cutoff,
  ]
  
  sig_df <- sig_df[order(sig_df[["padj"]], na.last = TRUE), ]
  
  write.table(
    sig_df,
    file = paste(
      opt[["output_prefix"]],
      paste0("deseq2.sig_genes.log2FC", lfc_cutoff),
      "tsv",
      sep = "."
    ),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

# Volcano plot
volcano <- EnhancedVolcano(
  comp.results,
  lab = comp.results[["SYMBOL"]],
  selectLab = head(comp.results[["SYMBOL"]], 25),
  x = "log2FoldChange",
  y = "padj",
  
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Differential expression",
  subtitle = opt[["contrast_string"]],
  caption = NULL,
  
  pointSize = 2.5,
  labSize = 4,
  
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  
  legendPosition = "right"
)

pdf(
  file = paste(opt[["output_prefix"]], "deseq2.volcano.pdf", sep = "."),
  width = 10,
  height = 8
)
print(volcano)
dev.off()

png(
  file = paste(opt[["output_prefix"]], "deseq2.volcano.png", sep = "."),
  width = 1800,
  height = 1600,
  res = 200
)
print(volcano)
dev.off()

# MA plot

# MA plot

ma_df <- comp.results[
  !is.na(comp.results[["baseMean"]]) &
    !is.na(comp.results[["log2FoldChange"]]),
]

ma_df[["baseMean"]] <- as.numeric(ma_df[["baseMean"]])
ma_df[["log2FoldChange"]] <- as.numeric(ma_df[["log2FoldChange"]])
ma_df[["padj"]] <- as.numeric(ma_df[["padj"]])

ma_df[["significant"]] <- !is.na(ma_df[["padj"]]) & ma_df[["padj"]] < 0.05

ma_plot <- ggplot(
  ma_df,
  aes(
    x = baseMean,
    y = log2FoldChange,
    color = significant
  )
) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-5, 5)) +
  labs(
    title = paste("MA plot:", opt[["contrast_string"]]),
    x = "Mean normalized count",
    y = "log2 fold change",
    color = "padj < 0.05"
  ) +
  theme_bw(base_size = 14)

pdf(
  file = paste(opt[["output_prefix"]], "deseq2.MA.pdf", sep = "."),
  width = 7,
  height = 6
)
print(ma_plot)
dev.off()

png(
  file = paste(opt[["output_prefix"]], "deseq2.MA.png", sep = "."),
  width = 1200,
  height = 1000,
  res = 180
)
print(ma_plot)
dev.off()

# GSEA analysis

# Hallmark gene sets
if (opt[["organism"]] == "human") {
  msig_species <- "Homo sapiens"
} else if (opt[["organism"]] == "mouse") {
  msig_species <- "Mus musculus"
}

hallmark_sets <- msigdbr(
  species = msig_species,
  category = "H"
)

hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]

# Ranked vector for GSEA
if (! "stat" %in% colnames(comp.results)) {
  comp.results[["stat"]] <-
    as.numeric(comp.results[["log2FoldChange"]]) /
    as.numeric(comp.results[["lfcSE"]])
}

rank_df <- comp.results[
  !is.na(comp.results[["SYMBOL"]]) &
    comp.results[["SYMBOL"]] != "" &
    !is.na(comp.results[["stat"]]) &
    is.finite(comp.results[["stat"]]),
]

rank_df <- rank_df[order(abs(rank_df[["stat"]]), decreasing = TRUE), ]
rank_df <- rank_df[!duplicated(rank_df[["SYMBOL"]]), ]

gene_ranks <- rank_df[["stat"]]
names(gene_ranks) <- rank_df[["SYMBOL"]]
gene_ranks <- sort(gene_ranks, decreasing = TRUE)


hallmark_gsea <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark_term2gene,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  verbose = FALSE
)

hallmark_gsea_df <- as.data.frame(hallmark_gsea)

if (nrow(hallmark_gsea_df) > 0) {
  
  hallmark_gsea_df <- hallmark_gsea_df[
    order(
      hallmark_gsea_df[["p.adjust"]],
      -abs(hallmark_gsea_df[["NES"]]),
      na.last = TRUE
    ),
  ]
  
  hallmark_gsea_df[["pathway"]] <- gsub(
    "_",
    " ",
    sub("^HALLMARK_", "", hallmark_gsea_df[["Description"]])
  )
  
  write.table(
    hallmark_gsea_df,
    file = paste(opt[["output_prefix"]], "hallmark_gsea.tsv", sep = "."),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
  
  hallmark_dotplot <- dotplot(hallmark_gsea, showCategory = 20) +
    ggtitle("Hallmark GSEA")
  
  pdf(
    file = paste(opt[["output_prefix"]], "hallmark_gsea.dotplot.pdf", sep = "."),
    width = 10,
    height = 7
  )
  print(hallmark_dotplot)
  dev.off()
  
  png(
    file = paste(opt[["output_prefix"]], "hallmark_gsea.dotplot.png", sep = "."),
    width = 1800,
    height = 1300,
    res = 200
  )
  print(hallmark_dotplot)
  dev.off()
  
  # Barplot: top 20 pathways by absolute normalized enrichment score
  hallmark_barplot_df <- hallmark_gsea_df[
    order(abs(hallmark_gsea_df[["NES"]]), decreasing = TRUE),
  ]
  
  hallmark_barplot_df <- head(hallmark_barplot_df, 20)
  
  hallmark_barplot_df[["pathway"]] <- factor(
    hallmark_barplot_df[["pathway"]],
    levels = rev(hallmark_barplot_df[["pathway"]])
  )
  
  hallmark_barplot <- ggplot(
    hallmark_barplot_df,
    aes(
      x = NES,
      y = pathway,
      fill = p.adjust
    )
  ) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Hallmark GSEA",
      x = "Normalized enrichment score",
      y = NULL,
      fill = "Adjusted p-value"
    ) +
    theme_bw(base_size = 14)
  
  pdf(
    file = paste(
      opt[["output_prefix"]],
      "hallmark_gsea.barplot.pdf",
      sep = "."
    ),
    width = 9,
    height = 7
  )
  print(hallmark_barplot)
  dev.off()
  
  png(
    file = paste(
      opt[["output_prefix"]],
      "hallmark_gsea.barplot.png",
      sep = "."
    ),
    width = 1600,
    height = 1300,
    res = 200
  )
  print(hallmark_barplot)
  dev.off()
  
} else {
  
  warning("No Hallmark GSEA terms returned; skipping GSEA plots.")
  
  write.table(
    hallmark_gsea_df,
    file = paste(opt[["output_prefix"]], "hallmark_gsea.tsv", sep = "."),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

# GO over-representation analysis on significant genes

go_gene_df <- comp.results[
  !is.na(comp.results[["ENTREZID"]]) &
    comp.results[["ENTREZID"]] != "" &
    !is.na(comp.results[["padj"]]) &
    comp.results[["padj"]] < 0.05 &
    !is.na(comp.results[["log2FoldChange"]]) &
    abs(comp.results[["log2FoldChange"]]) >= 1,
]

go_genes <- unique(as.character(go_gene_df[["ENTREZID"]]))

go_universe_df <- comp.results[
  !is.na(comp.results[["ENTREZID"]]) &
    comp.results[["ENTREZID"]] != "",
]

go_universe <- unique(as.character(go_universe_df[["ENTREZID"]]))

go_ontologies <- c("BP", "MF", "CC")

for (go_ont in go_ontologies) {
  
  if (length(go_genes) >= 5) {
    
    go_res <- enrichGO(
      gene = go_genes,
      universe = go_universe,
      OrgDb = anno_db,
      keyType = "ENTREZID",
      ont = go_ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 1000,
      readable = TRUE
    )
    
    go_df <- as.data.frame(go_res)
    
  } else {
    
    warning("Fewer than 5 significant genes for GO ", go_ont, "; writing empty table.")
    go_df <- data.frame()
    go_res <- NULL
  }
  
  if (nrow(go_df) > 0) {
    
    go_df <- go_df[
      order(
        go_df[["p.adjust"]],
        go_df[["pvalue"]],
        na.last = TRUE
      ),
    ]
    
    write.table(
      go_df,
      file = paste(opt[["output_prefix"]], paste0("go_", go_ont, "_ora.tsv"), sep = "."),
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    go_dotplot <- dotplot(go_res, showCategory = 20) +
      ggtitle(paste0("GO ", go_ont, " over-representation"))
    
    pdf(
      file = paste(opt[["output_prefix"]], paste0("go_", go_ont, "_ora.dotplot.pdf"), sep = "."),
      width = 10,
      height = 7
    )
    print(go_dotplot)
    dev.off()
    
    png(
      file = paste(opt[["output_prefix"]], paste0("go_", go_ont, "_ora.dotplot.png"), sep = "."),
      width = 1800,
      height = 1300,
      res = 200
    )
    print(go_dotplot)
    dev.off()
    
    go_barplot <- barplot(go_res, showCategory = 20) +
      ggtitle(paste0("GO ", go_ont, " over-representation"))
    
    pdf(
      file = paste(opt[["output_prefix"]], paste0("go_", go_ont, "_ora.barplot.pdf"), sep = "."),
      width = 9,
      height = 7
    )
    print(go_barplot)
    dev.off()
    
    png(
      file = paste(opt[["output_prefix"]], paste0("go_", go_ont, "_ora.barplot.png"), sep = "."),
      width = 1600,
      height = 1300,
      res = 200
    )
    print(go_barplot)
    dev.off()
    
  } else {
    
    write.table(
      go_df,
      file = paste(opt[["output_prefix"]], paste0("go_", go_ont, "_ora.tsv"), sep = "."),
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
  }
}

# See https://support.bioconductor.org/p/97676/

if (!is.null(opt\$transcript_lengths_file)){
  size_factors = estimateSizeFactorsForMatrix(counts(dds) / assays(dds)[["avgTxLength"]])
} else {
  size_factors = sizeFactors(dds)
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# if no explicit contrast_string was given, synthesize a two‐level default
if (! is_valid_string(opt\$contrast_string)) {
  opt\$contrast_string <- paste(opt\$target_level, opt\$reference_level, sep = "_vs_")
}

cat("Saving results for ", opt\$contrast_string, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results

if (! is.null(opt\$round_digits)){
  comp.results <- apply(data.frame(comp.results), 2, function(x) round(x, opt\$round_digits))
}
comp.results <- `colnames<-`(
  data.frame(
    gene_id = rownames(comp.results),
    comp.results,
    check.names = FALSE
  ),
  c(opt\$gene_id_col, colnames(comp.results))  # Setting all column names
)

write.table(
  comp.results,
  file = paste(opt\$output_prefix, 'deseq2.results.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Dispersion plot

png(
  file = paste(opt\$output_prefix, 'deseq2.dispersion.png', sep = '.'),
  width = 600,
  height = 600
)
plotDispEsts(dds)
dev.off()

pdf(
  file = paste(opt\$output_prefix, 'deseq2.dispersion.pdf', sep = '.')
)
plotDispEsts(dds)
dev.off()

# R object for other processes to use

saveRDS(dds, file = paste(opt\$output_prefix, 'dds.rld.rds', sep = '.'))

# Size factors

sf_df = data.frame(
  sample = names(size_factors),
  data.frame(size_factors, check.names = FALSE),
  check.names = FALSE
)
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
  sf_df,
  file = paste(opt\$output_prefix, 'deseq2.sizefactors.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Write specified matrices

normalised_matrix <- counts(dds, normalized = TRUE)
if (! is.null(opt\$round_digits)){
  normalised_matrix <- apply(normalised_matrix, 2, function(x) round(x, opt\$round_digits))
}
normalised_matrix <- `colnames<-`(
  data.frame(
    gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
    normalised_matrix,                # Other columns
    check.names = FALSE
  ),
  c(opt\$gene_id_col, colnames(normalised_matrix))  # Setting all column names
)

write.table(
  normalised_matrix,
  file = paste(opt\$output_prefix, 'normalised_counts.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# Note very limited rounding for consistency of results

for (vs_method_name in strsplit(opt[["vs_method"]], ",")[[1]]) {
  if (vs_method_name == "vst") {
    vs_mat <- assay(
      vst(
        dds,
        blind = opt[["vs_blind"]],
        nsub = opt[["vst_nsub"]]
      )
    )
    
  } else if (vs_method_name == "rlog") {
    vs_mat <- assay(
      rlog(
        dds,
        blind = opt[["vs_blind"]],
        fitType = opt[["fit_type"]]
      )
    )
  }
  
  if (!is.null(opt[["round_digits"]])) {
    vs_mat <- apply(vs_mat, 2, function(x) round(x, opt[["round_digits"]]))
  }
  
  vs_mat <- `colnames<-`(
    data.frame(
      gene_id = rownames(counts(dds)),
      vs_mat,
      check.names = FALSE
    ),
    c(opt[["gene_id_col"]], colnames(vs_mat))
  )
  
  write.table(
    vs_mat,
    file = paste(opt[["output_prefix"]], vs_method_name, "tsv", sep = "."),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

# Save model to file

write(model, file=paste(opt\$output_prefix, 'deseq2.model.txt', sep = '.'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r.version),
    paste('    bioconductor-deseq2:', deseq2.version)
  ),
  'versions.yml')

################################################
################################################
## RMARKDOWN REPORT                           ##
################################################
################################################


rmarkdown::render(
  input = "report_template.Rmd",
  output_file = paste(opt[["output_prefix"]], "deseq2.report.html", sep = "."),
  params = list(prefix = opt[["output_prefix"]]),
  quiet = FALSE
)