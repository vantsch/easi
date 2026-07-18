#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop(
        "Usage: isa_part2.R <switchlist_part1.rds> <pfam_isa.tsv>",
        call. = FALSE
    )
}

switchlist_file <- args[1]
pfam_file       <- args[2]

suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
})

message("Reading ISA Part 1 object: ", switchlist_file)
SwitchList <- readRDS(switchlist_file)

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
    quiet = TRUE
)

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

message("ISA Part 2 completed")