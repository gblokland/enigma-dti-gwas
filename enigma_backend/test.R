#!/usr/bin/env Rscript
# ENIGMA-DTI QC Combined Script (complete runnable version)
# Author: Gabriëlla Blokland
# Last update: 2025-11-11

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# ---------------- Setup / packages ----------------
needed <- c("argparse", "ggplot2", "dplyr", "tidyr", "stringr", "MBESS", "readr", "grDevices")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(MBESS)
  library(readr)
  library(grDevices)
})

# ---------------- Argument parser ----------------
parser <- ArgumentParser(description = "ENIGMA-DTI QC pipeline - combined")
parser$add_argument("--cohort", required = TRUE)
parser$add_argument("--covarFILE", required = TRUE)
parser$add_argument("--phenoFILE", required = TRUE)
parser$add_argument("--icvFILE", default = "ICV.csv")
parser$add_argument("--icvColumnHeader", default = "ICV")
parser$add_argument("--ageColumnHeader", default = "Age")
parser$add_argument("--sexColumnHeader", default = "Sex")
parser$add_argument("--maleIndicator", default = "1")
parser$add_argument("--CaseControlCohort", default = "1")
parser$add_argument("--affectedStatusColumnHeader", default = "AffectionStatus")
parser$add_argument("--affectedIndicator", default = "2")
parser$add_argument("--related", default = "0")
parser$add_argument("--rois", default = "GlobalAverage;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC")
parser$add_argument("--outDir", default = "./QC_ENIGMA/")
parser$add_argument("--outPDF", default = "ENIGMA_DTI_allROI_histograms.pdf")
parser$add_argument("--outTXT", default = "ENIGMA_DTI_allROI_stats.txt")
parser$add_argument("--eName", default = "ENIGMA_DTI_GWAS")

args <- parser$parse_args()
print(args)

# ---------------- Prepare ----------------
cohort <- args$cohort
covarFILE <- args$covarFILE
phenoFILE <- args$phenoFILE
icvFILE <- args$icvFILE
icvColumnHeader <- args$icvColumnHeader
ageColumnHeader <- args$ageColumnHeader
sexColumnHeader <- args$sexColumnHeader
maleIndicator <- args$maleIndicator
affectedStatusColumnHeader <- args$affectedStatusColumnHeader
affectedIndicator <- args$affectedIndicator
outDir <- args$outDir
outPDF <- args$outPDF
outTXT <- args$outTXT
eName <- args$eName
rois <- args$rois

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
message("Output dir: ", normalizePath(outDir))

# ---------------- Read data ----------------
read_table_auto <- function(path) {
  dat <- tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (is.null(dat)) {
    dat <- tryCatch(read_delim(path, delim = "\t", col_types = cols(.default = "c")),
                    error = function(e) read_csv(path, col_types = cols(.default = "c")))
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  }
  dat
}

covar <- read_table_auto(covarFILE)
pheno <- read_table_auto(phenoFILE)

if (!all(c("FID","IID") %in% names(covar)) || !all(c("FID","IID") %in% names(pheno))) {
  stop("Both covar and pheno must have FID and IID columns.")
}
Table <- merge(covar, pheno, by = c("FID", "IID"), all = TRUE)
if (file.exists(icvFILE)) {
  icv <- read.csv(icvFILE, header = TRUE)
  Table <- merge(Table, icv, by = c("FID", "IID"), all = TRUE)
  Table$ICV <- Table[[icvColumnHeader]]
}

# ---------------- Fix AffectionStatus ----------------
if (affectedStatusColumnHeader %in% names(Table)) {
  Table[[affectedStatusColumnHeader]] <- ifelse(Table[[affectedStatusColumnHeader]] == affectedIndicator, "Affected", "Control")
  Table[[affectedStatusColumnHeader]] <- factor(Table[[affectedStatusColumnHeader]], levels = c("Control", "Affected"))
} else {
  message("No affected status column found, skipping case-control coding.")
}

# ---------------- Extract DTI columns ----------------
dti_cols <- grep("^(FA|MD|AD|RD)_", names(Table), value = TRUE)
if (length(dti_cols) == 0) stop("No DTI columns found with FA_/MD_/AD_/RD_ prefix.")

TableLong <- Table %>%
  pivot_longer(cols = all_of(dti_cols), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Metric = str_extract(Variable, "^(FA|MD|AD|RD)"),
    ROI = str_remove(Variable, "^(FA|MD|AD|RD)_")
  )

# ---------------- Plot helpers ----------------
safe_ggsave <- function(plot, filename, width = 36, height = 24) {
  tryCatch({
    ggsave(filename = filename, plot = plot, width = width, height = height, units = "cm", dpi = 600)
    message("Saved: ", filename)
  }, error = function(e) message("Failed to save plot ", filename, ": ", e$message))
}

# ---------------- Plot per metric ----------------
for (metric in c("FA","MD","AD","RD")) {
  subdat <- TableLong %>% filter(Metric == metric)
  if (nrow(subdat) == 0) next
  p <- ggplot(subdat, aes(x = Value, fill = Table[[affectedStatusColumnHeader]])) +
    geom_histogram(bins = 30, colour = "black") +
    facet_wrap(~ROI, scales = "free_x") +
    theme_bw() +
    labs(title = paste0(cohort, " - ", metric), x = metric, y = "Count")
  safe_ggsave(p, file.path(outDir, paste0(cohort, "_", metric, "_histograms.pdf")))
}

# ---------------- Compute summary stats ----------------
summary_table <- TableLong %>%
  group_by(Metric, ROI, !!sym(affectedStatusColumnHeader)) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), N = sum(!is.na(Value)), .groups = "drop")
write.csv(summary_table, file.path(outDir, paste0(cohort, "_SummaryStats.csv")), row.names = FALSE)
message("Saved summary stats")

# ---------------- Compute Cohen's d ----------------
compute_cohens_d <- function(Table, metric, roi) {
  col <- paste0(metric, "_", roi)
  if (!(col %in% names(Table))) return(NA)
  x_aff <- as.numeric(Table[[col]][Table[[affectedStatusColumnHeader]] == "Affected"])
  x_con <- as.numeric(Table[[col]][Table[[affectedStatusColumnHeader]] == "Control"])
  if (length(x_aff) < 3 || length(x_con) < 3) return(NA)
  pooled_sd <- sqrt(((length(x_aff) - 1)*var(x_aff) + (length(x_con) - 1)*var(x_con)) / (length(x_aff)+length(x_con)-2))
  (mean(x_aff) - mean(x_con)) / pooled_sd
}

parsedROIs <- unlist(strsplit(rois, ";"))
results <- expand.grid(Metric = c("FA","MD","AD","RD"), ROI = parsedROIs, stringsAsFactors = FALSE)
results$d <- mapply(compute_cohens_d, Table = list(Table), metric = results$Metric, roi = results$ROI)
write.csv(results, file.path(outDir, paste0(cohort, "_CohensD.csv")), row.names = FALSE)
message("Saved Cohen's d values")

# ---------------- Done ----------------
message("QC pipeline completed successfully for cohort: ", cohort)
q(save = "no")
