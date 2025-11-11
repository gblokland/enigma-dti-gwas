#!/usr/bin/env Rscript
# ENIGMA-DTI QC Combined Script (complete)
# Author: Gabriella Blokland (merged & fixed)
# Last update: 2025-10-31
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# ---------------- Setup / packages ----------------
needed <- c("argparse", "ggplot2", "dplyr", "tidyr", "stringr", "MBESS", "readr")
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
})

# ---------------- Argument parser ----------------
parser <- ArgumentParser(description = "ENIGMA-DTI QC pipeline - combined")
parser$add_argument("--cohort", required = TRUE, help = "Cohort name")
parser$add_argument("--covarFILE", required = TRUE, help = "Covariate file (.covar, tab or csv)")
parser$add_argument("--phenoFILE", required = TRUE, help = "Phenotype file (.pheno, tab or csv)")
parser$add_argument("--icvFILE", required = TRUE, help = "ICV file (.csv)")
parser$add_argument("--icvColumnHeader", required = TRUE, default = "ICV", help = "Name of ICV column")
parser$add_argument("--ageColumnHeader", default = "Age", help = "Name of Age column")
parser$add_argument("--sexColumnHeader", default = "Sex", help = "Name of Sex column")
parser$add_argument("--maleIndicator", default = "1", help = "Value used for male")
parser$add_argument("--CaseControlCohort", default = "1", help = "CaseControlCohort yes(1) or no(0)")
parser$add_argument("--affectedStatusColumnHeader", default = "AffectionStatus", help = "AffectionStatus column name")
parser$add_argument("--affectedIndicator", default = "2", help = "Value indicating affected/case")
parser$add_argument("--related", default = "0", help = "Related cohort yes(1) or no(0)")
parser$add_argument("--rois", default = "GlobalAverage;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R", help = "Semicolon-separated list of ROIs (use same naming as pheno columns)")
parser$add_argument("--pheno_covar_dir", default = "./QC_ENIGMA/", help = "pheno/covar dir (not used heavily)")
parser$add_argument("--outDir", default = "./QC_ENIGMA/", help = "Output directory")
parser$add_argument("--outPDF", default = "ENIGMA_DTI_allROI_histograms.pdf", help = "output PDF name (multi hist)")
parser$add_argument("--outTXT", default = "ENIGMA_DTI_allROI_stats.txt", help = "output stats txt")
parser$add_argument("--eName", default = "ENIGMA_DTI_GWAS", help = "ENIGMA label")

args <- parser$parse_args()
cat("Inputs: \n")
print(args)

# assign
cohort <- args$cohort
covarFILE <- args$covarFILE
phenoFILE <- args$phenoFILE
icvFILE <- args$icvFILE
icvColumnHeader <- args$icvColumnHeader
ageColumnHeader <- args$ageColumnHeader
sexColumnHeader <- args$sexColumnHeader
maleIndicator <- args$maleIndicator
CaseControlCohort <- args$CaseControlCohort
affectedStatusColumnHeader <- args$affectedStatusColumnHeader
affectedIndicator <- args$affectedIndicator
related <- args$related
rois <- args$rois
pheno_covar_dir <- args$pheno_covar_dir
outDir <- args$outDir
outPDF <- args$outPDF
outTXT <- args$outTXT
eName <- args$eName

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
message("Output dir: ", normalizePath(outDir))

# ---------------- Read files (robust) ----------------
read_table_auto <- function(path) {
  # try read.table then readr fallback
  dat <- tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE, sep = "", comment.char = "", quote = "\"'"),
                  error = function(e) NULL)
  if (is.null(dat)) {
    dat <- tryCatch(readr::read_delim(path, delim = "\t", col_types = readr::cols(.default = "c")),
                    error = function(e) readr::read_csv(path, col_types = readr::cols(.default = "c")))
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  }
  return(dat)
}

message("Reading covariates from: ", covarFILE)
covar <- read_table_auto(covarFILE)
message("Reading phenotypes from: ", phenoFILE)
pheno <- read_table_auto(phenoFILE)

# Remove average columns from covar if present
covar <- covar[ , !(names(covar) %in% c("FA_GlobalAverage", "MD_GlobalAverage", "RD_GlobalAverage", "AD_GlobalAverage"))]

# replace literal "x" or "X" with NA across both tables for safety
for (tbl in list(covar, pheno)) {
  for (col in names(tbl)) {
    if (is.character(tbl[[col]])) {
      tbl[[col]][tbl[[col]] %in% c("x", "X")] <- NA
    }
  }
}
# put back (we mutated copies above, so reassign properly)
# (we re-read covar/pheno, easier to re-run replacement on original objects)
for (col in names(covar)) if (is.character(covar[[col]])) covar[[col]][covar[[col]] %in% c("x","X")] <- NA
for (col in names(pheno)) if (is.character(pheno[[col]])) pheno[[col]][pheno[[col]] %in% c("x","X")] <- NA

# Merge covar and pheno by FID/IID (keep all)
if (!all(c("FID","IID") %in% names(covar)) || !all(c("FID","IID") %in% names(pheno))) {
  stop("FID and IID must be columns in both covar and pheno files.")
}
Table <- merge(covar, pheno, by = c("FID", "IID"), all = TRUE)
message("Merged table rows: ", nrow(Table))

if (!is.na(icvFILE)) {
  icv <- read.csv("ICV.csv", header = TRUE)
}

# Merge ICV by FID/IID (keep all)
if (!all(c("FID","IID") %in% names(icv)) || !all(c("FID","IID") %in% names(Table))) {
  stop("FID and IID must be columns in both covar and pheno files.")
}
Table <- merge(icv, Table, by = c("FID", "IID"), all = TRUE)
message("Merged table rows: ", nrow(Table))
Table$ICV <- Table[,icvColumnHeader]

# Normalize group variable
Table$AffectionStatus <- trimws(as.character(Table$AffectionStatus))
Table$AffectionStatus <- tolower(Table$AffectionStatus)
Table$AffectionStatus <- ifelse(Table$AffectionStatus %in% c("affected", "1", "case"),
                                "Affected",
                                ifelse(Table$AffectionStatus %in% c("control", "0", "healthy"),
                                       "Control", NA))
Table$AffectionStatus <- factor(Table$AffectionStatus, levels = c("Control", "Affected"))

# Check
table(Table$AffectionStatus)


# ---------------- ROIs and labels ----------------
parsedROIs <- unlist(strsplit(rois, ";"))
parsedROIs <- parsedROIs[parsedROIs != ""]
Nrois <- length(parsedROIs)
message("Nrois = ", Nrois)
message("ROIs: ", paste(parsedROIs, collapse = ", "))

# var labels & colors
roiLabels <- gsub("\\.L$", " Left", gsub("\\.R$", " Right", parsedROIs))
# produce color vector
library(grDevices)
roiColors <- setNames(rainbow(Nrois), parsedROIs)

# ---------------- Standardize key covariates ----------------
# AffectionStatus: map affectedIndicator -> "Affected", others-> "Control"
if (affectedStatusColumnHeader %in% names(Table)) {
  Table[[affectedStatusColumnHeader]] <- as.character(Table[[affectedStatusColumnHeader]])
  Table[[affectedStatusColumnHeader]][Table[[affectedStatusColumnHeader]] == as.character(affectedIndicator)] <- "Affected"
  Table[[affectedStatusColumnHeader]][Table[[affectedStatusColumnHeader]] != "Affected"] <- "Control"
  Table[[affectedStatusColumnHeader]] <- factor(Table[[affectedStatusColumnHeader]], levels = c("Affected","Control"))
} else {
  message("Warning: affected status column '", affectedStatusColumnHeader, "' not found in Table.")
}

# Sex: try to standardize to "Male"/"Female"
if (sexColumnHeader %in% names(Table)) {
  sex_vals <- unique(na.omit(Table[[sexColumnHeader]]))
  # attempt familiar codings
  Table[[sexColumnHeader]] <- as.character(Table[[sexColumnHeader]])
  if (all(sex_vals %in% c("Male","Female","M","F","1","2","0"))) {
    # map numeric codes if needed (common: 1=Male, 2=Female; or 0=Female,1=Male)
    Table[[sexColumnHeader]][Table[[sexColumnHeader]] %in% c("M","m","male")] <- "Male"
    Table[[sexColumnHeader]][Table[[sexColumnHeader]] %in% c("F","f","female")] <- "Female"
    # numeric strings
    if (all(sex_vals %in% c("1","2"))) {
      # assume maleIndicator provided -> maleIndicator value means male
      Table[[sexColumnHeader]][Table[[sexColumnHeader]] == as.character(maleIndicator)] <- "Male"
      Table[[sexColumnHeader]][Table[[sexColumnHeader]] != "Male"] <- "Female"
    } else if (all(sex_vals %in% c("0","1"))) {
      # common 0=Female,1=Male
      Table[[sexColumnHeader]][Table[[sexColumnHeader]] == "1"] <- "Male"
      Table[[sexColumnHeader]][Table[[sexColumnHeader]] == "0"] <- "Female"
    }
    Table[[sexColumnHeader]] <- factor(Table[[sexColumnHeader]], levels = c("Male","Female"))
  } else {
    message("Unusual sex values detected: ", paste(sex_vals, collapse = ", "), " — left as-is.")
  }
} else {
  message("Sex column '", sexColumnHeader, "' not found.")
}

# ---------------- Build long DTI table ----------------
dti_cols <- grep("^(FA|MD|AD|RD)_", names(Table), value = TRUE)
if (length(dti_cols) == 0) {
  warning("No DTI columns found with prefixes FA_, MD_, AD_, or RD_.")
  TableLong <- data.frame() # empty
} else {
  TableLong <- Table %>%
    pivot_longer(cols = all_of(dti_cols), names_to = "Variable", values_to = "Value")
  # parse metric & tract & side
  TableLong <- TableLong %>%
    mutate(
      Metric = str_extract(Variable, "^(FA|MD|AD|RD)"),
      TractSide = stringr::str_remove(Variable, "^(FA|MD|AD|RD)_"),
      Side = case_when(
        str_detect(TractSide, "\\.L$") ~ "L",
        str_detect(TractSide, "\\.R$") ~ "R",
        TRUE ~ NA_character_
      ),
      Tract = str_remove(TractSide, "\\.L$|\\.R$"),
      parsedROI = ifelse(is.na(Side), Tract, paste0(Tract, ".", Side))
    )
  message("TableLong rows: ", nrow(TableLong), " (DTI values long)")
}

# ---------------- Numeric ranges for plotting ----------------
FA_range <- c(0, 1)
MD_range <- c(0.5e-3, 1.0e-3)
RD_range <- c(0.3e-3, 0.8e-3)
AD_range <- c(1.0e-3, 1.75e-3)
range_map <- list(FA = FA_range, MD = MD_range, RD = RD_range, AD = AD_range)

# ---------------- Functions: plotting & stats ----------------
safe_ggsave <- function(plot, filename, width = 36, height = 24) {
  ggsave(filename = filename, plot = plot, width = width, height = height, units = "cm", dpi = 600)
  message("Saved: ", filename)
}

plot_multi_histogram <- function(Table, metric, outdir = outDir, ncol = 8) {
  metric <- toupper(metric)
  metric_long <- TableLong %>% filter(Metric == metric)
  if (nrow(metric_long) == 0) {
    message("No data for metric: ", metric)
    return(invisible(NULL))
  }
  # set x limits if defined in range_map
  x_limits <- range_map[[metric]]
  p <- ggplot(metric_long, aes(x = Value, fill = .data[[affectedStatusColumnHeader]])) +
    geom_histogram(aes(y = after_stat(count)), colour = "black", bins = 30) +
    facet_wrap(~parsedROI, ncol = ncol, scales = "free_x") +
    theme_bw() +
    labs(x = metric, y = "Number of Subjects", fill = "AffectionStatus", title = paste0(cohort, " - ", metric))
  if (!is.null(x_limits)) p <- p + xlim(x_limits)
  fn <- file.path(outdir, paste0(cohort, "_histogram_multi-panel_", metric, "_by_AffectionStatus.pdf"))
  safe_ggsave(p, fn)
}

plot_age_histograms <- function(Table, outdir = outDir) {
  if (!("Age" %in% names(Table))) {
    message("No Age column found.")
    return(invisible(NULL))
  }
  p_all <- ggplot(Table, aes(x = Age, fill = .data[[affectedStatusColumnHeader]])) +
    geom_histogram(aes(y = after_stat(count)), colour = "black", bins = 30) +
    theme_bw() + labs(x = "Age (Years)", y = "Number of Subjects", title = paste0(cohort, " - Age"))
  safe_ggsave(p_all, file.path(outdir, paste0(cohort, "_histogram_Age.pdf")))
}

plot_icv_checks <- function(Table, outdir = outDir) {
  if (!("ICV" %in% names(Table))) {
    message("ICV not present - skipping ICV checks.")
    return(invisible(NULL))
  }
  # histogram
  p_hist <- ggplot(Table, aes(x = ICV, fill = .data[[sexColumnHeader]])) +
    geom_histogram(aes(y = after_stat(count)), bins = 30, colour = "black") +
    theme_bw() + labs(title = paste0(cohort, " - ICV distribution by Sex"))
  safe_ggsave(p_hist, file.path(outdir, paste0(cohort, "_ICV_hist_by_Sex.pdf")), width = 20, height = 15)
  # boxplot
  if (sexColumnHeader %in% names(Table)) {
    p_box <- ggplot(Table, aes(x = .data[[sexColumnHeader]], y = ICV, fill = .data[[sexColumnHeader]])) +
      geom_boxplot() + theme_bw() + labs(title = paste0(cohort, " - ICV by Sex"))
    safe_ggsave(p_box, file.path(outdir, paste0(cohort, "_ICV_box_by_Sex.pdf")), width = 12, height = 9)
  }
}

# ---------------- Summary table: mean ± SD by metric x group ----------------
compute_summary_table <- function(TableLong, outdir = outDir) {
  if (nrow(TableLong) == 0) {
    warning("No long DTI data to summarise.")
    return(data.frame())
  }
  summary_table <- TableLong %>%
    group_by(Metric, parsedROI, !!sym(affectedStatusColumnHeader)) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      N = sum(!is.na(Value)),
      .groups = "drop"
    ) %>%
    mutate(Mean_SD = sprintf("%.6g \u00B1 %.6g", Mean, SD)) %>%
    rename(AffectionStatus = !!sym(affectedStatusColumnHeader)) %>%
    arrange(Metric, parsedROI, AffectionStatus)
  # save
  out_csv <- file.path(outdir, paste0(cohort, "_Summary_MeanSD_byMetricGroup.csv"))
  write.csv(summary_table, out_csv, row.names = FALSE)
  message("Saved summary table: ", out_csv)
  return(summary_table)
}

# ---------------- Cohen's d per metric x ROI ----------------
compute_cohens_d <- function(Table, parsedROIs, roiLabels, roiColors, outdir = outDir) {
  results <- data.frame()
  # counts
  if (!(affectedStatusColumnHeader %in% names(Table))) {
    message("Affection status column missing. Skipping Cohen's d.")
    return(results)
  }
  n_CON <- sum(Table[[affectedStatusColumnHeader]] == "Control", na.rm = TRUE)
  n_AFF <- sum(Table[[affectedStatusColumnHeader]] == "Affected", na.rm = TRUE)
  if (n_CON <= 5 || n_AFF <= 5) {
    message("Not enough subjects per group for Cohen's d (n_CON=", n_CON, ", n_AFF=", n_AFF, ").")
    return(results)
  }

  for (metric in c("FA","MD","AD","RD")) {
    for (i in seq_along(parsedROIs)) {
      colname <- paste0(metric, "_", parsedROIs[i])
      if (!(colname %in% names(Table))) {
        # skip if missing
        next
      }
      x_aff <- as.numeric(Table[[colname]][Table[[affectedStatusColumnHeader]] == "Affected"])
      x_con <- as.numeric(Table[[colname]][Table[[affectedStatusColumnHeader]] == "Control"])
      x_aff <- x_aff[!is.na(x_aff)]; x_con <- x_con[!is.na(x_con)]
      if (length(x_aff) < 2 || length(x_con) < 2) {
        # insufficient data
        row <- data.frame(metric = metric,
                          variable = colname,
                          varlabel = roiLabels[i],
                          colour = roiColors[i],
                          cohens_d_smd = NA_real_,
                          cohens_d_cilb = NA_real_,
                          cohens_d_ciub = NA_real_,
                          cohens_d_se = NA_real_,
                          N_aff = length(x_aff),
                          N_con = length(x_con),
                          stringsAsFactors = FALSE)
        results <- bind_rows(results, row)
        next
      }
      # pooled sd
      pooled_sd <- sqrt((sum((x_aff - mean(x_aff))^2) + sum((x_con - mean(x_con))^2)) / (length(x_aff) + length(x_con) - 2))
      d <- (mean(x_aff) - mean(x_con)) / pooled_sd
      # t-test and MBESS CI (try-catch if MBESS fails)
      tt <- tryCatch(t.test(x_aff, x_con), error = function(e) NULL)
      ci <- tryCatch({
        if (!is.null(tt)) MBESS::ci.smd(ncp = tt$statistic, n.1 = length(x_aff), n.2 = length(x_con), conf.level = 0.95)
        else list(smd = d, Lower.Conf.Limit.smd = NA, Upper.Conf.Limit.smd = NA)
      }, error = function(e) list(smd = d, Lower.Conf.Limit.smd = NA, Upper.Conf.Limit.smd = NA))
      se <- if (!is.na(ci$Upper.Conf.Limit.smd) && !is.na(ci$Lower.Conf.Limit.smd)) (ci$Upper.Conf.Limit.smd - ci$Lower.Conf.Limit.smd) / 3.92 else NA_real_
      row <- data.frame(metric = metric,
                        variable = colname,
                        varlabel = roiLabels[i],
                        colour = roiColors[i],
                        cohens_d_smd = ci$smd,
                        cohens_d_cilb = ci$Lower.Conf.Limit.smd,
                        cohens_d_ciub = ci$Upper.Conf.Limit.smd,
                        cohens_d_se = se,
                        N_aff = length(x_aff),
                        N_con = length(x_con),
                        stringsAsFactors = FALSE)
      results <- bind_rows(results, row)
    } # end roi
  } # end metric

  # Save results
  out_csv <- file.path(outdir, paste0(cohort, "_cohensd_results.csv"))
  write.csv(results, out_csv, row.names = FALSE)
  message("Saved Cohen's d results: ", out_csv)
  return(results)
}

# ---------------- Run the steps ----------------
message("Running QC pipeline...")

# Age histograms and group breakdowns
plot_age_histograms(Table, outDir)

# ICV checks
plot_icv_checks(Table, outDir)

# Per-metric multi-panel histograms
for (m in c("FA","MD","AD","RD")) {
  plot_multi_histogram(Table, m, outdir = outDir)
}

# Build & save summary table from long data
summary_table <- compute_summary_table(TableLong, outDir)

# Compute Cohen's d and produce plot if results exist
cohens_d_results <- compute_cohens_d(Table, parsedROIs, roiLabels, roiColors, outDir)

if (nrow(cohens_d_results) > 0) {
  # ensure ordering & factor levels
  cohens_d_results <- cohens_d_results %>%
    mutate(varlabel = factor(varlabel, levels = unique(varlabel))) %>%
    arrange(metric, varlabel)

  # Plot Cohen's d (faceted by metric)
  p_cd <- ggplot(cohens_d_results, aes(x = varlabel, y = cohens_d_smd, fill = metric)) +
    geom_col(colour = "black") +
    geom_errorbar(aes(ymin = cohens_d_smd - cohens_d_se, ymax = cohens_d_smd + cohens_d_se), width = 0.25) +
    facet_wrap(~metric, ncol = 2, scales = "free_y") +
    coord_flip() +
    geom_hline(yintercept = 0) +
    labs(x = "", y = "Cohen's d ± SE", title = paste0(cohort, ": Cohen's d (AFF vs CON)")) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"), axis.text.y = element_text(size = 6))
  safe_ggsave(p_cd, file.path(outDir, paste0(cohort, "_cohensd_bar_graph_AFF-CON.png")), width = 18, height = 24)
} else {
  message("No Cohen's d results to plot.")
}

message("QC pipeline finished. Outputs in: ", normalizePath(outDir))
