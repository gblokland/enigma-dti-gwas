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
parser$add_argument("--icvFILE", required = TRUE, default = "ICV.csv", help = "ICV file (.csv)")
parser$add_argument("--icvColumnHeader", required = TRUE, default = "ICV", help = "Name of ICV column")
parser$add_argument("--ageColumnHeader", default = "Age", help = "Name of Age column")
parser$add_argument("--sexColumnHeader", default = "Sex", help = "Name of Sex column")
parser$add_argument("--maleIndicator", default = "1", help = "Value used for male")
parser$add_argument("--CaseControlCohort", default = "1", help = "CaseControlCohort yes(1) or no(0)")
parser$add_argument("--affectedStatusColumnHeader", default = "AffectionStatus", help = "AffectionStatus column name")
parser$add_argument("--affectedIndicator", default = "2", help = "Value indicating affected/case")
parser$add_argument("--related", default = "0", help = "Related cohort yes(1) or no(0)")
parser$add_argument("--rois", default = "GlobalAverage;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R", help = "Semicolon-separated list of ROIs (use same naming as pheno columns)")
parser$add_argument("--outDir", default = "./QC_ENIGMA/", help = "Output directory")
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
outDir <- args$outDir
eName <- args$eName

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
}
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
  icv <- read.csv(icvFILE, header = TRUE)
}

# Merge ICV by FID/IID (keep all)
if (!all(c("FID","IID") %in% names(icv)) || !all(c("FID","IID") %in% names(Table))) {
  stop("FID and IID must be columns in both covar and pheno files.")
}
Table <- merge(icv, Table, by = c("FID", "IID"), all = TRUE)
message("Merged table rows: ", nrow(Table))
if (icvColumnHeader %in% names(Table)) {
  Table$ICV <- Table[,icvColumnHeader]
}
# Normalize group variable
Table$AffectionStatus <- trimws(as.character(Table$AffectionStatus))
Table$AffectionStatus <- tolower(Table$AffectionStatus)
Table$AffectionStatus <- ifelse(Table$AffectionStatus %in% c("affected", "1", "case"), "Affected",
                                ifelse(Table$AffectionStatus %in% c("control", "0", "healthy"), "Control", NA))
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
#safe_ggsave <- function(plot, filename, width = 36, height = 24) {
safe_ggsave <- function(plot, filename, width = 24, height = 18) {
  tryCatch({
    ggsave(filename, plot = plot, width = width, height = height, dpi = 300)
    message("Saved: ", filename)
  }, error = function(e) message("ggsave failed for ", filename, ": ", e$message))
}

plot_single_histograms <- function(TableLong, metric, outdir = outDir) {
  metric <- toupper(metric)
  metric_long <- TableLong %>% filter(Metric == metric)
  
  if (nrow(metric_long) == 0) {
    message("No data for metric: ", metric)
    return(invisible(NULL))
  }
  
  # x limits if defined
  x_limits <- range_map[[metric]]
  
  # output PDF 1
  fn <- file.path(outdir, paste0(cohort, "_", eName, "_ROIs_histograms_", metric, "_by_AffectionStatus.pdf"))
  pdf(fn, width = 12, height = 8) # one page per ROI
  for (roi in unique(metric_long$parsedROI)) {
    roi_data <- metric_long %>% filter(parsedROI == roi)
    p <- ggplot(roi_data, aes(x = Value, fill = .data[[affectedStatusColumnHeader]])) +
      geom_histogram(aes(y = after_stat(count)), colour = "black", bins = 30) +
      theme_bw() +
      labs(
        x = metric,
        y = "Number of Subjects",
        fill = "AffectionStatus",
        title = paste0(cohort, ": ", metric, " - ", roi)
      ) +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    if (!is.null(x_limits)) p <- p + xlim(x_limits)
    print(p)  # write this page to PDF
  }
  dev.off()
  
  # output PDF 2
  fn <- file.path(outdir, paste0(cohort, "_", eName, "_ROIs_histograms_", metric, ".pdf"))
  pdf(fn, width = 12, height = 8) # one page per ROI
  for (roi in unique(metric_long$parsedROI)) {
    roi_data <- metric_long %>% filter(parsedROI == roi)
    p <- ggplot(roi_data, aes(x = Value)) +
      geom_histogram(aes(y = after_stat(count)), colour = "black", fill = "lightblue", bins = 30) +
      theme_bw() +
      labs(
        x = metric,
        y = "Number of Subjects",
        title = paste0(cohort, ": ", metric, " - ", roi)
      ) +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    if (!is.null(x_limits)) p <- p + xlim(x_limits)
    print(p)  # write this page to PDF
  }
  dev.off()
  
  message("Saved single-ROI histograms for ", metric, " to ", fn)
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
  p <- ggplot(metric_long, aes(x = Value)) +
    geom_histogram(aes(y = after_stat(count)), colour = "black", fill = "lightblue", bins = 30) +
    facet_wrap(~parsedROI, ncol = ncol, scales = "free_x") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    labs(x = metric, y = "Number of Subjects", title = paste0(cohort, ": ", metric))
  if (!is.null(x_limits)) p <- p + xlim(x_limits)
  fn <- file.path(outdir, paste0(cohort, "_", eName, "_ROIs_histograms_multi-panel_", metric, ".pdf"))
  safe_ggsave(p, fn)
  
  p <- ggplot(metric_long, aes(x = Value, fill = .data[[affectedStatusColumnHeader]])) +
    geom_histogram(aes(y = after_stat(count)), colour = "black", bins = 30) +
    facet_wrap(~parsedROI, ncol = ncol, scales = "free_x") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    labs(x = metric, y = "Number of Subjects", fill = "AffectionStatus", title = paste0(cohort, ": ", metric))
  if (!is.null(x_limits)) p <- p + xlim(x_limits)
  fn <- file.path(outdir, paste0(cohort, "_", eName, "_ROIs_histograms_multi-panel_", metric, "_by_AffectionStatus.pdf"))
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
  safe_ggsave(p_all, file.path(outdir, paste0(cohort, "_", eName, "_Age_histogram.pdf")))
}

plot_icv_checks <- function(Table, outdir = outDir) {
  if (!("ICV" %in% names(Table))) {
    message("ICV not present - skipping ICV checks.")
    return(invisible(NULL))
  }
  if (sexColumnHeader %in% names(Table)) {
    tmpTable <- Table[complete.cases(Table[, c(sexColumnHeader, "ICV")]),]  # remove NAs
    # histogram
    p_hist <- ggplot(tmpTable, aes(x = ICV, fill = .data[[sexColumnHeader]])) +
      geom_histogram(aes(y = after_stat(count)), bins = 30, colour = "black") +
      theme_bw() + 
      theme(
        axis.title = element_text(size = 16),   # axis titles (x and y)
        axis.text = element_text(size = 12),     # axis tick labels
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5) # bold & centered
      ) +
      labs(title = paste0(cohort, " - ICV distribution by Sex"), x = expression("Intracranial Volume in mm"^3)) +
      scale_fill_manual(values = c("blue", "red"))  # Optional: Customize the colors
    safe_ggsave(p_hist, file.path(outdir, paste0(cohort, "_", eName, "_ICV_hist_by_Sex.pdf")), width = 12, height = 10)
    # boxplot
    p_box <- ggplot(tmpTable, aes(x = .data[[sexColumnHeader]], y = ICV, fill = .data[[sexColumnHeader]])) +
      geom_boxplot() + 
      theme_bw() + 
      theme(
        axis.title = element_text(size = 16),   # axis titles (x and y)
        axis.text = element_text(size = 12),     # axis tick labels
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5) # bold & centered
      ) +
      labs(title = paste0(cohort, " - ICV by Sex"), y = expression("Intracranial Volume in mm"^3)) +
      scale_fill_manual(values = c("blue", "red"))  # Optional: Customize the colors
    safe_ggsave(p_box, file.path(outdir, paste0(cohort, "_", eName, "_ICV_box_by_Sex.pdf")), width = 12, height = 9)
  }
}

plot_pc <- function(covar, outdir = outDir, pc_cols, phenotype_col = NULL) {
  if (!all(pc_cols %in% colnames(covar))) stop("Some PCs not found in covar")
  
  # --- Define output file ---
  outPDF <- file.path(outdir, paste0(cohort, "_", eName, "_PC_plots.pdf"))
  
  # --- Open PDF device ---
  pdf(outPDF, width = 10, height = 8)
  
  # --- 1. Histograms per PC ---
  for(pc in pc_cols) {
    hist(covar[[pc]], breaks = 30, col = "lightblue", main = paste("Histogram:", pc),
         xlab = pc)
  }
  
  # --- 2. Boxplots of all PCs ---
  boxplot(covar[pc_cols], las = 2, col = "lightgreen", main = "Boxplots of PCs")
  
  # --- 3. Pairwise scatter plots ---
  if(length(pc_cols) >= 2){
    pairs(covar[pc_cols], main="Pairwise PC Scatter Plots", pch=19,
          col=if(!is.null(phenotype_col) && phenotype_col %in% colnames(covar)) {
            as.factor(covar[[phenotype_col]])
          } else "blue")
  }
  
  # --- 4. Outlier detection ---
  outlier_matrix <- apply(covar[pc_cols], 1, function(x) abs(x) > 6)
  outlier_flags <- apply(outlier_matrix, 1, any)
  cat("Number of samples flagged as PC outliers (>6 SD):", sum(outlier_flags), "\n")
  if(sum(outlier_flags) > 0){
    cat("Outlier sample indices:\n")
    print(which(outlier_flags))
  }
  
  # --- 5. Correlation heatmap ---
  library(ggplot2)
  library(reshape2)
  pc_cor <- cor(covar[pc_cols])
  melted_cor <- melt(pc_cor)
  ggplot(melted_cor, aes(Var1, Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
    theme_minimal() +
    labs(title="Correlation Heatmap of PCs", x="", y="") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  dev.off()
  message("PC QC plots saved to: ", outPDF)
}

# ---------------- Summary table: mean ± SD by metric x group ----------------
compute_summary_table <- function(TableLong, outdir = outDir) {
  if (nrow(TableLong) == 0) {
    warning("No long DTI data to summarise.")
    return(data.frame())
  }
  
  # Compute summary by AffectionStatus
  summary_by_status <- TableLong %>%
    group_by(Metric, parsedROI, !!sym(affectedStatusColumnHeader)) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      N = sum(!is.na(Value)),
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compute overall summary (no AffectionStatus grouping)
  summary_all <- TableLong %>%
    group_by(Metric, parsedROI) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      N = sum(!is.na(Value)),
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(!!sym(affectedStatusColumnHeader) := "All")
  
  # Combine both tables
  summary_table <- bind_rows(summary_by_status, summary_all) %>%
    mutate(Mean_SD = sprintf("%.6g \u00B1 %.6g", Mean, SD)) %>%
    rename(AffectionStatus = !!sym(affectedStatusColumnHeader)) %>%
    arrange(Metric, parsedROI, AffectionStatus)
  
  # Save output
  out_csv <- file.path(outdir, paste0(cohort, "_", eName, "_ROIs_SummaryStats.csv"))
  write.csv(summary_table, out_csv, row.names = FALSE)
  message("Saved summary table: ", out_csv)
  
  return(summary_table)
}


# Define a function to process data for a specific group
generate_stats_and_plots <- function(data, group_label, covariate) {
  # Extract column and coerce to numeric
  DATA <- suppressWarnings(as.numeric(as.vector(data[[covariate]])))
  DATA <- DATA[!is.na(DATA)]  # remove NAs
  
  # If empty or all NA, skip gracefully
  if (length(DATA) == 0) {
    warning(paste("No valid numeric data for", covariate, "in group", group_label))
    return(invisible(NULL))
  }
  mu <- mean(DATA)
  sdev <- sd(DATA)
  N <- length(DATA)
  minV <- min(DATA)
  maxV <- max(DATA)
  minSubj <- as.character(data[which(data[[covariate]] == minV)[1], 1])
  maxSubj <- as.character(data[which(data[[covariate]] == maxV)[1], 1])
  minO <- which(DATA < mu - 5 * sdev)
  maxO <- which(DATA > mu + 5 * sdev)
  outliers <- if (length(minO) + length(maxO) > 0) {
    paste(
      "Outliers (5-sd):",
      paste(as.character(data[c(minO, maxO), 1]), collapse = ",")
    )
  } else {
    "None"
  }
  stats <- c(cohort, covariate, group_label, N, mu, sdev, 
             minV, maxV, minSubj, maxSubj, outliers)
  #print(t(as.matrix(stats)))
  
  outPDF <- file.path(outDir, paste0(cohort, "_", eName, "_Age_histograms.pdf"))
  outTXT <- file.path(outDir, paste0(cohort, "_", eName, "_Age_stats.txt"))
  
  write.table(
    t(as.matrix(stats)),
    file = outTXT,
    append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  # Safe histogram (guarded)
  tryCatch({
    hist(
      DATA, breaks = 20,
      main = paste0(cohort, ": ", covariate, " (", group_label, ")"),
      xlab = covariate, col = "lightblue"
    )
  }, error = function(e) {
    warning(paste("Skipping histogram for", covariate, ":", e$message))
  })
  # Safe histogram (guarded)
  tryCatch({
    hist(
      DATA, breaks = 20,
      main = paste0(cohort, ": ", covariate, " (", group_label, ")"),
      xlab = covariate, col = "lightblue", xlim = c(0, 100)
    )
  }, error = function(e) {
    warning(paste("Skipping histogram with xlim = c(0, 100) for", covariate, ":", e$message))
  })
  #hist(DATA, breaks = 20, main = paste0(cohort, ": ", "Age in Years"), xlab = "Age in Years", col = "lightgray")
  #hist(DATA, breaks = 20, main = paste0(cohort, ": ", "Age in Years"), xlab = "Age in Years", col = "lightgray", xlim = c(0, 100))
  #hist(DATA, breaks = 20, main = paste0(cohort, ": ", "Age in Years"), xlab = "Age in Years", col = "lightblue")
  #hist(DATA, breaks = 20, main = paste0(cohort, ": ", "Age in Years"), xlab = "Age in Years", col = "lightblue", xlim = c(0, 100))
}

# Function to generate group size plots
generate_group_size_plots <- function(data, group_var, group_label) {
  group_counts <- table(data[[group_var]])
  barplot(
    group_counts,
    main = paste0(
      cohort,
      ": Group sizes by ",
      group_var,
      " (",
      group_label,
      ")"
    ),
    col = "lightgreen",
    xlab = group_var,
    ylab = "Count (N)"
  )
}

# Function to create overlapping histograms
generate_overlapping_histograms <- function(data, group_var, covariate, group_label) {
  groups <- unique(data[[group_var]])
  
  # Generate a palette excluding red
  color_palette <- c(
    "blue",
    "green",
    "orange",
    "purple",
    "cyan",
    "yellow",
    "pink",
    "brown"
  )
  colors <- color_palette[seq_along(groups)]
  
  hist_list <- list()
  valid_groups <- character()
  
  # Build histograms for each group
  for (group in groups) {
    group_data <- data[data[[group_var]] == group, covariate]
    
    # Convert to numeric safely
    group_data <- suppressWarnings(as.numeric(as.vector(group_data)))
    
    # Skip if empty or all NA
    if (length(group_data) == 0 || all(is.na(group_data))) {
      warning(paste0("Skipping group '", group, "' for covariate '", covariate, "' — no valid data."))
      next
    }
    
    # Only include non-empty valid groups
    hist_list[[length(hist_list) + 1]] <- hist(
      group_data,
      breaks = 20,
      plot = FALSE
    )
    valid_groups <- c(valid_groups, group) # Track valid groups
  }
  
  # Check if there are valid histograms
  if (length(hist_list) > 0) {
    # Set up the base plot using the first valid histogram
    plot(
      hist_list[[1]],
      col = adjustcolor(colors[1], alpha.f = 0.5),
      main = paste0(cohort, ": ", covariate, " by ", group_var, " (", group_label, ")"),
      xlab = covariate,
      xlim = range(unlist(lapply(hist_list, function(h) h$breaks))),
      ylim = range(unlist(lapply(hist_list, function(h) h$counts))),
      border = colors[1]
    ) # Set border color to match bar color
    #border = NA)  # No border for transparency
    
    # Add additional histograms
    if (length(hist_list) > 1) {
      for (i in 2:length(hist_list)) {
        plot(
          hist_list[[i]],
          col = adjustcolor(colors[i], alpha.f = 0.5),
          add = TRUE,
          border = colors[i]
        ) # Set border color to match bar color
        #border = NA) # No border for transparency
      }
    }
    
    # Add a legend for the valid groups
    legend(
      "topright",
      legend = valid_groups,
      fill = adjustcolor(colors[1:length(hist_list)], alpha.f = 0.5)
    )
  } else {
    warning(paste("No valid data for overlapping histograms by", group_var, "in", group_label))
  }
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
  out_csv <- file.path(outdir, paste0(cohort, "_", eName, "_Cohensd_results.csv"))
  write.csv(results, out_csv, row.names = FALSE)
  message("Saved Cohen's d results: ", out_csv)
  return(results)
}

# ---------------- Run the steps ----------------
message("Running QC pipeline...")

# Per-metric multi-panel histograms
for (m in c("FA","MD","AD","RD")) {
  plot_multi_histogram(Table, m, outdir = outDir)
  plot_single_histograms(TableLong, metric = m, outdir = outDir)
}

# Build & save summary table from long data
summary_table <- compute_summary_table(TableLong, outDir)

# Age histograms and group breakdowns
plot_age_histograms(Table, outDir)

# ICV checks
plot_icv_checks(Table, outDir)

# PC checks
# optional: phenotype column for coloring scatter plots
plot_pc(Table, outDir, pc_cols = pc_columns, phenotype_col = "CaseControl")

# Open PDF for plots
pdf(file = paste0(outDir, "/", cohort, "_", eName, "_Age_histograms.pdf"))
write(
  "Cohort\tCovariate\tGroup\tNumberIncluded\tMean\tStandDev\tMinValue\tMaxValue\tMinSubject\tMaxSubject\t5StDev_Off",
  file = paste0(outDir, "/", cohort, "_", eName, "_Age_stats.txt")
)

# Process Age for the entire cohort
if ("Age" %in% colnames(Table)) {
  generate_stats_and_plots(Table, "All", "Age")
  
  # Split by Sex if the variable exists
  if ("Sex" %in% colnames(Table)) {
    for (sex in unique(Table$Sex)) {
      subset_sex <- Table[Table$Sex == sex, ]
      generate_stats_and_plots(subset_sex, paste0("Sex: ", sex), "Age")
    }
    generate_group_size_plots(Table, "Sex", "All")
    generate_overlapping_histograms(Table, "Sex", "Age", "All")
  }
  
  # Split by AffectionStatus if the variable exists
  if ("AffectionStatus" %in% colnames(Table)) {
    for (status in unique(Table$AffectionStatus)) {
      subset_status <- Table[Table$AffectionStatus == status, ]
      generate_stats_and_plots(subset_status, paste0("AffectionStatus: ", status), "Age")
    }
    generate_group_size_plots(Table, "AffectionStatus", "All")
    generate_overlapping_histograms(Table, "AffectionStatus", "Age", "All")
  }
  
  # Split by both Sex and AffectionStatus if both variables exist
  if ("Sex" %in% colnames(Table) && "AffectionStatus" %in% colnames(Table)) {
    for (sex in unique(Table$Sex)) {
      for (status in unique(Table$AffectionStatus)) {
        subset_combined <- Table[
          Table$Sex == sex & Table$AffectionStatus == status,
        ]
        if (nrow(subset_combined) > 0) {
          generate_stats_and_plots(subset_combined, paste0("Sex: ", sex, ", AffectionStatus: ", status), "Age")
        }
      }
    }
    generate_group_size_plots(Table, "Sex", "Grouped by Sex and AffectionStatus")
    generate_group_size_plots(Table,"AffectionStatus","Grouped by Sex and AffectionStatus")
    generate_overlapping_histograms(Table,"Sex","Age","Grouped by Sex and AffectionStatus")
    generate_overlapping_histograms(Table,"AffectionStatus","Age","Grouped by Sex and AffectionStatus")
  }
} else {
  cat("Covariate not found in the table:", "Age", "\n")
}

# Close PDF device
dev.off()

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
  safe_ggsave(p_cd, file.path(outDir, paste0(cohort, "_", eName, "_Cohensd_bar_graph_AFF-CON.png")), width = 18, height = 24)
} else {
  message("No Cohen's d results to plot.")
}

message("QC pipeline finished. Outputs in: ", normalizePath(outDir))
