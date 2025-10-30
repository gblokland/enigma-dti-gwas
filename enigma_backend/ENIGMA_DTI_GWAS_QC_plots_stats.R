#!/usr/bin/env Rscript
# ===============================================================
# ENIGMA-DTI QC Combined Script
#   Performs data loading, basic QC, histograms, and Cohen's d
#   Author: (Your name)
#   Date: (auto)
# ===============================================================

# ------------------ SETUP ------------------
suppressPackageStartupMessages({
  if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("MBESS", quietly = TRUE)) install.packages("MBESS")
})

library(argparse)
library(ggplot2)
library(dplyr)
library(MBESS)

# ------------------ ARGUMENT PARSER ------------------
parser <- ArgumentParser(description = "ENIGMA-DTI-GWAS QC pipeline")
parser$add_argument("--cohort", required = TRUE, help = "Cohort name")
parser$add_argument("--covarFILE", required = TRUE, help = "Covariate file")
parser$add_argument("--phenoFILE", required = TRUE, help = "Phenotype file")
parser$add_argument("--outDir", default = ".", help = "Output directory")
args <- parser$parse_args()

cohort <- args$cohort
covarFILE <- args$covarFILE
phenoFILE <- args$phenoFILE
outDir <- args$outDir

# ------------------ LOAD DATA ------------------
message("Loading data...")
covar <- read.table(covarFILE, header = TRUE)
pheno <- read.table(phenoFILE, header = TRUE)
Table <- merge(covar, pheno, by = c("FID", "IID"))
message("Data loaded: ", nrow(Table), " subjects")

# ------------------ BASIC RANGES ------------------
FA_range <- c(0, 1)
MD_range <- c(0, 0.002)
AD_range <- c(0, 0.0025)
RD_range <- c(0, 0.002)

# ------------------ PARSED VARIABLES ------------------
parsedROIs <- unique(gsub(".*_", "", grep("FA_|MD_|AD_|RD_", names(Table), value = TRUE)))
roiLabels <- parsedROIs
roiColors <- rep("steelblue", length(parsedROIs))

# ===============================================================
# 1. HISTOGRAM FUNCTIONS
# ===============================================================
plot_histogram <- function(data, variable, metric = variable) {
  p <- ggplot(data, aes(x = .data[[variable]], fill = AffectionStatus)) +
    geom_histogram(aes(y = after_stat(count)), colour = "black", bins = 30) +
    xlab(variable) + ylab("Number of Subjects") +
    theme_bw() +
    theme(legend.position = "top")
  plot_fd <- file.path(outDir, paste0(cohort, "_histogram_", variable, ".pdf"))
  ggsave(plot_fd, plot = p, width = 20, height = 15, units = "cm")
}

plot_multi_histogram <- function(Table, metric) {
  metric_cols <- grep(paste0("^", metric, "_"), names(Table), value = TRUE)
  TableLongMetric <- reshape(
    Table,
    varying = metric_cols,
    v.names = "Value",
    timevar = "parsedROIs",
    times = gsub(paste0(metric, "_"), "", metric_cols),
    idvar = c("FID", "IID", "AffectionStatus"),
    direction = "long"
  )
  p <- ggplot(TableLongMetric, aes(x = Value, fill = AffectionStatus)) +
    geom_histogram(aes(y = after_stat(count)), bins = 30, colour = "black", size = 0.2) +
    facet_wrap(~parsedROIs, ncol = 8, scales = "free_x") +
    xlab(metric) + ylab("Number of Subjects") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 8, face = "bold"))
  plot_fd <- file.path(outDir, paste0(cohort, "_multi_hist_", metric, ".pdf"))
  ggsave(plot_fd, plot = p, height = 24, width = 36, units = "cm", dpi = 600)
}

# ===============================================================
# 2. COHEN'S D CALCULATION
# ===============================================================
compute_cohens_d <- function(Table, parsedROIs, roiLabels, roiColors, cohort, outDir) {
  cohens_d_results <- data.frame()
  n_CON <- nrow(Table[Table$AffectionStatus == "Control", ])
  n_AFF <- nrow(Table[Table$AffectionStatus == "Affected", ])

  if (n_CON > 5 & n_AFF > 5) {
    for (metric in c("FA", "MD", "AD", "RD")) {
      metric_cols <- grep(paste0("^", metric, "_"), names(Table), value = TRUE)
      for (i in seq_along(metric_cols)) {
        group_aff <- Table[Table$AffectionStatus == "Affected", metric_cols[i]]
        group_con <- Table[Table$AffectionStatus == "Control", metric_cols[i]]

        if (all(is.na(group_aff)) || all(is.na(group_con))) next

        pooled_sd <- sqrt(
          (sum((group_con - mean(group_con, na.rm = TRUE))^2, na.rm = TRUE) +
            sum((group_aff - mean(group_aff, na.rm = TRUE))^2, na.rm = TRUE)) /
            (length(na.omit(group_con)) + length(na.omit(group_aff)) - 2)
        )

        TT <- t.test(group_aff, group_con)
        N1 <- length(na.omit(group_aff))
        N2 <- length(na.omit(group_con))
        ci <- MBESS::ci.smd(ncp = TT$statistic, n.1 = N1, n.2 = N2, conf.level = .95)

        cohens_d_results <- rbind(cohens_d_results, data.frame(
          metric = metric,
          variable = metric_cols[i],
          varlabel = gsub(paste0(metric, "_"), "", metric_cols[i]),
          colour = roiColors[i],
          cohens_d_smd = ci$smd,
          cohens_d_cilb = ci$Lower.Conf.Limit.smd,
          cohens_d_ciub = ci$Upper.Conf.Limit.smd,
          cohens_d_se = (ci$Upper.Conf.Limit.smd - ci$Lower.Conf.Limit.smd) / 3.92
        ))
      }
    }

    # Save results
    write.table(cohens_d_results, file.path(outDir, paste0(cohort, "_cohens_d_results.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Plot by metric
    p <- ggplot(cohens_d_results, aes(x = varlabel, y = cohens_d_smd, fill = metric)) +
      geom_bar(stat = "identity", colour = "black") +
      geom_errorbar(aes(ymin = cohens_d_smd - cohens_d_se,
                        ymax = cohens_d_smd + cohens_d_se),
                    width = 0.25) +
      coord_flip() +
      facet_wrap(~metric, ncol = 2, scales = "free_y") +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
      geom_hline(yintercept = 0) +
      xlab("") + ylab("Cohen's d ± SE") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8, face = "bold"))
    ggsave(file.path(outDir, paste0(cohort, "_cohensd_bar_graph_AFF-CON.png")),
           plot = p, height = 24, width = 18, units = "cm", dpi = 600)
  } else {
    message("Not enough subjects per group for Cohen's d.")
  }

  return(cohens_d_results)
}

# ===============================================================
# 3. RUN ALL STEPS
# ===============================================================
message("Running QC plots and Cohen’s d...")

# Example QC plots
if ("ICV" %in% names(Table)) {
  plot_histogram(Table, "ICV")
  plot_multi_histogram(Table, "FA")
}

# Compute effect sizes
cohens_d_results <- compute_cohens_d(Table, parsedROIs, roiLabels, roiColors, cohort, outDir)

message("All analyses complete.")
# ===============================================================
