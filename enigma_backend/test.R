#### ENIGMA DTI ###
####
### This is a function to print out plots/stats for Quality Control of ENIGMA-DTI Covariates
#############################################################################################
### Author: Gabriella Blokland, based on script from Neda Jahanshad / Derrek Hibar
### Last update October 2025 (fixed and hardened Nov 2025)
### Questions or Comments:
### enigma.dtigenetics@gmail.com
#############################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check and install argparse if not already installed
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "QC script for ENIGMA DTI")
parser$add_argument("--cohort",  required = TRUE,  help = "Cohort name")
parser$add_argument("--covarFILE",  default = "${COHORT}_enigma_dti_gwas.covar",  required = TRUE,  help = "Input text file containing covariates")
parser$add_argument("--phenoFILE",  default = "${COHORT}_enigma_dti_gwas.pheno",  required = TRUE,  help = "Input text file containing phenotypes (MRI/DTI)")
parser$add_argument("--ageColumnHeader",  default = "Age",  help = "Name of the Age variable")
parser$add_argument("--sexColumnHeader",  default = "Sex",  help = "Name of the Sex variable")
parser$add_argument("--maleIndicator",  default = "1",  help = "value used for male individuals")
parser$add_argument("--CaseControlCohort",  default = "1",  help = "CaseControlCohort yes(1) or no(0)")
parser$add_argument("--affectedStatusColumnHeader",  default = "AffectionStatus",  help = "Column header for Affection Status (case-control status)")
parser$add_argument("--affectedIndicator",  default = "2",  help = "value used for affected individuals")
parser$add_argument("--related",  default = "0",  help = "related cohort yes(1) or no(0)")
parser$add_argument("--rois",  default = "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R",  help = "Semicolon-separated list of ROIs")
parser$add_argument("--pheno_covar_dir",  default = "./QC_ENIGMA/",  help = "Output directory")
parser$add_argument("--outDir",  default = "./QC_ENIGMA/",  help = "Output directory")
parser$add_argument("--outPDF",  default = "ENIGMA_DTI_Age_histograms.pdf",  help = "Output PDF for histograms")
parser$add_argument("--outTXT",  default = "ENIGMA_DTI_Age_stats.txt",  help = "Output TXT for statistics")
parser$add_argument("--eName",  default = "ENIGMA_DTI_GWAS",  help = "Output label for ENIGMA project")

args <- parser$parse_args()

# Extract arguments
cohort <- args$cohort
covarFILE <- args$covarFILE
phenoFILE <- args$phenoFILE
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

# Replace ${COHORT} in filenames (if user left the template string)
covarFILE <- gsub("\\$\\{COHORT\\}", cohort, covarFILE)
phenoFILE <- gsub("\\$\\{COHORT\\}", cohort, phenoFILE)

# Create output directory if it doesn't exist
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# Read input file (stringsAsFactors = FALSE to avoid factor surprises)
if (!file.exists(covarFILE)) {
  stop("Covariate file not found: ", covarFILE)
}
Table <- read.table(covarFILE, header = TRUE, stringsAsFactors = FALSE, sep = "", comment.char = "", quote = "\"'", strip.white = TRUE)
colTable <- names(Table)
message("Columns in covariate file: ", paste(colTable, collapse = ", "))

# Replace literal "x" or "X" with NA (works across column types)
Table[Table == "x" | Table == "X"] <- NA

# Remove rows with NA values in ALL columns? original used complete.cases(Table)
# Keep behaviour: remove rows with ANY NA
Table <- Table[complete.cases(Table), ]

if (nrow(Table) == 0) {
  stop("No complete cases remain after removing NAs.")
}

# Prepare output TXT header
header_line <- paste0("Cohort\tCovariate\tGroup\tNumberIncluded\tMean\tStandDev\tMinValue\tMaxValue\tMinSubject\tMaxSubject\t5StDev_Off")
writeLines(header_line, con = file.path(outDir, outTXT))

# Open PDF for plots
pdf(file = file.path(outDir, outPDF))

# Helper: safe numeric conversion
to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# Define a function to process data for a specific group
generate_stats_and_plots <- function(data, group_label, covariate) {
  if (!(covariate %in% colnames(data))) {
    warning("Covariate not present: ", covariate)
    return(invisible(NULL))
  }

  DATA_raw <- data[[covariate]]
  DATA <- to_num(DATA_raw)
  # drop NAs
  valid_idx <- which(!is.na(DATA))
  DATA <- DATA[valid_idx]

  if (length(DATA) == 0) {
    warning("No numeric data for covariate ", covariate, " in group ", group_label)
    return(invisible(NULL))
  }

  mu <- mean(DATA)
  sdev <- sd(DATA)
  N <- length(DATA)
  minV <- min(DATA)
  maxV <- max(DATA)

  # find subject IDs (first column assumed to be FID/IID or similar)
  subj_col <- 1
  # map to original data rows
  # locate subject where covariate equals min/max (use first match)
  full_vals <- to_num(as.character(data[[covariate]]))
  min_idx_all <- which(!is.na(full_vals) & full_vals == minV)
  max_idx_all <- which(!is.na(full_vals) & full_vals == maxV)
  minSubj <- if (length(min_idx_all) > 0) as.character(data[min_idx_all[1], subj_col]) else NA
  maxSubj <- if (length(max_idx_all) > 0) as.character(data[max_idx_all[1], subj_col]) else NA

  minO <- which(DATA < mu - 5 * sdev)
  maxO <- which(DATA > mu + 5 * sdev)
  outlier_subjects <- character()
  if (length(minO) + length(maxO) > 0) {
    # translate indices back to subject ids
    # valid_idx holds indices in original `data`
    all_out_idx <- valid_idx[c(minO, maxO)]
    outlier_subjects <- as.character(data[all_out_idx, subj_col])
    outliers <- paste("Outliers (5-sd):", paste(outlier_subjects, collapse = ","))
  } else {
    outliers <- "None"
  }

  stats <- c(
    cohort,
    covariate,
    group_label,
    N,
    sprintf("%.4f", mu),
    sprintf("%.4f", sdev),
    sprintf("%.4f", minV),
    sprintf("%.4f", maxV),
    minSubj,
    maxSubj,
    outliers
  )

  write.table(
    t(as.matrix(stats)),
    file = file.path(outDir, outTXT),
    append = TRUE,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
  )

  # Plot histograms
  hist(DATA,
       breaks = 20,
       main = paste0(cohort, ": ", covariate, " (", group_label, ")"),
       xlab = covariate,
       col = "lightblue"
  )

  # second histogram with fixed xlim 0-100 if meaningful (ignore if data outside)
  xlim_fixed <- c(0, 100)
  if (min(DATA, na.rm = TRUE) < xlim_fixed[2] && max(DATA, na.rm = TRUE) > xlim_fixed[1]) {
    hist(DATA,
         breaks = 20,
         main = paste0(cohort, ": ", covariate, " (", group_label, ") limited 0-100"),
         xlab = covariate,
         col = "lightblue",
         xlim = xlim_fixed
    )
  }
}

# Function to generate group size plots
generate_group_size_plots <- function(data, group_var, group_label) {
  if (!(group_var %in% colnames(data))) {
    warning("Group var not in table: ", group_var)
    return(invisible(NULL))
  }
  group_counts <- table(data[[group_var]])
  barplot(
    group_counts,
    main = paste0(cohort, ": Group sizes by ", group_var, " (", group_label, ")"),
    xlab = group_var,
    ylab = "Count (N)",
    col = "lightgreen"
  )
}

# Function to create overlapping histograms
generate_overlapping_histograms <- function(data, group_var, covariate, group_label) {
  if (!(group_var %in% colnames(data))) {
    warning("Group var not in table: ", group_var)
    return(invisible(NULL))
  }
  if (!(covariate %in% colnames(data))) {
    warning("Covariate not in table: ", covariate)
    return(invisible(NULL))
  }

  groups <- unique(data[[group_var]])
  groups <- groups[!is.na(groups)]
  if (length(groups) == 0) {
    warning("No groups present for ", group_var)
    return(invisible(NULL))
  }

  # palette excluding red
  color_palette <- c("blue", "green", "orange", "purple", "cyan", "yellow", "pink", "brown")
  colors <- rep(color_palette, length.out = length(groups))

  hist_list <- list()
  valid_groups <- character()

  for (i in seq_along(groups)) {
    group <- groups[i]
    group_rows <- data[data[[group_var]] == group, , drop = FALSE]
    group_data <- to_num(group_rows[[covariate]])
    group_data <- group_data[!is.na(group_data)]
    if (length(group_data) > 0) {
      hist_list[[length(hist_list) + 1]] <- hist(group_data, breaks = 20, plot = FALSE)
      valid_groups <- c(valid_groups, as.character(groups[i]))
    }
  }

  if (length(hist_list) == 0) {
    warning("No valid numeric data for overlapping histograms by ", group_var)
    return(invisible(NULL))
  }

  # set up a plotting range
  xlim_range <- range(unlist(lapply(hist_list, function(h) h$breaks)), na.rm = TRUE)
  ylim_range <- range(unlist(lapply(hist_list, function(h) h$counts)), na.rm = TRUE)

  # plot first histogram (as base)
  plot(hist_list[[1]],
       main = paste0(cohort, ": Histograms of ", covariate, " by ", group_var, " (", group_label, ")"),
       xlab = covariate,
       xlim = xlim_range,
       ylim = ylim_range,
       col = adjustcolor(colors[1], alpha.f = 0.5),
       border = colors[1]
  )

  if (length(hist_list) > 1) {
    for (i in 2:length(hist_list)) {
      plot(hist_list[[i]],
           add = TRUE,
           col = adjustcolor(colors[i], alpha.f = 0.5),
           border = colors[i]
      )
    }
  }

  legend("topright", legend = valid_groups, fill = adjustcolor(colors[1:length(hist_list)], alpha.f = 0.5), cex = 0.8)
}

### ---- Now process Age using the configured header names ----

age_col <- ageColumnHeader
sex_col <- sexColumnHeader
aff_col <- affectedStatusColumnHeader

if (age_col %in% colnames(Table)) {
  generate_stats_and_plots(Table, "All", age_col)

  # Split by Sex if the variable exists
  if (sex_col %in% colnames(Table)) {
    for (sex in unique(Table[[sex_col]])) {
      subset_sex <- Table[Table[[sex_col]] == sex, , drop = FALSE]
      if (nrow(subset_sex) > 0) {
        generate_stats_and_plots(subset_sex, paste0(sex_col, ": ", sex), age_col)
      }
    }
    generate_group_size_plots(Table, sex_col, "All")
    generate_overlapping_histograms(Table, sex_col, age_col, "All")
  }

  # Split by AffectionStatus if the variable exists
  if (aff_col %in% colnames(Table)) {
    for (status in unique(Table[[aff_col]])) {
      subset_status <- Table[Table[[aff_col]] == status, , drop = FALSE]
      if (nrow(subset_status) > 0) {
        generate_stats_and_plots(subset_status, paste0(aff_col, ": ", status), age_col)
      }
    }
    generate_group_size_plots(Table, aff_col, "All")
    generate_overlapping_histograms(Table, aff_col, age_col, "All")
  }

  # Split by both Sex and AffectionStatus if both variables exist
  if (sex_col %in% colnames(Table) && aff_col %in% colnames(Table)) {
    for (sex in unique(Table[[sex_col]])) {
      for (status in unique(Table[[aff_col]])) {
        subset_combined <- Table[Table[[sex_col]] == sex & Table[[aff_col]] == status, , drop = FALSE]
        if (nrow(subset_combined) > 0) {
          generate_stats_and_plots(subset_combined,
                                   paste0(sex_col, ": ", sex, ", ", aff_col, ": ", status),
                                   age_col)
        }
      }
    }
    generate_group_size_plots(Table, sex_col, "Grouped by Sex and AffectionStatus")
    generate_group_size_plots(Table, aff_col, "Grouped by Sex and AffectionStatus")
    generate_overlapping_histograms(Table, sex_col, age_col, "Grouped by Sex and AffectionStatus")
    generate_overlapping_histograms(Table, aff_col, age_col, "Grouped by Sex and AffectionStatus")
  }

} else {
  cat("Covariate not found in the table:", age_col, "\n")
}

# Close PDF device
dev.off()

message("QC script finished. Outputs written to: ", outDir)
