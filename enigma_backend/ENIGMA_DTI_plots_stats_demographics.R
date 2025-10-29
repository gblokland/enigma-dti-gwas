####
####  ENIGMA DTI ###
####
### This is a function to print out plots/stats for Quality Control of ENIGMA-DTI Covariates
#############################################################################################
### Written by Gabriella Blokland, based on script from Neda Jahanshad / Derrek Hibar
### Last update October 2025
### Questions or Comments:
### enigma.dtigenetics@gmail.com
#############################################################################################

# Check and install argparse if not already installed
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "QC script for ENIGMA DTI")
parser$add_argument(
  "--cohort",
  required = TRUE,
  help = "Cohort name"
)
parser$add_argument(
  "--covarFILE",
  default = "${COHORT}_enigma_dti_gwas.covar",
  required = TRUE,
  help = "Input text file containing covariates"
)
parser$add_argument(
  "--phenoFILE",
  default = "${COHORT}_enigma_dti_gwas.pheno",
  required = TRUE,
  help = "Input text file containing phenotypes (MRI/DTI)"
)
parser$add_argument(
  "--ageColumnHeader",
  default = "Age",
  help = "Name of the Age variable"
)
parser$add_argument(
  "--sexColumnHeader",
  default = "Sex",
  help = "Name of the Sex variable"
)
parser$add_argument(
  "--maleIndicator",
  default = "1",
  help = "value used for male individuals"
)
parser$add_argument(
  "--CaseControlCohort",
  default = "1",
  help = "CaseControlCohort yes(1) or no(0)"
)
parser$add_argument(
  "--affectedStatusColumnHeader",
  default = "AffectionStatus",
  help = "Column header for Affection Status (case-control status)"
)
parser$add_argument(
  "--affectedIndicator",
  default = "2",
  help = "value used for affected individuals"
)
parser$add_argument(
  "--related",
  default = "0",
  help = "related cohort yes(1) or no(0)"
)
parser$add_argument(
  "--rois",
  default = "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R",
  help = "Semicolon-separated list of ROIs"
)
parser$add_argument(
  "--pheno_covar_dir",
  default = "./QC_ENIGMA/",
  help = "Output directory"
)
parser$add_argument(
  "--outDir",
  default = "./QC_ENIGMA/",
  help = "Output directory"
)
parser$add_argument(
  "--outPDF",
  default = "ENIGMA_DTI_Age_histograms.pdf",
  help = "Output PDF for histograms"
)
parser$add_argument(
  "--outTXT",
  default = "ENIGMA_DTI_Age_stats.txt",
  help = "Output TXT for statistics"
)
parser$add_argument(
  "--eName",
  default = "ENIGMA_DTI_GWAS",
  help = "Output label for ENIGMA project"
)

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

# Create output directory if it doesn't exist
dir.create(outDir, showWarnings = FALSE)

# Read input CSV file
Table <- read.csv(CSVfile, header = TRUE)
colTable <- names(Table)
print(colTable)

# Replace "x" or "X" values with NA
for (m in seq_along(colTable)) {
  ind <- which(Table[, m] == "x")
  ind2 <- which(Table[, m] == "X")
  Table[ind, m] <- NA
  Table[ind2, m] <- NA
}

# Remove rows with NA values
Table <- Table[complete.cases(Table), ]

# Open PDF for plots
pdf(file = paste0(outDir, "/", outPDF))
write(
  "Cohort\tCovariate\tGroup\tNumberIncluded\tMean\tStandDev\tMinValue\tMaxValue\tMinSubject\tMaxSubject\t5StDev_Off",
  file = paste0(outDir, "/", outTXT)
)

# Define a function to process data for a specific group
generate_stats_and_plots <- function(data, group_label, covariate) {
  DATA <- as.numeric(as.vector(data[[covariate]]))
  mu <- mean(DATA)
  sdev <- sd(DATA)
  N <- length(DATA)
  minV <- min(DATA)
  maxV <- max(DATA)
  minSubj <- as.character(data[which(DATA == minV), 1])
  maxSubj <- as.character(data[which(DATA == maxV), 1])
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
  stats <- c(
    cohort,
    covariate,
    group_label,
    N,
    mu,
    sdev,
    minV,
    maxV,
    minSubj,
    maxSubj,
    outliers
  )
  write.table(
    t(as.matrix(stats)),
    file = paste0(outDir, "/", outTXT),
    append = TRUE,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  hist(
    DATA,
    breaks = 20,
    main = paste0(cohort, ": ", covariate, " (", group_label, ")"),
    xlab = covariate,
    col = "lightblue"
  )
  hist(
    DATA,
    breaks = 20,
    main = paste0(cohort, ": ", covariate, " (", group_label, ")"),
    xlab = covariate,
    col = "lightblue",
    xlim = c(0, 100)
  )
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
generate_overlapping_histograms <- function(
  data,
  group_var,
  covariate,
  group_label
) {
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
    group_data <- as.numeric(as.vector(data[
      data[[group_var]] == group,
      covariate
    ]))
    if (length(group_data) > 0) {
      # Only include non-empty groups
      hist_list[[length(hist_list) + 1]] <- hist(
        group_data,
        breaks = 20,
        plot = FALSE
      )
      valid_groups <- c(valid_groups, group) # Track valid groups
    }
  }

  # Check if there are valid histograms
  if (length(hist_list) > 0) {
    # Set up the base plot using the first valid histogram
    plot(
      hist_list[[1]],
      col = adjustcolor(colors[1], alpha.f = 0.5),
      main = paste0(cohort, ": Overlapping Histograms of ", covariate, " by ", group_var, " (", group_label, ")"),
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
      generate_stats_and_plots(
        subset_status,
        paste0("AffectionStatus: ", status),
        "Age"
      )
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
          generate_stats_and_plots(
            subset_combined,
            paste0("Sex: ", sex, ", AffectionStatus: ", status),
            "Age"
          )
        }
      }
    }
    generate_group_size_plots(
      Table,
      "Sex",
      "Grouped by Sex and AffectionStatus"
    )
    generate_group_size_plots(
      Table,
      "AffectionStatus",
      "Grouped by Sex and AffectionStatus"
    )
    generate_overlapping_histograms(
      Table,
      "Sex",
      "Age",
      "Grouped by Sex and AffectionStatus"
    )
    generate_overlapping_histograms(
      Table,
      "AffectionStatus",
      "Age",
      "Grouped by Sex and AffectionStatus"
    )
  }
} else {
  cat("Covariate not found in the table:", "Age", "\n")
}

# Close PDF device
dev.off()
