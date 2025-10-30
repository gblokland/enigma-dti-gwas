#
#  ENIGMA DTI GWAS
#
# This is a function to print out plots and stats for Quality Control of ENIGMA-DTI measures
#############################################################################################
# Written by Gabriella Blokland, based on script from Neda Jahanshad / Derrek Hibar
# Last update September 2025
# Questions or Comments:
# enigma.dtigenetics@gmail.com
#############################################################################################

## to run:
## ${Rbin} --no-save --slave --args ${1} ${2} ...   ${8} <  ./.R
## R --no-save --slave --args ${covarFILE} ${phenoFILE} ${ageColumnHeader} ${sexColumnHeader} ${maleIndicator} \
## ${CaseControlCohort} ${affectedStatusColumnHeader} ${affectedIndicator} ${related} ${pheno_covar_dir} ${ALL_ROIS} ${eName} <  ${run_directory}/ENIGMA_DTI_phenotype_and_covariate_sumstats.R

options(repos = c(CRAN = "https://cloud.r-project.org/"))

#Check and install argparse if not already installed
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
suppressMessages(suppressWarnings(library(argparse)))

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
  default = "ENIGMA_DTI_allROI_histograms.pdf",
  help = "Output PDF for histograms"
)
parser$add_argument(
  "--outTXT",
  default = "ENIGMA_DTI_allROI_stats.txt",
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

# Read covariates CSV file
covarTable <- read.table(covarFILE, header = TRUE)
colCovarTable <- names(covarTable)
print(colCovarTable)

# Create output directory if it doesn't exist
dir.create(outDir, showWarnings = FALSE)

# Read phenotypes CSV file
phenoTable <- read.table(phenoFILE, header = TRUE)
colTable <- names(phenoTable)
print(colTable)

# Replace "x" or "X" values with NA
for (m in seq_along(colTable)) {
  ind <- which(phenoTable[, m] == "x")
  ind2 <- which(phenoTable[, m] == "X")
  phenoTable[ind, m] <- NA
  phenoTable[ind2, m] <- NA
}

# Merge covar and pheno tables
Table <- merge(covarTable, phenoTable, by = c("FID", "IID"), all = TRUE)

# Remove rows with NA values
Table <- Table[complete.cases(Table), ]

# Parse ROIs
parsedROIs <- unlist(strsplit(rois, ";"))
Nrois <- length(parsedROIs)
writeLines(paste0("Nrois = ", Nrois))

# Function to calculate Cohen's d
calculate_cohens_d <- function(group1, group2) {
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  sd1 <- sd(group1)
  sd2 <- sd(group2)
  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  d <- (mean1 - mean2) / pooled_sd
  return(d)
}

if (Nrois > 0) {
  dev.new()
  pdf(file = paste0(outDir, "/", outPDF))
  write(
    "Phenotype\tNumberIncluded\tMean\tStandDev\tMinValue\tMaxValue\tMinSubject\tMaxSubject\t5StDev_Off\tCohensD",
    file = paste0(outDir, "/", outTXT)
  )

  for (x in seq_len(Nrois)) {
    for (metric in c("FA", "MD", "AD", "RD")) {
      #ROI <- parsedROIs[x]
      ROI <- paste0(metric, "_", parsedROIs[x])

      if (ROI %in% colnames(Table)) {
        writeLines(paste0("Current ROI = ", ROI))

        DATA <- as.numeric(as.vector(Table[, ROI]))
        print(head(DATA))

        mu <- mean(DATA)
        writeLines(paste0("Mean = ", mu))
        sdev <- sd(DATA)
        writeLines(paste0("SD = ", sdev))
        N <- length(DATA)
        writeLines(paste0("N = ", N))

        minV <- min(DATA)
        writeLines(paste0("Min = ", minV))
        maxV <- max(DATA)
        writeLines(paste0("Max = ", maxV))

        i <- which(DATA == minV)
        minSubj <- as.character(Table[i, 1])
        writeLines(paste0("minSubj = ", minSubj))

        j <- which(DATA == maxV)
        maxSubj <- as.character(Table[j, 1])
        writeLines(paste0("maxSubj = ", maxSubj))

        minO <- which(DATA < mu - 5 * sdev)
        writeLines(paste0("minO = ", minO))
        maxO <- which(DATA > mu + 5 * sdev)
        writeLines(paste0("maxO = ", maxO))

        outliers <- if (length(minO) + length(maxO) > 0) {
          paste("Outliers (5-sd):", paste(as.character(Table[c(minO, maxO), 1]), collapse = ","))
        } else {
          "None"
        }

        # Check if both groups exist in the data
        if (all(c("1", "2") %in% unique(Table$AffectionStatus))) {
          ### Case-Control group
          #group_var <- Table$group  # Assuming 'group' column exists and has values 'case' and 'control'
          #group_case <- DATA[group_var == "case"]
          #group_control <- DATA[group_var == "control"]
          group_var <- Table$AffectionStatus # Assuming 'group' column exists and has values 'case' and 'control'
          group_case <- DATA[group_var == "2"]
          group_control <- DATA[group_var == "1"]
          ### Calculate Cohen's d for case-control comparison
          cohens_d <- calculate_cohens_d(group_case, group_control)
          writeLines(paste0("Cohen's d = ", cohens_d))
        } else {
          # If one group is missing, set result to NA
          cohens_d <- NA
          message("Only one group in AffectionStatus, Cohen's d set to NA")
        }
        stats <- c(
          ROI,
          N,
          mu,
          sdev,
          minV,
          maxV,
          minSubj,
          maxSubj,
          outliers,
          cohens_d
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

        hist(DATA, breaks = 20, main = paste0(cohort, ": ", ROI))
      } else {
        cat("ROI not found in the table:", ROI, "\n")
      }
    }
  }
  dev.off()
}
