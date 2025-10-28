#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                        %%%  ENIGMA DTI %%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% This is a function to print out plots and stats for Quality Control of ENIGMA-DTI measures 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% Written by Gabriella Blokland, based on script from Neda Jahanshad / Derrek Hibar
#%% Last update September 2025
#%% Questions or Comments??
#%% enigma.dtigenetics@gmail.com
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Check and install argparse if not already installed
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "QC script for ENIGMA DTI")
parser$add_argument("--site", required = TRUE, help = "Site name")
parser$add_argument("--outD", default = "./QC_ENIGMA/", help = "Output directory")
parser$add_argument("--CSVfile", default = "combinedROItable.csv", help = "Input CSV file")
parser$add_argument("--rois", default = "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R", help = "Semicolon-separated list of ROIs")
parser$add_argument("--outPDF", default = "ENIGMA_DTI_allROI_histograms.pdf", help = "Output PDF for histograms")
parser$add_argument("--outTXT", default = "ENIGMA_DTI_allROI_stats.txt", help = "Output TXT for statistics")

args <- parser$parse_args()

# Extract arguments
site <- args$site
outD <- args$outD
CSVfile <- args$CSVfile
rois <- args$rois
outPDF <- args$outPDF
outTXT <- args$outTXT

# Create output directory if it doesn't exist
dir.create(outD, showWarnings = FALSE)

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
  pdf(file = paste0(outD,"/",outPDF))
  write("Structure\tNumberIncluded\tMean\tStandDev\tMinValue\tMaxValue\tMinSubject\tMaxSubject\t5StDev_Off\tCohensD", file = paste0(outD,"/",outTXT))
  
  for (x in seq_len(Nrois)) {
    ROI <- parsedROIs[x]
    
    if (ROI %in% colnames(Table)) {
      writeLines(paste0("Current ROI = ", ROI))
      
      DATA <- as.numeric(as.vector(Table[,ROI]))
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
      }
      
      # Case-Control group
      group_var <- Table$group  # Assuming 'group' column exists and has values 'case' and 'control'
      group_case <- DATA[group_var == "case"]
      group_control <- DATA[group_var == "control"]
      
      # Calculate Cohen's d for case-control comparison
      cohens_d <- calculate_cohens_d(group_case, group_control)
      
      stats <- c(ROI, N, mu, sdev, minV, maxV, minSubj, maxSubj, outliers, cohens_d)
      write.table(t(as.matrix(stats)), file = paste0(outD,"/",outTXT), append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
      
      hist(DATA, breaks = 20, main = paste0(site, ": ", ROI))
      
    } else {
      cat("ROI not found in the table:", ROI, "\n")
    }
  }
  dev.off()
}

