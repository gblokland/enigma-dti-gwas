### ENIGMA-DTI-GWAS 2025
### This is a function to print out plots/stats for Quality Control 
### and Basic Statistics of ENIGMA-DTI Phenotypes and Covariates
###########################################################################
### Author: Gabriella Blokland
### Last update: October 2025
###########################################################################

options(repos = c(CRAN = "https://cloud.r-project.org/"))

#Check and install argparse if not already installed
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
suppressMessages(suppressWarnings(library(argparse)))
if (!requireNamespace("gtools", quietly = TRUE)) {
  install.packages("gtools")
}
suppressMessages(suppressWarnings(library(gtools)))
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
suppressMessages(suppressWarnings(library(ggplot2)))
if (!requireNamespace("MBESS", quietly = TRUE)) {
  install.packages("MBESS")
}
suppressMessages(suppressWarnings(library(MBESS)))
if (!requireNamespace("reshape", quietly = TRUE)) {
  install.packages("reshape")
}
suppressMessages(suppressWarnings(library(reshape)))
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
suppressMessages(suppressWarnings(library(dplyr)))
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
suppressMessages(suppressWarnings(library(stringr)))


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

covar <- read.table(covarFILE, header = TRUE)
covar <- covar[ , !(names(covar) %in% c("FA_AverageFA", "MD_AverageMD", "RD_AverageRD", "AD_AverageAD"))]
pheno <- read.table(phenoFILE, header = TRUE)
Table <- merge(covar, pheno, by = c("FID", "IID"), all = TRUE)
#print(head(Table))

# Parse ROIs
parsedROIs <- unlist(strsplit(rois, ";"))
Nrois <- length(parsedROIs)

writeLines(paste0("Nrois = ", Nrois))
writeLines(paste0("ROIs: ", parsedROIs))

roiLabels <- gsub("\\.L", " Left", gsub("\\.R", " Right", parsedROIs))
roiColors <- setNames(rainbow(length(parsedROIs)), parsedROIs)

Table$AffectionStatus <- factor(Table$AffectionStatus, 
  levels = c(1, 0), 
  labels = c("Affected", "Control"))

# Identify varying vs non-varying columns
varying_cols <- grep("^(FA|MD|AD|RD)_", names(Table), value = TRUE)
id_cols <- setdiff(names(Table), varying_cols)

# Reshape from wide → long
TableLong <- reshape(
  Table,
  varying = varying_cols,
  v.names = "Value",
  timevar = "Variable",
  times = varying_cols,
  idvar = id_cols,
  direction = "long"
)
#print(head(TableLong))
TableLong <- TableLong %>%
  mutate(
    Metric = str_extract(Variable, "^(FA|MD|AD|RD)"),
    Tract = str_extract(Variable, "(?<=_)[A-Z]+"),
    Side = str_extract(Variable, "(L|R)$"),
    parsedROIs = paste0(Tract, ifelse(is.na(Side), "", paste0(".", Side)))
  )

#Histograms of phenotypes by AffectionStatus
for (metric in c("FA", "MD", "AD", "RD")) {
TableLongMetric <- TableLong[grepl(paste0("^", metric), TableLong$Variable), ]
rownames(TableLongMetric) <- 1:nrow(TableLongMetric)
print(head(TableLongMetric))

# FA range (unitless)
FA_range <- c(0, 1)
# MD, RD, AD ranges (×10^-3 mm²/s)
MD_range <- c(0.5e-3, 1.0e-3)
RD_range <- c(0.3e-3, 0.8e-3)
AD_range <- c(1.0e-3, 1.75e-3)

# Map metric to range
range_map <- list(
  FA = FA_range,
  MD = MD_range,
  RD = RD_range,
  AD = AD_range
)

# Select the correct range
x_limits <- range_map[[metric]]

#dev.new() #Causes empty Rplots.pdf to get written
p <- ggplot(TableLongMetric, aes(fill = AffectionStatus)) +
  geom_histogram(
    aes(x = Value, y = after_stat(count)),
    colour = "black",
    bins = 30
  ) +
  facet_wrap(~parsedROIs, ncol = 8, scales = "free_x") +
  xlim(x_limits) +
  theme_set(theme_bw()) +
  xlab("Phenotype Value") +
  ylab("Number of Subjects") +
  ggtitle(metric) +
  theme(
    strip.background = element_rect(fill = "grey80", colour = "grey80"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "black", angle = 0),
    panel.background = element_rect(fill = "white", colour = "grey80"),
    panel.border = element_rect(colour = "grey80"),
    plot.title = element_text(angle = 0, hjust = 0.5, colour = "black", size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 10),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = "top",
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  #scale_fill_brewer(palette = "Spectral") +
  #scale_x_continuous(limits = x_limits, breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(expand = c(0, 0))
plot_fd <- paste0(outDir, cohort, "_histogram_multi-panel_", metric, "_by_AffectionStatus.pdf")
ggsave(filename = plot_fd, p, height = 24, width = 36, units = "cm", scale = 1, dpi = 600)
#dev.off()
#graphics.off()

}


#Histogram of Age range by AffectionStatus
#dev.new() #Causes empty Rplots.pdf to get written
p <- ggplot(Table, aes(fill = AffectionStatus)) +
  geom_histogram(
    aes(x = Age, y = after_stat(count)),
    colour = "black",
    linewidth = 0.2,
    bins = 30
  ) +
  xlim(range(Table[, "Age"])) +
  theme_set(theme_bw()) +
  xlab("Age (Years)") +
  ylab("Number of Subjects") +
  ggtitle("Age") +
  theme(
    strip.background = element_rect(fill = "grey80", colour = "grey80"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "black", angle = 0),
    panel.background = element_rect(fill = "white", colour = "grey80"),
    panel.border = element_rect(colour = "grey80"),
    plot.title = element_text(angle = 0, hjust = 0.5, colour = "black", size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 10),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = "top",
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  #scale_fill_brewer(palette = "Spectral") +
  #scale_x_continuous(limits = x_limits, breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(expand = c(0, 0))
plot_fd <- paste0(outDir, cohort, "_histogram_Age.pdf")
ggsave(filename = plot_fd, p, height = 24, width = 36, units = "cm", scale = 1, dpi = 600)
#dev.off()
#graphics.off()

#Single histogram:
plot_histogram <- function(df, feature) {
  plt <- ggplot(df, aes(x = eval(parse(text = feature)))) +
    geom_histogram(
      aes(y = ..density..),
      alpha = 0.7,
      fill = "#33AADE",
      color = "black"
    ) +
    geom_density(alpha = 0.3, fill = "red") +
    geom_vline(
      aes(xintercept = mean(eval(parse(text = feature)))),
      color = "black",
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(x = feature, y = "Density")
  print(plt)
}

#Overlapping histograms
#https://stackoverflow.com/questions/6957549/overlaying-histograms-with-ggplot2-in-r
plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(
    df,
    aes(x = eval(parse(text = feature)), fill = eval(parse(text = label_column)))
  ) +
    geom_histogram(
      alpha = 0.7,
      position = "identity",
      aes(y = ..density..),
      color = "black"
    ) +
    geom_density(alpha = 0.7) +
    geom_vline(
      aes(xintercept = mean(eval(parse(text = feature)))),
      color = "black",
      linetype = "dashed",
      slinewidth = 1
    ) +
    labs(x = feature, y = "Density")
  plt + guides(fill = guide_legend(title = label_column))
}

#https://www.geeksforgeeks.org/how-to-export-multiple-plots-to-pdf-in-r/
# Open pdf file
#pdf(file = paste0(outDir, cohort, ".pdf"))

## create a 2X2 grid
#par( mfrow= c(2,2) ) #Only needed when plotting to same page; if this is left out each plot is plotted on a different page

# Check intracranial volume to see whether sex assignments are correct 
# -> ICV should be larger for men.
# draw plots
if ("ICV" %in% names(Table)) {
  for (phenotype in c("ICV")) {
    plot_histogram(Table, phenotype)
    plot_multi_histogram(Table, phenotype, 'Sex')
    plot_multi_histogram(Table, phenotype, 'AffectionStatus')
  }

#if (is.numeric(Table$Sex)) {
#  Table$Sex <- factor(
#    Table$Sex,
#    levels = c(0, 1),
#    labels = c("Female", "Male")
#  )
#}

# Check unique values in the 'Sex' column
print(unique(Table$Sex))

# Clean and standardize 'Sex' variable based on detected coding
if (all(c("Male", "Female") %in% unique(Table$Sex))) {
  # Already properly coded
  Table$Sex <- factor(
    Table$Sex,
    levels = c("Male", "Female"),
    labels = c("Male", "Female")
  )
} else if (all(c(1, 2) %in% unique(Table$Sex))) {
  # Recode 1 = Male, 2 = Female
  Table$Sex <- factor(Table$Sex, levels = c(1, 2), labels = c("Male", "Female"))
  message("Re-coded Sex: 1 → Male, 2 → Female")
} else if (all(c(0, 1) %in% unique(Table$Sex))) {
  # Recode 0 = Female, 1 = Male (common convention)
  Table$Sex <- factor(Table$Sex, levels = c(1, 0), labels = c("Male", "Female"))
  message("Re-coded Sex: 1 → Male, 0 → Female")
} else if (all(c("M", "F") %in% unique(Table$Sex))) {
  # Recode F = Female, M = Male (common convention)
  Table$Sex <- factor(
    Table$Sex,
    levels = c("M", "F"),
    labels = c("Male", "Female")
  )
  message("Re-coded Sex: M → Male, F → Female")
} else {
  # Unexpected values
  message(
    "Sex variable contains unexpected values — please check before recoding."
  )
}

# Create the plot
p <- ggplot(Table, aes(x = Sex, y = ICV, fill = Sex)) +
  geom_boxplot() +
  labs(
    title = "Intracranial Volume by Sex",
    x = "Sex",
    y = "Intracranial Volume (ICV)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "red")) # Optional: Customize the colors
# Save as PNG
ggsave(
  paste0(cohort, "_ICV_by_Sex.png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 600
)

} #end if ("ICV" %in% names(Table)) {

# Descriptives: mean, median, standard deviation, min, max for each covariate

# Cohen’s d for case-control differences:
#   Is the difference visible? Decrease in FA, Increase in MD, AD, RD?
#   Is the difference in the correct direction, suggesting case-control status is coded correctly?

### Cohen's d - standardized mean difference
### Cohen's d can be calculated as the difference between the means divided by the pooled SD
### pooled SD = SQRT((SUM((subject_value - mean_value)^2)) / (total number of observations - number of AffectionStatuss))

# Initialize results dataframe
cohens_d_results <- data.frame()

# Count group sizes
n_CON <- sum(Table$AffectionStatus == "Control", na.rm = TRUE)
n_AFF <- sum(Table$AffectionStatus == "Affected", na.rm = TRUE)

# Run analysis only if both groups have sufficient data
if (n_CON > 5 & n_AFF > 5) {

  # Loop over metrics and ROIs
  for (metric in c("FA", "MD", "AD", "RD")) {
    for (i in seq_along(parsedROIs)) {
      col_name <- paste0(metric, "_", parsedROIs[i])
      
      if (!all(is.na(Table[[col_name]]))) {
        group_aff <- Table$AffectionStatus == "Affected"
        group_con <- Table$AffectionStatus == "Control"
        
        x_aff <- Table[group_aff, col_name]
        x_con <- Table[group_con, col_name]
        
        # Pooled SD
        pooled_sd <- sqrt(
          (sum((x_aff - mean(x_aff, na.rm = TRUE))^2, na.rm = TRUE) +
             sum((x_con - mean(x_con, na.rm = TRUE))^2, na.rm = TRUE)) /
            (length(na.omit(x_aff)) + length(na.omit(x_con)) - 2)
        )
        
        # Cohen's d
        d <- (mean(x_aff, na.rm = TRUE) - mean(x_con, na.rm = TRUE)) / pooled_sd
        
        # t-test
        TT <- t.test(x_aff, x_con)
        N1 <- length(na.omit(x_aff))
        N2 <- length(na.omit(x_con))
        
        # Confidence intervals and SE (ci.smd from MBESS)
        ci <- ci.smd(ncp = TT$statistic, n.1 = N1, n.2 = N2, conf.level = 0.95)
        
      } else {
        d <- NA
        ci <- list(smd = NA, Lower.Conf.Limit.smd = NA, Upper.Conf.Limit.smd = NA)
      }
      
      # Combine results into a row
      row <- data.frame(
        colour = roiColors[i],
        metric = metric,
        variable = col_name,
        varlabel = roiLabels[i],
        cohens_d_smd = ci$smd,
        cohens_d_cilb = ci$Lower.Conf.Limit.smd,
        cohens_d_ciub = ci$Upper.Conf.Limit.smd,
        cohens_d_se = (ci$Upper.Conf.Limit.smd - ci$Lower.Conf.Limit.smd)/3.92
      )
      
      cohens_d_results <- rbind(cohens_d_results, row)
    }
  }
  
  # Factor levels
  cohens_d_results$variable <- factor(cohens_d_results$variable, levels = cohens_d_results$variable)
  cohens_d_results$colour <- factor(cohens_d_results$colour, levels = unique(cohens_d_results$colour))
  cohens_d_results$comparison <- "AFF-CON"
  cohens_d_results$cohort <- cohort
  
  # Save results as csv
  write.csv(cohens_d_results,
          file = paste0(outDir, cohort, "_cohensd_results.csv"),
          row.names = FALSE)

cohens_d_results <- cohens_d_results %>%
  mutate(varlabel = factor(varlabel, levels = unique(varlabel[order(variable)])))

  # Plot
  p <- ggplot(cohens_d_results, aes(x = varlabel, y = cohens_d_smd, fill = colour)) +
    #geom_bar(stat = "identity", colour = "black", position = "dodge") +
    geom_col(colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = cohens_d_smd - cohens_d_se,
                      ymax = cohens_d_smd + cohens_d_se),
                  width = 0.25, position = position_dodge(0.9)) +
    coord_flip() +
    facet_wrap(~metric, ncol = 2, scales = "free_y") +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    #scale_x_discrete(labels = cohens_d_results$varlabel) +
    geom_hline(yintercept = 0) +
    xlab("") + ylab("Cohen's d ± SE") +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "white", colour = "grey90"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 6, face = "bold"),
      axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 6, face = "bold"),
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      legend.position = "none"
    ) +
    annotate("text", x = nrow(cohens_d_results)/4 - 2, y = 0.7, label = paste("n AFF =", n_AFF), size = 3.5) +
    annotate("text", x = nrow(cohens_d_results)/4 - 5, y = 0.7, label = paste("n CON =", n_CON), size = 3.5)
  
  # Save
  plot_fd <- paste0(outDir, cohort, "_cohensd_bar_graph_AFF-CON.png")
  ggsave(plot_fd, p, height = 24, width = 18, units = "cm", dpi = 600)
}
