#ENIGMA-DTI-GWAS 2025
#Author: Gabriella Blokland
#Last Modified: September 2025
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

roiLabels <- gsub("_", " ", parsedROIs)
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

dev.new()
p <- ggplot(TableLongMetric, aes(fill = AffectionStatus)) +
  geom_histogram(
    aes(x = Value, y = after_stat(count)),
    colour = "black",
    bins = 30
  ) +
  facet_wrap(~parsedROIs, ncol = 8, scales = "free_x") +
  xlim(range(0,1)) +
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
dev.off()
graphics.off()

}


#Histogram of Age range by AffectionStatus
dev.new()
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
dev.off()
graphics.off()

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
pdf(file = paste0(outDir, cohort, ".pdf"))

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

# Initialize results vector
cohens_d_results <- vector()

# Subset table for the current cohort - do we need this if there are Dummy variables?
#Table <- Table_orig[Table_orig$cohort == cohort, ]

# Count group sizes
n_CON <- nrow(Table[Table$AffectionStatus == "Control", ])
n_AFF <- nrow(Table[Table$AffectionStatus == "Affected", ])

# Run analysis only if both groups have sufficient data
if (n_CON > 5 & n_AFF > 5) {
  cohens_d <- vector()
  cohens_d_smd <- vector()
  cohens_d_cilb <- vector()
  cohens_d_ciub <- vector()
  cohens_d_se <- vector()
  variable <- vector()
  varlabel <- vector()
  colour <- vector()

  for (i in seq_along(parsedROIs)) {
    for (metric in c("FA", "MD", "AD", "RD")) {
    if (!all(is.na(Table[, paste0(metric, "_", parsedROIs[i])]))) {
      # Compute Cohen's d manually
      group_aff <- Table[Table$AffectionStatus == "Affected", paste0(metric, "_", parsedROIs[i])]
      group_con <- Table[Table$AffectionStatus == "Control", paste0(metric, "_", parsedROIs[i])]

      pooled_sd <- sqrt(
        (sum((group_con - mean(group_con, na.rm = TRUE))^2, na.rm = TRUE) +
          sum((group_aff - mean(group_aff, na.rm = TRUE))^2, na.rm = TRUE)) /
          (length(na.omit(group_con)) + length(na.omit(group_aff)) - 2)
      )

      cohens_d[i] <- (mean(group_aff, na.rm = TRUE) - mean(group_con, na.rm = TRUE)) / pooled_sd

      # t-test
      TT <- t.test(group_aff, group_con)
      N1 <- length(na.omit(group_aff))
      N2 <- length(na.omit(group_con))

      # Confidence intervals and SE
      ci <- ci.smd(ncp = TT$statistic, n.1 = N1, n.2 = N2, conf.level = .95)
      cohens_d_smd[i] <- ci$smd
      cohens_d_cilb[i] <- ci$Lower.Conf.Limit.smd
      cohens_d_ciub[i] <- ci$Upper.Conf.Limit.smd
      cohens_d_se[i] <- (cohens_d_ciub[i] - cohens_d_cilb[i]) / 3.92
    } else {
      cohens_d[i] <- NA
      cohens_d_smd[i] <- NA
      cohens_d_cilb[i] <- NA
      cohens_d_ciub[i] <- NA
      cohens_d_se[i] <- NA
    }

    variable[i] <- paste0(metric, "_", parsedROIs[i])
    varlabel[i] <- roiLabels[i]
    colour[i] <- roiColors[i]
    }
  }

  # Combine results into a dataframe
  cohens_d <- data.frame(
    colour,
    variable,
    varlabel,
    cohens_d_smd,
    cohens_d_cilb,
    cohens_d_ciub,
    cohens_d_se
  )

  cohens_d$variable <- factor(
    as.character(cohens_d$variable),
    levels = parsedROIs,
    labels = parsedROIs
  )
  cohens_d$colour <- factor(
    as.character(cohens_d$colour),
    levels = unique(roiColors),
    labels = unique(roiColors)
  )
  cohens_d$comparison <- paste("AFF-CON_", "All", sep = "")
  cohens_d$cohort <- cohort

  cohens_d_results <- rbind(cohens_d_results, cohens_d)

  # Plot
  dev.new()
  par(mar = c(0, 0, 0, 0))
  ggplot(cohens_d, aes(fill = colour)) +
    geom_bar(
      aes_string(x = "variable", y = "cohens_d_smd"),
      colour = "black",
      stat = "identity",
      position = "dodge"
    ) +
    geom_errorbar(
      aes(
        x = variable,
        ymin = cohens_d_smd - cohens_d_se,
        ymax = cohens_d_smd + cohens_d_se
      ),
      stat = "identity",
      position = "dodge",
      width = .25
    ) +
    xlab("") +
    ylab("Cohen's d Effect Size ± Standard Error") +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "white", colour = "grey90"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 10, face = "bold"),
      axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 10, face = "bold"),
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      legend.position = "none",
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.y = element_line(colour = "black"),
      axis.title.x = element_text(angle = 0, hjust = 0.5, colour = "black")
    ) +
    scale_x_discrete(expand = c(0, 0), labels = cohens_d$varlabel) +
    scale_y_continuous(limits = c(-1, 1), breaks = round(seq(-1, 1, 0.25), 2)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    annotate("text", x = nvar - 0.5, y = 0.7, label = paste("n AFF =", n_AFF)) +
    annotate("text", x = nvar - 1.5, y = 0.7, label = paste("n CON =", n_CON))

  #plot_fd <- paste0(outDir, cohort, "_cohensd_bar_graph_AFF-CON_", correction, ".png")
  plot_fd <- paste0(outDir, cohort, "_cohensd_bar_graph_AFF-CON.png")
  ggsave(filename = plot_fd, height = 24, width = 18, units = "cm", scale = 1, dpi = 600)
  dev.off()
  graphics.off()
}
