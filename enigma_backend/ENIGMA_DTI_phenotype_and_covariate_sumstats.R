#ENIGMA-DTI-GWAS 2025
#Author: Gabriella Blokland
#Last Modified: September 2025
#############################################################################################

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
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
suppressMessages(suppressWarnings(library(MBESS)))

# # get command line args (excluding defaults like --file)
# args <- commandArgs(trailingOnly = TRUE)
# # first argument = cohort name
# cohort <- args[1]
# # 2nd argument = output directory
# outDir <- args[2]
# # 3rd argument = Semicolon-separated list of ROIs
# rois <- "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R"

# Set up argument parser
parser <- ArgumentParser(description = "QC script for ENIGMA DTI")
parser$add_argument(
  "--cohort",
  required = TRUE,
  help = "Cohort name"
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

#Table <- read.csv("Table.csv") #For initial testing

cohort <- "DMS" # <-- change this as needed
rois <- "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R"
outDir <- "./QC_ENIGMA/"
eName <- "ENIGMA-DTI-GWAS"

covarFILE <- paste0(cohort, "_enigma_dti_gwas.covar")
phenoFILE <- paste0(cohort, "_enigma_dti_gwas.pheno")

covar <- read.table(covarFILE, header = TRUE)
pheno <- read.table(phenoFILE, header = TRUE)
Table <- merge(covar, pheno, by = c("FID", "IID"), all = TRUE)

# Parse ROIs
parsedROIs <- unlist(strsplit(rois, ";"))
Nrois <- length(parsedROIs)
writeLines(paste0("Nrois = ", Nrois))


roiLabels <- gsub("_", " ", parsedROIs)
roiColors <- setNames(rainbow(length(parsedROIs)), parsedROIs)


#Reshape from wide to long
TableLong <- reshape2(Table, direction = "long")

#Histograms of phenotypes by AffectionStatus
for (metric in c("FA", "MD", "AD", "RD")) {
TableLongMetric <- TableLong[grepl(metric, TableLong$variable), ]

dev.new()
p <- ggplot(TableLongMetric, aes(fill = AffectionStatus)) +
  geom_histogram(
    aes_string(x = value, y = "..count.."),
    colour = "black",
    size = 0.2,
    bins = 30
  ) +
  facet_wrap(~parsedROIs, ncol = 8, scales = "free_x") +
  xlim(range(Table[, paste0(metric, "_", parsedROIs[i])])) +
  theme_set(theme_bw()) +
  xlab(axis_label) +
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


#Histogram of age range
dev.new()
p <- ggplot(Table, aes(fill = AffectionStatus)) +
  geom_histogram(
    aes_string(x = Age, y = "..count.."),
    colour = "black",
    size = 0.2,
    bins = 30
  ) +
  xlim(range(Table[, "Age"])) +
  theme_set(theme_bw()) +
  xlab(axis_label) +
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
      size = 1
    ) +
    labs(x = feature, y = "Density")
  print(plt)
}

#Overlapping histograms
#https://stackoverflow.com/questions/6957549/overlaying-histograms-with-ggplot2-in-r
plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(
    df,
    aes(
      x = eval(parse(text = feature)),
      fill = eval(parse(text = label_column))
    )
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
      size = 1
    ) +
    labs(x = feature, y = "Density")
  plt + guides(fill = guide_legend(title = label_column))
}

#https://www.geeksforgeeks.org/how-to-export-multiple-plots-to-pdf-in-r/
# Open pdf file
pdf(file = paste0(cohort, ".pdf"))

## create a 2X2 grid
#par( mfrow= c(2,2) ) #Only needed when plotting to same page; if this is left out each plot is plotted on a different page

# draw plots
for (phenotype in c("ICV")) {
  plot_histogram(Table, phenotype)
  plot_multi_histogram(Table, phenotype, 'Sex')
  plot_multi_histogram(Table, phenotype, 'AffectionStatus')
}

if (is.numeric(Table_orig$SEX)) {
  Table_orig$SEX <- factor(
    Table_orig$SEX,
    levels = c(0, 1),
    labels = c("Female", "Male")
  )
}

# Check intracranial volume to see whether sex assignments are correct -> ICV should be larger for men.
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


# Descriptives: mean, median, standard deviation, min, max for each covariate

# Cohen’s d for case-control differences:
#   Is the difference visible? Decrease in FA, Increase in MD, AD, RD?
#   Is the difference in the correct direction, suggesting case-control status is coded correctly?

### Cohen's d - standardized mean difference
### Cohen's d can be calculated as the difference between the means divided by the pooled SD
### pooled SD = SQRT((SUM((subject_value - mean_value)^2)) / (total number of observations - number of AffectionStatuss))

# Initialize results vector
cohens_d_results <- vector()

# Subset table for the current cohort
Table <- Table_orig[Table_orig$cohort == cohort, ]

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
    if (!all(is.na(Table[, parsedROIs[i]]))) {
      # Compute Cohen's d manually
      group_aff <- Table[Table$AffectionStatus == "Affected", parsedROIs[i]]
      group_con <- Table[Table$AffectionStatus == "Control", parsedROIs[i]]

      pooled_sd <- sqrt(
        (sum((group_con - mean(group_con, na.rm = TRUE))^2, na.rm = TRUE) +
          sum((group_aff - mean(group_aff, na.rm = TRUE))^2, na.rm = TRUE)) /
          (length(na.omit(group_con)) + length(na.omit(group_aff)) - 2)
      )

      cohens_d[i] <- (mean(group_aff, na.rm = TRUE) -
        mean(group_con, na.rm = TRUE)) /
        pooled_sd

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

    variable[i] <- parsedROIs[i]
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
