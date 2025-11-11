
#R script that generates box plots for intracranial volume (ICV) by Sex, which will help visualize if the Sex variable is coded correctly. 
#The code assumes you have a Table frame with columns for ICV (intracranial volume) and Sex (the Sex variable, where it may be coded as "Male"/"Female", "1"/"2", or any other pattern).

# Set working directory for output
setwd("~/enigma/DTIgenetics")

# Load necessary libraries (install if not found)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Sample Table frame (replace this with your actual Table)
# Assuming 'ICV' is the intracranial volume and 'Sex' is the Sex variable.
# 'Sex' can be coded as "Male"/"Female", "1"/"2", or "M"/"F", etc.
#Table <- data.frame(
#  ICV = c(1400, 1600, 1550, 1500, 1350, 1700, 1450, 1600, 1550, 1450, 1750, 1650, 1800),
#  Sex = c("Male", "Female", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Female", "Male", "Female", "Male")
#)
Table <- read.csv("pheno_covar/ICV.csv", header = TRUE)
covars <- read.table("pheno_covar/DMS_enigma_dti_gwas.covar", header = TRUE)
Table <- merge(Table, covars, by = c("FID", "IID"), all = TRUE))
icvColumnHeader <- "segVOL_EstimatedTotalIntraCranialVol"
Table$ICV <- Table[,icvColumnHeader]

# Check unique values in the 'Sex' column
print(unique(Table$Sex))

# Clean and standardize 'Sex' variable based on detected coding
if (all(c("Male", "Female") %in% unique(Table$Sex))) {
  # Already properly coded
  Table$Sex <- factor(Table$Sex, levels = c("Male", "Female"), labels = c("Male", "Female"))
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
  Table$Sex <- factor(Table$Sex, levels = c("M", "F"), labels = c("Male", "Female"))
  message("Re-coded Sex: M → Male, F → Female")
} else {
  # Unexpected values
  message("Sex variable contains unexpected values — please check before recoding.")
}

# Create the plot
p <- ggplot(Table, aes(x = Sex, y = ICV, fill = Sex)) +
  geom_boxplot() +
  labs(title = "Intracranial Volume by Sex", x = "Sex", y = "Intracranial Volume (ICV)") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "red"))  # Optional: Customize the colors
# Show plot
p
# Save as PNG
ggsave(paste0(cohort, "_ICV_by_Sex.png"), plot = p, width = 6, height = 4, dpi = 300)
