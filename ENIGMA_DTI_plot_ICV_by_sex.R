
#R script that generates box plots for intracranial volume (ICV) by sex, which will help visualize if the sex variable is coded correctly. 
#The code assumes you have a data frame with columns for ICV (intracranial volume) and sex (the sex variable, where it may be coded as "Male"/"Female", "1"/"2", or any other pattern).

# Set working directory for output
setwd("~/enigma/DTIgenetics")

# Load necessary libraries (install if not found)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Sample data frame (replace this with your actual data)
# Assuming 'ICV' is the intracranial volume and 'sex' is the sex variable.
# 'sex' can be coded as "Male"/"Female", "1"/"2", or "M"/"F", etc.
data <- data.frame(
  ICV = c(1400, 1600, 1550, 1500, 1350, 1700, 1450, 1600, 1550, 1450, 1750, 1650, 1800),
  sex = c("Male", "Female", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Female", "Male", "Female", "Male")
)

# Check unique values in the 'sex' column
print(unique(data$sex))

# Clean and standardize 'sex' variable if necessary
# For example, you can convert to "Male" and "Female" if coded differently
data$sex <- factor(data$sex, levels = c("Male", "Female"))

# Create the plot
p <- ggplot(data, aes(x = sex, y = ICV, fill = sex)) +
  geom_boxplot() +
  labs(
    title = "Intracranial Volume by Sex",
    x = "Sex",
    y = "Intracranial Volume (ICV)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "pink"))  # Optional: Customize the colors
# Show plot
p
# Save as PNG
ggsave("COHORTNAME_ICV_by_sex.png", plot = p, width = 6, height = 4, dpi = 300)
