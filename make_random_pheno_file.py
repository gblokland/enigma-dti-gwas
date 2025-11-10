import pandas as pd
import numpy as np

np.random.seed(42)
n = 100  # number of subjects

# IDs
FID = np.arange(1, n+1)
IID = np.arange(1001, 1001+n)

# Base structure
tracts = [
    "ACR", "ALIC", "BCC", "CC", "CGC", "CGH", "CR", "CST", "EC", "FX", "FX.ST",
    "FXST", "GCC", "IC", "IFO", "PCR", "PLIC", "PTR", "RLIC", "SCC", "SCR", "SFO",
    "SLF", "SS", "UNC"
]

# Metric templates
metrics = ["FA", "MD", "RD", "AD"]

# Generate all column names (both hemispheres and mean)
columns = []
for m in metrics:
    for t in tracts:
        columns.extend([
            f"{m}_{t}",
            f"{m}_{t}.L",
            f"{m}_{t}.R"
        ])
    columns.append(f"{m}_Average{m}")

# Function to generate plausible random DTI values per metric
def generate_metric_values(metric, size):
    if metric == "FA":
        return np.random.uniform(0.2, 0.8, size)
    elif metric == "MD":
        return np.random.uniform(0.0005, 0.0015, size)
    elif metric == "RD":
        return np.random.uniform(0.0003, 0.0010, size)
    elif metric == "AD":
        return np.random.uniform(0.0006, 0.0015, size)

# Build dataframe
data = {"FID": FID, "IID": IID}

for col in columns:
    metric = col.split("_")[0]  # FA, MD, RD, or AD
    data[col] = generate_metric_values(metric, n)

df = pd.DataFrame(data)

# Save as tab-separated .pheno file
df.to_csv("DMS_enigma_dti_gwas.pheno", sep="\t", index=False)
print("✅ File saved as DMS_enigma_dti_gwas.pheno")
print(df.head(3))
