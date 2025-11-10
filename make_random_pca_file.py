#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

# --- Usage ---
if len(sys.argv) < 3:
    sys.exit("Usage: python make_random_pca_file.py <n_subjects> <output_file>")

n_subjects = int(sys.argv[1])
#n = 100  # number of subjects

out_file = sys.argv[2]

# --- Settings ---
n_pcs = 20
np.random.seed(42)

# --- Generate fake IDs ---
FID = np.arange(1, n_subjects+1)
IID = np.arange(1001, 1001+n_subjects)

# --- Generate random PCs ---
pcs = np.random.normal(loc=0.0, scale=0.03, size=(n_subjects, n_pcs))

# --- Combine ---
columns = ["FID", "IID"] + [f"PC{i+1}" for i in range(n_pcs)]
df_out = pd.DataFrame(np.column_stack([FID, IID, pcs]), columns=columns)

# --- Save ---
df_out.to_csv(out_file, sep="\t", index=False, header=True)
print(f"✅ Wrote fake eigenvec file: {out_file}")
print(f"   Columns: {', '.join(columns)}")
print(f"   Subjects: {n_subjects}")
