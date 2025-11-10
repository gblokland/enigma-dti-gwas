import pandas as pd
import numpy as np

np.random.seed(42)  # reproducibility
n = 100  # number of subjects

# IDs
FID = np.arange(1, n+1)
IID = np.arange(1001, 1001+n)

# Demographics
Sex = np.random.choice([0,1], size=n)  # 0=Female, 1=Male
Age = np.random.randint(18, 81, size=n)
AgeCsq = Age ** 2
AgexSex = Age * Sex
AgeCsqxSex = AgeCsq * Sex
AffectionStatus = np.random.choice([0,1], size=n)

# PCs
PCs = {f'PC{i}': np.random.normal(0,1,n) for i in range(1,11)}

# DTI metrics
FA_AverageFA = np.random.uniform(0.4,0.7,n)
MD_AverageMD = np.random.uniform(0.0005,0.0015,n)
RD_AverageRD = np.random.uniform(0.0003,0.001,n)
AD_AverageAD = np.random.uniform(0.0006,0.0015,n)

# Dummy MR variables
Dummy1_MR = np.random.randint(0,5,n)
Dummy2_MR = np.random.randint(0,5,n)

# Assemble DataFrame
df = pd.DataFrame({
    'FID': FID,
    'IID': IID,
    'Sex': Sex,
    'Age': Age,
    'AgeCsq': AgeCsq,
    'AgexSex': AgexSex,
    'AgeCsqxSex': AgeCsqxSex,
    'AffectionStatus': AffectionStatus,
    **PCs,
    'FA_AverageFA': FA_AverageFA,
    'MD_AverageMD': MD_AverageMD,
    'RD_AverageRD': RD_AverageRD,
    'AD_AverageAD': AD_AverageAD,
    'Dummy1_MR': Dummy1_MR,
    'Dummy2_MR': Dummy2_MR
})

# Save as .covar file
df.to_csv('DMS_enigma_dti_gwas.covar', index=False, sep='\t')  # tab-separated for GWAS covar
print("File saved as DMS_enigma_dti_gwas.covar")
