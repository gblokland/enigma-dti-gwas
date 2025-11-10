#!/bin/bash
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                        %%%  ENIGMA DTI %%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#QC of extracted DTI phenotypes after merging with the genetic data
#%% This is a wrapper script for the R function to print out images for 
#%% Quality Control of DTI_ENIGMA FA images with TBSS (FSL) skeletons overlaid as well as JHU atlas ROIs
#%% Please QC your images to make sure they are correct FA maps and oriented and aligned properly
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% Written by Gabriella Blokland
#%% Last update December 2024
#%% Questions or Comments??
#%% enigma.dtigenetics@gmail.com
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COHORT=DMS
echo $COHORT
#ENIGMA_DTI_GWAS_dir=~/enigmaDTI/enigma_dti_gwas
ENIGMA_DTI_GWAS_dir=~/enigma/DTIgenetics
echo $ENIGMA_DTI_GWAS_dir
cd $ENIGMA_DTI_GWAS_dir

pheno_covar_dir=$ENIGMA_DTI_GWAS_dir/pheno_covar

#cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_call_QC_script.sh ./
#cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL-GB.R ./
#cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL.R ./
#cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics.R ./
#cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics-simple-GB.R ./


#wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_call_QC_script.sh"
#wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL-GB.R"
#wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL.R"
#wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics.R"
#wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics-simple-GB.R"

rm -rf SCRIPTS; git clone git@github.com:gblokland/enigma-dti-gwas.git SCRIPTS
run_directory=$HOME/enigma/DTIgenetics/SCRIPTS


csvFILE_FA="${pheno_covar_dir}/combinedROItable_FA.csv"
csvFILE_MD="${pheno_covar_dir}/combinedROItable_MD.csv"
csvFILE_AD="${pheno_covar_dir}/combinedROItable_AD.csv"
csvFILE_RD="${pheno_covar_dir}/combinedROItable_RD.csv"
#localfamFILE="$HOME/enigma/DTIgenetics/${datafile}_QC4.fam"
localfamFILE="$HOME/enigma/DTIgenetics/DMS_EUR_GB_20251031_QC4.fam"
outFolder=$pheno_covar_dir


R --no-save --slave --args ${csvFILE_FA} ${csvFILE_MD} ${csvFILE_AD} ${csvFILE_RD} ${localfamFILE} ${outFolder} < ${run_directory}/enigma_backend/ENIGMA_DTI_pheno_merge.R

## Specify the inputs:
csvFILE="${pheno_covar_dir}/Covariates.csv" # modify if different
#localfamFILE="$HOME/enigma/DTIgenetics/${datafile}_QC4.fam" # Do we need another sanity check to make sure numbers of subjects match up?
localfamFILE="$HOME/enigma/DTIgenetics/DMS_EUR_GB_20251031_QC4.fam"
combinedROItableFILE="${pheno_covar_dir}/Phenotypes.csv" # Needed to get global DTI metrics
#pcaFILE="${pheno_covar_dir}/${datafile}_PCACovariates.eigenvec" # modify if the path is different
pcaFILE="${pheno_covar_dir}/DMS_PCACovariates.eigenvec" # modify if the path is different
pcaFILE=$pcaFILE #specified location above
#ageColumnHeader="Age" #Column header for your age covariate (Case sensitive); modify if different
ageColumnHeader="Age_MRI" #Column header for your age covariate (Case sensitive); modify if different
#sexColumnHeader="Sex" #Column header for your sex covariate (Case sensitive); modify if different
sexColumnHeader="SEX" #Column header for your sex covariate (Case sensitive); modify if different
maleIndicator="male" #maleIndicator="M" #What is the indicator for males in the sex column (M? 1?); modify if different. ### Please be EXTRA careful if using a numeric indicator: Triple-check that the indicator code is correct! CaseControlCohort="0" #Does your cohort have a case-control or a case-only design? (mark 0=population-based or control-only cohort and 1=case-control or case-only cohort);
CaseControlCohort="1"
affectionStatusColumnHeader="AffectionStatus"; #Column header for your affection status covariate (Case sensitive); modify if different
affectedIndicator="1" #affectedIndicator="A" #What is the indicator for affected (i.e. cases, patients) individuals?; modify if different
related="0" #0 for unrelated cohorts; 1 for related cohorts
outDir="${ENIGMA_DTI_GWAS_dir}/pheno_covar"
run_dir="${run_directory}"
eName="ENIGMA_DTI_GWAS"

echo "Inputs":
echo ${csvFILE} ${localfamFILE} ${pcaFILE} ${combinedROItableFILE} ${ageColumnHeader} ${sexColumnHeader} ${maleIndicator} ${CaseControlCohort} ${affectionStatusColumnHeader} ${affectedIndicator} ${related} ${outDir} ${run_dir} ${eName} ${run_directory}


R --no-save --slave --args ${csvFILE} ${localfamFILE} ${pcaFILE} ${combinedROItableFILE} ${ageColumnHeader} ${sexColumnHeader} ${maleIndicator} ${CaseControlCohort} ${affectionStatusColumnHeader} ${affectedIndicator} ${related} ${outDir} ${run_directory} ${eName} <  ${run_directory}/enigma_backend/ENIGMA_DTI_create_formatted_files.R


###


ls ${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/

Rscript ${ENIGMA_DTI_GWAS_dir}/SCRIPTS/enigma_backend/ENIGMA_DTI_plots_All_inclCohensd.R \
--cohort ${COHORT} \
--covarFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.covar" \
--phenoFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.pheno" \
--rois "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R" \
--outDir "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/" \
--outPDF "${COHORT}_ENIGMA_DTI_allROI_histograms.pdf" \
--outTXT "${COHORT}_ENIGMA_DTI_allROI_stats.txt"

Rscript ${ENIGMA_DTI_GWAS_dir}/SCRIPTS/enigma_backend/ENIGMA_DTI_plots_stats_demographics.R \
--cohort ${COHORT} \
--covarFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.covar" \
--phenoFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.pheno" \
--outDir "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/" \
--outPDF "${COHORT}_ENIGMA_DTI_Age_histograms.pdf" \
--outTXT "${COHORT}_ENIGMA_DTI_Age_stats.txt"

Rscript ${ENIGMA_DTI_GWAS_dir}/SCRIPTS/enigma_backend/ENIGMA_DTI_phenotype_and_covariate_sumstats.R \
--cohort ${COHORT} \
--covarFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.covar" \
--phenoFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.pheno" \
--outDir "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/" \


Rscript ${ENIGMA_DTI_GWAS_dir}/SCRIPTS/enigma_backend/ENIGMA_DTI_GWAS_QC_plots_stats.R \
--cohort ${COHORT} \
--covarFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.covar" \
--phenoFILE "${ENIGMA_DTI_GWAS_dir}/pheno_covar/${COHORT}_enigma_dti_gwas.pheno" \
--outDir "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/"



