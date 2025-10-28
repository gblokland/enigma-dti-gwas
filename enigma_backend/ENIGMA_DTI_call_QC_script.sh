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

Cohort=DMS
echo $Cohort
ENIGMA_DTI_GWAS_dir=~/enigmaDTI/enigma_dti_gwas
echo $ENIGMA_DTI_GWAS_dir
cd $ENIGMA_DTI_GWAS_dir

cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_call_QC_script.sh ./
cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL-GB.R ./
cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL.R ./
cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics.R ./
cp ~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics-simple-GB.R ./


wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_call_QC_script.sh"
wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL-GB.R"
wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_ALL.R"
wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics.R"
wget "~/enigma.dtigenetics\@gmail.com\ -\ Google\ Drive/My\ Drive/ENIGMA_DTI_GWAS_genetics_taskforce/ENIGMA_DTI_association_protocols/ENIGMA_DTI_plots_stats_demographics-simple-GB.R"


#Rscript ENIGMA_DTI_plots_qc_script.R \
Rscript ENIGMA_DTI_plots_ALL-GB.R \
--site ${Cohort} \
--outD "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/" \
--CSVfile "${ENIGMA_DTI_GWAS_dir}/combinedROItable_FA.csv" \
--rois "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR.L;ACR.R;ALIC.L;ALIC.R;CGC.L;CGC.R;CGH.L;CGH.R;CR.L;CR.R;CST.L;CST.R;EC.L;EC.R;FX.ST.L;FX.ST.R;IC.L;IC.R;IFO.L;IFO.R;PCR.L;PCR.R;PLIC.L;PLIC.R;PTR.L;PTR.R;RLIC.L;RLIC.R;SCR.L;SCR.R;SFO.L;SFO.R;SLF.L;SLF.R;SS.L;SS.R;UNC.L;UNC.R" \
--outPDF "${Cohort}_ENIGMA_DTI_allROI_histograms.pdf" \
--outTXT "${Cohort}_ENIGMA_DTI_allROI_stats.txt"

ls ${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/

#--rois "AverageFA;BCC;GCC;SCC;CC;CGC;CGH;CR;EC;FX;FXST;IC;IFO;PTR;SFO;SLF;SS;UNC;CST;ACR;ALIC;PCR;PLIC;RLIC;SCR;ACR-L;ACR-R;ALIC-L;ALIC-R;CGC-L;CGC-R;CGH-L;CGH-R;CR-L;CR-R;CST-L;CST-R;EC-L;EC-R;FX/ST-L;FX/ST-R;IC-L;IC-R;IFO-L;IFO-R;PCR-L;PCR-R;PLIC-L;PLIC-R;PTR-L;PTR-R;RLIC-L;RLIC.-R;SCR-L;SCR-R;SFO-L;SFO-R;SLF-L;SLF-R;SS-L;SS-R;UNC-L;UNC-R" \



Rscript ENIGMA_DTI_plots_stats_demographics.R \
--site ${Cohort} \
--outD "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/" \
--CSVfile "${ENIGMA_DTI_GWAS_dir}/Covariates.csv" \
--outPDF "${Cohort}_ENIGMA_DTI_Age_histograms.pdf" \
--outTXT "${Cohort}_ENIGMA_DTI_Age_stats.txt"


Rscript ENIGMA_DTI_plots_stats_demographics-simple-GB.R \
--site ${Cohort} \
--outD "${ENIGMA_DTI_GWAS_dir}/QC_ENIGMA/" \
--CSVfile "${ENIGMA_DTI_GWAS_dir}/Covariates.csv" \
--outPDF "${Cohort}_ENIGMA_DTI_Age_histograms.pdf" \
--outTXT "${Cohort}_ENIGMA_DTI_Age_stats.txt"

