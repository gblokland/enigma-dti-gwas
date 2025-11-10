# * Written for ENIGMA-DTI-GWAS by Gabriella Blokland, Nina Roth Mota
# * Inspired by createDatPed.R and createDatPed_flexible_files_plasticity.R
# * Last edit 2025-10-31, by gblokland

###### Basic INFO

###################################################################
## to run:
# ${Rbin} --no-save --slave --args ${1} ${2} ...   ${10} <  ./.R
# R --no-save --slave --args ${csvFILE} ${localfamFILE} ${pcaFILE} ${combinedROItableFILE} ${ageColumnHeader} ${sexColumnHeader} ${maleIndicator} ${CaseControlCohort} ${affectionStatusColumnHeader} ${affectedIndicator} ${related} ${outDir} ${run_dir} ${eName} <  ${run_directory}/ENIGMA_DTI_create_formatted_files.R
###################################################################

# 14 INPUTS
######  1. path to your Covariates.csv file
######  2. path to your local *.fam file
######  3. path to your *_PCACovariates.eigenvec file
######  4. path to your combinedROItableFILE (merged output from TBSS pipeline) file: Phenotypes.csv
######  5. the column header for your age covariate
######  6. the column header for your sex covariate
######  7. what is the indicator for males in the sex column (M? 1? 2? ... )
######  8. does your dataset contain patients (0 for no, 1 for yes)
######  9. the column header for your AffectionStatus covariate
######  10. what is the indicator for affected individuals in the sex column (A? 1? 2? ... )
######  11. does your dataset contain related individuals (0 for no, 1 for yes)
######  12. output directory to write formatted files into (doesn't have to exist yet!)
######  list of ROIs to run
######  13. run_directory = SCRIPTS directory 
######  14. eName; used in output
######  15. Cohort Name; used in output

####################################################################

options(stringsAsFactors = FALSE)

# Output files will include .covar and .pheno files with all necessary covariates
# standard covariates include age, sex, age^2, age-x-sex, age^2-x-sex, and 4 MDS components
# for DTI additional covariates, depending on model, are GlobalAverageFA, GlobalAverageMD, GlobalAverageAD, GlobalAverageRD
###### sex and interaction terms not included if single-sex study
###### AffectionStatus included if patients are included
###### AffectionStatus2 etc. included if multiple patient groups are present
###### Any site specific covariates, including imaging site (if multiple) should be included.
###### Each site needs its own dummy covariate column
###### Each patient group needs its own dummy covariate column
####################################################################

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check and install argparse if not already installed
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Formatting script for ENIGMA DTI GWAS")

args <- commandArgs(trailingOnly = TRUE)
print(args)

####
csvFILE <- args[1]
localfamFILE <- args[2] # IS THIS NEEDED HERE? DO WE NEED TO CHECK MATCHING SUBJECTS AGAIN?
#Another function that greps from a fam file and the covariates and phenotype file to check whether the merge between genotype and phenotype worked, and no data are lost or mismatched that shouldn’t be. FID and IID columns should match (same ID format) and there shouldn’t be (many) more subjects in one file compared to the other.
####
pcaFILE <- args[3]
combinedROItableFILE <- args[4]
ageColumnHeader <- args[5]
sexColumnHeader <- args[6]

maleIndicator <- args[7]

if (is.na(as.numeric(maleIndicator)) == "TRUE") {
  maleIndicator=maleIndicator
} else {
  maleIndicator=as.numeric(maleIndicator)
}

CaseControlCohort <- args[8]  ## have a column where all healthy are marked as 0s and all patients as 1

if (is.na(as.numeric(CaseControlCohort)) == "TRUE") {
  CaseControlCohort=CaseControlCohort
} else {
  CaseControlCohort=as.numeric(CaseControlCohort)
}

affectionStatusColumnHeader <- args[9] 

affectedIndicator <- args[10] ## what is the number or letter or string used to identify an affected individual?

related <- args[11]  ## does your dataset contain related individuals (0 for no, 1 for yes)

if (is.na(as.numeric(related)) == "TRUE") {
  related <- related
} else {
  related <- as.numeric(related)
}

outDir <- args[12]
dir.create(outDir, showWarnings = FALSE)

#ALL_ROIS <- args[11]
#ALL_ROIS <- as.character(parse(text=ALL_ROIS))
ALL_ROIS <- c("ACR","ACR.L","ACR.R","ALIC","ALIC.L","ALIC.R","GlobalAverage","BCC","CC","CGC","CGC.L","CGC.R","CGH","CGH.L","CGH.R","CR","CR.L","CR.R","CST","CST.L","CST.R","EC","EC.L","EC.R","FX","FX.ST.L","FX.ST.R","FXST","GCC","IC","IC.L","IC.R","IFO","IFO.L","IFO.R","PCR","PCR.L","PCR.R","PLIC","PLIC.L","PLIC.R","PTR","PTR.L","PTR.R","RLIC","RLIC.L","RLIC.R","SCC","SCR","SCR.L","SCR.R","SFO","SFO.L","SFO.R","SLF","SLF.L","SLF.R","SS","SS.L","SS.R","UNC","UNC.L","UNC.R")
ALL_ROIS <- c(paste0("FA_", ALL_ROIS), 
              paste0("MD_", ALL_ROIS), 
              paste0("AD_", ALL_ROIS), 
              paste0("RD_", ALL_ROIS))

run_dir <- args[13]

eName <- args[14]

cohort <- args[15]

output_format <- "plink"

paste0(outDir, "/", eName, "_RUN_NOTES.txt")

#Open a text file for writing (by subsequent commands)
zz <- file(paste0(outDir, "/", eName, "_RUN_NOTES.txt"), "w")

####################################################################

source(paste0(run_dir, "/enigma_backend/ENIGMA_functions_DTI.R"))

####################################################################

#ALL_IDS=c("FID", "IID", "MID", "PID", "Sex");
ALL_IDS=c("FID", "IID", "MID", "PID");

Nids=length(ALL_IDS)
Nrois=length(ALL_ROIS)

####################################################################

###Read in the covariates file
covar <- data.frame(read.csv(csvFILE, colClasses = "character")) 

###Read in the ancestry principal components
pca <- read.table(pcaFILE, colClasses = "character")  #Read in the ancestry principal components
pca$SOL <- NULL; #Remove the “SOL” column in the MDS components since this is not a covariate to be included
colnames(pca) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

###Read in fam file
fam <- read.table(localfamFILE, colClasses = "character")
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
fam <- fam[,c("FID", "IID")]

###Read in the phenotypes file
pheno <- data.frame(read.csv(combinedROItableFILE, colClasses = "character")) 

#pheno$IID = pheno[,1] #This just renames a column for easier merging
pheno$subjectID = NULL
pheno$SubID_base = NULL
pheno$SubID_fup = NULL

###Merge the covariates and ROIs
cat('Merging your Phenotypes and Covariates Files...')
InfoFile <- merge(covar, pheno, by = c("FID", "IID"), all = TRUE)
cat('Done\n')

missing="";
l_missing=0;
for (i in 1:Nrois) {
  columnnames = colnames(InfoFile);
  if (length(InfoFile[,which(columnnames==ALL_ROIS[i])]) == 0) {
    missing=paste(missing,ALL_ROIS[i])
    l_missing=l_missing+1;
  }
}

if (l_missing > 0) {
  stop(paste("ERROR: You are missing the following ROIs:", missing ,". Please re-run latest ROI script!",sep=""))
}

###Merge the PCA and other covariates
cat('Merging your PCA files with your Phenotypes and Covariates Files...')
merged_temp <- merge(pca, InfoFile, by = c("FID", "IID"))
#print(head(merged_temp))
cat('Done\n')

### make sure subject names match up!
numsubjects = length(merged_temp[,1])
if (numsubjects==0){
  stop("ERROR: Please make sure your subjectID's in your phenotype csv file are the same as those listed in your PCA file")
}

### make sure there are no duplicates
dups <- duplicated(merged_temp$IID)
Ldups <- length(which(dups=="TRUE"))  ##**##** depending on R version this (& below) may have to be which(dups,"TRUE")
cat('    There are ',Ldups,' duplicate subjects.\n')
if (Ldups > 0){
  merged_temp <- merged_temp[-which(dups=="TRUE"),]
  print('    Duplicates have been removed.\n')
}

### make sure there are no subjects with missing covariates
cat("Subjects before incomplete cases are removed:", nrow(merged_temp), "\n")
# Remove rows with missing covariates
covariate_cols <- c(ageColumnHeader, sexColumnHeader, affectionStatusColumnHeader, "PC1")
merged_temp <- merged_temp[complete.cases(merged_temp[, covariate_cols]), ]
# Optional: report how many were dropped
cat("Subjects retained:", nrow(merged_temp), "\n")

numsubjects <- length(merged_temp$IID)

writeLines(paste0('VERSION: ',format(Sys.Date(),"%m/%d/%Y")), con=zz, sep="\n")

# check whether we should also write files for subjects below AND above 18
possible_subsets <- list(c(1,0,0), c("all","child","adult"), c(list(merged_temp,NULL,NULL)))

age <- as.numeric(merged_temp[,ageColumnHeader])
#print(age)
#min(age) < 18 & max(age >= 18)
if (min(age) < 18 & max(age >= 18)) {
  cat("Separating data file into children and adults\n")
  writeLines(paste('Separating data file into children and adults'), con=zz, sep="\n")
  num_child=dim(merged_temp[which(as.numeric(merged_temp[,ageColumnHeader]) < 18),])[1]
  possible_subsets[[3]][[2]]  <- merged_temp[which(as.numeric(merged_temp[,ageColumnHeader]) < 18),]
  possible_subsets[[1]][[2]] <- 1
  num_adult=dim(merged_temp[which(as.numeric(merged_temp[,ageColumnHeader]) >= 18),])[1]
  possible_subsets[[1]][[3]] <- 1
  possible_subsets[[3]][[3]] <- merged_temp[which(as.numeric(merged_temp[,ageColumnHeader]) >= 18),]
}	

for (s in c(1:3)) {
  if (possible_subsets[[1]][[s]] == 1) {
    suffix=possible_subsets[[2]][[s]]
    # select datafile 
    merged_temp <- possible_subsets[[3]][[s]]	
    cat('\n')	
    cat(sprintf("PROCESSING %s SUBJECTS",suffix))
    writeLines(paste(''))
    writeLines(sprintf('PROCESSING %s SUBJECTS',suffix), con=zz, sep="\n")
    columnnames = colnames(merged_temp)
    
    ### if no related individuals, create dummy paternal and maternal IDs
    # otherwise break to make sure these are entered
    if (related==0) {
      writeLines(paste('STUDY DESIGN: There are no related individuals.'), con=zz, sep="\n")
    } 
    if (related==1) {
      writeLines(paste('STUDY DESIGN: There are related individuals.'), con=zz, sep="\n")
    } 
    if ( (length(merged_temp[,which(columnnames=="MID")])==0 ) || (length(merged_temp[,which(columnnames=="PID")])	==0) ) {
      p=matrix(0, nrow=dim(merged_temp)[1], ncol=1)
      merged_temp$PID=p
      merged_temp$MID=p
      cat("Adding in artifical MID and PID\n")
    }
    
    ### Find age and sex columns, center, and create new age^2, age-x-sex and age^2-x-sex columns
    columnnames <- colnames(merged_temp)
    age=as.numeric(merged_temp[,which(columnnames==ageColumnHeader)])
    age_mean=mean(age)
    age_md=median(age)
    min_age=min(age)
    max_age=max(age)
    
    ageC=(age-age_mean)
    ageCsq=ageC*ageC
    
    ### Sanity check for missingness of sex variable
    #merged_temp[,sexColumnHeader]
    bad_vals <- merged_temp[,sexColumnHeader][!merged_temp[,sexColumnHeader] %in% c("M", "F") | is.na(merged_temp[,sexColumnHeader])]
    unique(bad_vals)
    
    ### What happens if there is missing sex? It assigns as female.
    sex=merged_temp[,which(columnnames==sexColumnHeader)]
    males=which(sex==maleIndicator)
    sexC=matrix(0,nrow=dim(merged_temp)[1],ncol=1)
    sexC[males]<- -0.5  
    sexC[-males]<- 0.5
    
    #Recode sex to match the needs of Plink: 1-2 coding
    StandardSex=data.frame("Sex"=sexC)
    StandardSex[which(sexC==-.5),]<- 1
    StandardSex[which(sexC==.5),]<- 2
    
    #make sure this also happens for case-control status
    
    Nm=length(which(sexC==-.5))
    Nf=length(which(sexC==.5))
    
    ## Print some basic stats on age and sex to the RUN_NOTES.txt
    writeLines(paste(''))
    writeLines(paste('Distribution stats for age in the selected group:'), con=zz, sep="\n")
    writeLines(paste('		mean:',age_mean), con=zz, sep="\n")
    writeLines(paste('		median:',age_md), con=zz, sep="\n")
    writeLines(paste('		min:',min_age), con=zz, sep="\n")
    writeLines(paste('		max:',max_age), con=zz, sep="\n")
    writeLines(paste('There are',Nm,'males in this study.'), con=zz, sep="\n")
    writeLines(paste('There are',Nf,'females in this study.'), con=zz, sep="\n")
    
    ##################### GB added #########################
    # TO DO:
    ## Print some basic stats on other covariates to the RUN_NOTES.txt
    writeLines(paste(''))
    writeLines(paste('Distribution stats for age in the selected group:'), con=zz, sep="\n")
    writeLines(paste('		mean:',age_mean), con=zz, sep="\n")
    writeLines(paste('		median:',age_md), con=zz, sep="\n")
    writeLines(paste('		min:',min_age), con=zz, sep="\n")
    writeLines(paste('		max:',max_age), con=zz, sep="\n")
    writeLines(paste('There are',Nm,'males in this study.'), con=zz, sep="\n")
    writeLines(paste('There are',Nf,'females in this study.'), con=zz, sep="\n")
    ##################### end GB added #########################
    
    age_sexC=age*sexC
    ageCsq_sexC=ageCsq*sexC
    
    ### Do not include sex or sex interaction variables if population is all M or all F or age stuff if all the same age
    if (sd(sexC) ==0) {
      cat("	WARNING: It appears this study (or subgroup) is of a single sex. If this is not the case, please check your Covariates.csv file")
      writeLines(paste('  WARNING: It appears this study (or subgroup) is of a single sex. If this is not the case, please check your Covariates.csv file.'), con=zz, sep="\n")
      age_sexC=age_sexC*0
      ageCsq_sexC=ageCsq_sexC*0
    }
    
    if (sd(ageC) ==0) {
      cat("	WARNING: It appears this study (or subgroup) is of a single age group. If this is not the case, please check your Covariates.csv file")
      writeLines(paste('  WARNING: It appears this study (or subgroup) is of a single age group. If this is not the case, please check your Covariates.csv file.'), con=zz, sep="\n")
      ageCsq=ageCsq*0
      ageC_sexC=ageC_sexC*0
      ageCsq_sexC=ageCsq_sexC*0
    }
    
    ## set columns as variables
    
    merged_temp_rest=as.data.frame(merged_temp)
    columnnames = colnames(merged_temp_rest);
    merged_temp_rest=merged_temp_rest[,-which(columnnames==ageColumnHeader)]
    columnnames = colnames(merged_temp_rest);
    merged_temp_rest=merged_temp_rest[,-which(columnnames==sexColumnHeader)]
    columnnames = colnames(merged_temp_rest);
    merged_temp_rest=merged_temp_rest[,-which(columnnames=="fup_duration")]
    columnnames = colnames(merged_temp_rest);
    merged_temp_rest=merged_temp_rest[,-which(columnnames=="FID")]
    columnnames = colnames(merged_temp_rest);
    merged_temp_rest=merged_temp_rest[,-which(columnnames=="IID")]
    columnnames = colnames(merged_temp_rest);
    if (length(merged_temp_rest[,which(columnnames=="MID")]) > 0) {
      merged_temp_rest=merged_temp_rest[,-which(columnnames=="MID")]
    }
    columnnames = colnames(merged_temp_rest);
    if (length(merged_temp_rest[,which(columnnames=="PID")]) > 0) {
      merged_temp_rest=merged_temp_rest[,-which(columnnames=="PID")]
    }
    
    # columnnames = colnames(merged_temp_rest);
    # if (length(merged_temp_rest[,which(columnnames=="zygosity")]) > 0) {
    # merged_temp_rest=merged_temp_rest[,-which(columnnames=="zygosity")]
    # }
    columnnames = colnames(merged_temp_rest);
    
    VarNames=colnames(merged_temp_rest)
    print(VarNames)
    
    FullInfoFile=cbind(merged_temp[,c('FID','IID','MID','PID')],StandardSex,sexC,age,ageCsq,age_sexC,ageCsq_sexC,merged_temp_rest)
    
    VarNames=colnames(FullInfoFile)
    print(VarNames)
    numsubjects = length(FullInfoFile$IID);
    
    nVar=dim(FullInfoFile)[2]
    Nset=Nids+Nrois
    nCov=nVar-Nset ### all the covariates that are left after removal of genetic-family columns and followup duration 
    #FullInfoFile=moveMe(FullInfoFile,ALL_ROIS,"after","Sex")
    FullInfoFile[,sexColumnHeader] <- NULL
    
    numsubjects = length(FullInfoFile$IID);
    
    if (sum(is.na(FullInfoFile$FID)) > 0 ) {
      FullInfoFile=FullInfoFile[-which(is.na(FullInfoFile$FID)),]
      merged_temp_rest=merged_temp_rest[-which(is.na(FullInfoFile$FID)),]
    }
    
    numsubjects = length(FullInfoFile$IID);
    VarNames=names(FullInfoFile)
    columnnames = colnames(FullInfoFile);
    
    drp=0
    
    # Remove covariates with sd = 0 keeping the patient columns if they exist
    
    for (l in (Nset+1):length(VarNames)) {
      columnnames = colnames(FullInfoFile);
      ###if (sd(as.numeric(FullInfoFile[,which(columnnames==VarNames[l])]),na.rm=T)==0) {
      if (!is.na(sd_val <- sd(as.numeric(FullInfoFile[, which(columnnames == VarNames[l])]), na.rm = TRUE)) && sd_val > 0) {
        if (CaseControlCohort==1 && length(grep("AffectionStatus",VarNames[l])) > 0) {  # No need to take into account possible multiple diagnoses, in this case AffectionStatus itself will have variance.
          next
        } else if (CaseControlCohort!=0 && VarNames[l] == CaseControlCohort ) {
          next
        }
        else {
          cat(paste('The standard deviation of column', VarNames[l], 'is zero. Therefore, the column will be removed.\n'))
          columnnames = colnames(FullInfoFile)
          #FullInfoFile=FullInfoFile[,-which(columnnames==VarNames[l])]
          FullInfoFile[,VarNames[l]] <- NULL  # remove the column safely
          drp=drp+1
        }
      }
    }	
    cat('Done\n')
    
    ## if diseases exist, make one .dat file without the AffectionStatus covariates
    ###### if when removing patients, all healthy individuals have the same value for a covariate, remove that as a covariate too
    
    FullInfoFile_healthy=FullInfoFile
    FullInfoFile_irrespective=FullInfoFile
    FullInfoFile_all_disease_corrected=FullInfoFile
    patients_covars=NULL
    
    if (CaseControlCohort!=0) {
      
      # remove patient columns for full group analysis 
      
      VarNames=colnames(FullInfoFile_irrespective)
      patients_covars=grep("AffectionStatus",VarNames)
      for (l in 1:length(patients_covars)) {
        FullInfoFile_irrespective=FullInfoFile_irrespective[,-which(colnames(FullInfoFile_irrespective)==VarNames[patients_covars[l]])]
      }
      
      VarNames=colnames(FullInfoFile_irrespective)
      columnnames = colnames(FullInfoFile_irrespective)
      
      for (l in (Nset+1):length(VarNames)) {
        #if (sd(FullInfoFile_irrespective[,which(columnnames==VarNames[l])])==0) {
        if (!is.na(sd_val <- sd(as.numeric(FullInfoFile_irrespective[, which(columnnames == VarNames[l])]), na.rm = TRUE)) && sd_val > 0) {
          cat(paste('For all individuals combined, the standard deviation of column', VarNames[l], 'is zero. Therefore, the column will be removed.\n'))
          FullInfoFile_irrespective=FullInfoFile_irrespective[,-which(colnames(FullInfoFile_irrespective)==VarNames[l])]
        }
      }
      cat('Done\n')
      
      # remove covariates without variance in the healthy group - TO FIX
      
      ##FullInfoFile_healthy <- subset(FullInfoFile_healthy,rowSums(data.frame(sapply(FullInfoFile_healthy[,patients_covars],as.numeric))[,,drop=FALSE])==0);
      
      FullInfoFile_healthy <- FullInfoFile_healthy[complete.cases(FullInfoFile_healthy[, patients_covars]), ]
      
      if (nrow(FullInfoFile_healthy) > 0) {
        VarNames=colnames(FullInfoFile_healthy)
        columnnames = colnames(FullInfoFile_healthy);
        for (l in (Nset+1):length(VarNames)){
          columnnames = colnames(FullInfoFile_healthy);
          #if (sd(FullInfoFile_healthy[,which(columnnames==VarNames[l])])==0) {
          if (!is.na(sd_val <- sd(as.numeric(FullInfoFile_healthy[, which(columnnames == VarNames[l])]), na.rm = TRUE)) && sd_val > 0) {
            cat(paste('For healthy individuals only, the standard deviation of column', VarNames[l], 'is zero. Therefore, the column will be removed.\n'))
            FullInfoFile_healthy=FullInfoFile_healthy[,-which(columnnames==VarNames[l])]
          }
        }
      }
      cat('Done\n')
      
    }
    cat('CHECKPOINT1\n')
    
    nVar=dim(FullInfoFile)[2] 
    nVar_healthy=dim(FullInfoFile_healthy)[2]
    nCov_healthy=nVar_healthy-Nset 
    nSub_healthy=dim(FullInfoFile_healthy)[1]
    nSub_patients=numsubjects - nSub_healthy
    
    nVar_irrespective=dim(FullInfoFile_irrespective)[2]
    nCov_irrespective=nVar_irrespective-Nset
    nSub_irrespective=dim(FullInfoFile_irrespective)[1]
    # only used when patients and controls are present 
    nCov_all_disease_corr=nCov_irrespective + length(patients_covars)
    nVar_disease_corr=dim(FullInfoFile_all_disease_corrected)[2]
    
    cat('CHECKPOINT2\n')
    
    # Ensure Nset and nVar are consistent with current data
    nVar <- ncol(FullInfoFile)
    if (Nset >= nVar) {
      warning("Nset is greater than or equal to number of columns in FullInfoFile.")
    }
    
    #test=NULL
    test <- data.frame(Variable = character(),
                       Mean = numeric(),
                       SD = numeric(),
                       Min = numeric(),
                       Max = numeric(),
                       stringsAsFactors = FALSE)
    
    #for (val in (Nset+1):nVar) {
    for (val in 5:nVar) {
      #A=as.numeric(FullInfoFile[,val])
      #if ( length(which(is.na(A))) > 0) {
      #  A=A[-which(is.na(A))]
      #}
      #test=rbind(test,cbind(colnames(FullInfoFile)[val],
      #                      (mean(as.numeric(A))),
      #                      (sd(as.numeric(A))),
      #                      (min(as.numeric(A))),
      #                      (max(as.numeric(A)))))
      
      # Defensive check
      if (val > ncol(FullInfoFile)) next
      
      A <- suppressWarnings(as.numeric(FullInfoFile[[val]]))
      A <- A[!is.na(A)]
      
      if (length(A) == 0) next  # skip empty columns
      
      test <- rbind(test, data.frame(
        Variable = colnames(FullInfoFile)[val],
        Mean = mean(A),
        SD = sd(A),
        Min = min(A),
        Max = max(A),
        stringsAsFactors = FALSE
      ))
      
    }
    header=(c('Covariate','Mean','SD','Min','Max'))
    write.table(test,file=zz,quote=F,col.names=header,row.names=FALSE,sep = "\t");
    
    cat('Step 3: Now write the formatted output files.\n')
    
    print(colnames(FullInfoFile))
    print(head(FullInfoFile))
    
    if (output_format== "plink") {
      cat('Writing PLINK formatted phenotype/covariate files.\n')
      write.table(FullInfoFile[, c("FID", "IID", ALL_ROIS)], paste0(outDir, "/", cohort, "_enigma_dti_gwas.pheno"),quote=F,col.names=F,row.names=F);
      write.table(FullInfoFile[, c("FID", "IID", "StandardSex","sexC","age","ageCsq","age_sexC","ageCsq_sexC", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")], paste0(outDir, "/", cohort, "_enigma_dti_gwas.covar"),quote=F,col.names=T,row.names=F);
    }
    
    if (output_format== "saige") {
      cat('Writing SAIGE formatted phenotype/covariate files.\n')
      write.table(FullInfoFile[, c("FID", "IID", ALL_ROIS)], paste0(outDir, "/", cohort, "_enigma_dti_gwas.pheno.saige.txt"),quote=F,col.names=F,row.names=F);
      write.table(FullInfoFile[, c("FID", "IID", "StandardSex","sexC","age","ageCsq","age_sexC","ageCsq_sexC", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")], paste0(outDir, "/", cohort, "_enigma_dti_gwas.covar.saige.txt"),quote=F,col.names=T,row.names=F);
    }
    
    if (output_format== "PedDat") {
      ################ now print out the .dat and the .ped files
      #
      # Output names have been hard-coded for formatting of follow-up scripts.
      #
      # Write out ped file only if the number of subjects is larger than the number of covariates
      
      #### print out one file irrespective of disease -- without covarying for AffectionStatus 
      #### this is all output in healthy and patient-groups-only cohorts
      
      if (numsubjects >= nCov_irrespective) {
        cat('    There are',numsubjects,'subjects irrespective of disease \n')
        cat('    There are',nCov_irrespective,'covariates for all subjects irrespective of disease\n')
        writeLines(paste('    There are',numsubjects,'subjects.'),con=zz,sep="\n")
        writeLines(paste('    There are',nCov_irrespective,'covariates for all subjects irrespective of disease.'),con=zz,sep="\n")
        writeLines(paste('     -', cbind(colnames(FullInfoFile_irrespective)[(Nset+1):nVar_irrespective])),con=zz,sep="\n")
        write.table(cbind(c(rep("T",Nrois),rep("C",nCov_irrespective)),c(colnames(FullInfoFile_irrespective)[(Nids+1):(nVar_irrespective)])),paste(outDir,"/ENIGMA_",eName,"_DATfile_",suffix,"_irrespective.dat",sep=""),col.names=F,row.names=F,quote=F);
        write.table(FullInfoFile_irrespective,paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_irrespective.ped",sep=""),quote=F,col.names=F,row.names=F);
        write.table(FullInfoFile_irrespective,paste(outDir,"/ENIGMA_",eName,"_PEDfile_wColNames_",suffix,"_irrespective.tbl",sep=""),quote=F,col.names=T,row.names=F);
        write.table(colnames(FullInfoFile_irrespective),paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_irrespective.header",sep=""),quote=F,col.names=F,row.names=F);
      } else {
        cat('	Not enough subjects for the ',suffix,' group, no files written.\n')
        writeLines(paste('    Not enough subjects for the ',suffix,' group, no files written.'),con=zz,sep="\n")
      }
      
      ####
      #### when both patients and controls are present, print out files for healthy only, and one including covariates for disease 
      if (nSub_healthy > 0 && nSub_patients > 0) {
        
        cat('    There are',nSub_irrespective,'subjects\n')
        cat('    There are',nCov_all_disease_corr,'covariates for all subjects correcting for disease\n')
        writeLines(paste('    There are',nSub_irrespective,'subjects.'),con=zz,sep="\n")
        writeLines(paste('    There are',nCov_all_disease_corr,'covariates for all subjects correcting for disease.'),con=zz,sep="\n")
        if (nSub_irrespective >= nCov_all_disease_corr) {
          writeLines(paste('     -', cbind(colnames(FullInfoFile_all_disease_corrected)[(Nset+1):nVar_disease_corr])),con=zz,sep="\n")
          write.table(FullInfoFile_all_disease_corrected,paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_mixedHD.ped",sep=""),quote=F,col.names=F,row.names=F);
          write.table(FullInfoFile_all_disease_corrected,paste(outDir,"/ENIGMA_",eName,"_PEDfile_wColNames_",suffix,"_mixedHD.tbl",sep=""),quote=F,col.names=T,row.names=F);
          write.table(colnames(FullInfoFile_all_disease_corrected),paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_mixedHD.header",sep=""),quote=F,col.names=F,row.names=F);
          write.table(cbind(c(rep("T",Nrois),rep("C",nCov_all_disease_corr)),c(colnames(FullInfoFile_all_disease_corrected)[(Nids+1):nVar_disease_corr])),paste(outDir,"/ENIGMA_",eName,"_DATfile_",suffix,"_mixedHD.dat",sep=""),col.names=F,row.names=F,quote=F);
        } else {
          cat('	Not enough subjects for the ',suffix,' group when controlling for disease, no files written.\n')
          writeLines(paste('    Not enough subjects for the ',suffix,' group when controlling for disease, no files written.'),con=zz,sep="\n")
        }
        
        cat('    There are',nSub_healthy,'healthy subjects\n')
        cat('    There are',nCov_healthy,'covariates for all healthy subjects\n')
        writeLines(paste('    There are',nSub_healthy,'healthy subjects.'),con=zz,sep="\n")
        writeLines(paste('    There are',nCov_healthy,'covariates for all healthy subjects.'),con=zz,sep="\n")
        if (nSub_healthy >= nCov_healthy) {
          writeLines(paste('     -', colnames(FullInfoFile_healthy)[(Nset+1):nVar_healthy]),con=zz,sep="\n")
          write.table(FullInfoFile_healthy,paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_healthy.ped",sep=""),quote=F,col.names=F,row.names=F);
          write.table(FullInfoFile_healthy,paste(outDir,"/ENIGMA_",eName,"_PEDfile_wColNames_",suffix,"_healthy.tbl",sep=""),quote=F,col.names=T,row.names=F);
          write.table(colnames(FullInfoFile_healthy),paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_healthy.header",sep=""),quote=F,col.names=F,row.names=F);
          write.table(cbind(c(rep("T",Nrois),rep("C",nCov_healthy)),c(colnames(FullInfoFile_healthy)[(Nids+1):nVar_healthy])),paste(outDir,"/ENIGMA_",eName,"_DATfile_",suffix,"_healthy.dat",sep=""),col.names=F,row.names=F,quote=F);
        } else {
          cat('Not enough subjects for the ',suffix,'healthy group, no files written.\n')
          writeLines(paste('Not enough subjects for the ',suffix,' healthy group, no files written.'),con=zz,sep="\n")
        }
        
        cat('    There are',nSub_patients,'disease subjects\n')
        cat('    There are',nCov_patients,'covariates for all disease subjects\n')
        writeLines(paste('    There are',nSub_patients,'disease subjects.'),con=zz,sep="\n")
        writeLines(paste('    There are',nCov_patients,'covariates for all disease subjects.'),con=zz,sep="\n")
        if (nSub_patients <= nCov_healthy) { #TO CHECK
          writeLines(paste('     -', colnames(FullInfoFile_patients)[(Nset+1):nVar_patients]),con=zz,sep="\n")
          write.table(FullInfoFile_patients,paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_disease.ped",sep=""),quote=F,col.names=F,row.names=F);
          write.table(FullInfoFile_patients,paste(outDir,"/ENIGMA_",eName,"_PEDfile_wColNames_",suffix,"_disease.tbl",sep=""),quote=F,col.names=T,row.names=F);
          write.table(colnames(FullInfoFile_patients),paste(outDir,"/ENIGMA_",eName,"_PEDfile_",suffix,"_disease.header",sep=""),quote=F,col.names=F,row.names=F);
          write.table(cbind(c(rep("T",Nrois),rep("C",nCov_patients)),c(colnames(FullInfoFile_patients)[(Nids+1):nVar_patients])),paste(outDir,"/ENIGMA_",eName,"_DATfile_",suffix,"_disease.dat",sep=""),col.names=F,row.names=F,quote=F);
        } else {
          cat('Not enough subjects for the ',suffix,'disease group, no files written.\n')
          writeLines(paste('Not enough subjects for the ',suffix,' disease group, no files written.'),con=zz,sep="\n")
        }
        
      } 
    } #end if (output_format== "PedDat")
  }
}

cat('****** DONE ****** Files have been written to ',outDir,'\n')
writeLines(paste('****** DONE ****** Files have been written. '),con=zz,sep="\n")
close(zz)


#The output file should be named: ${COHORT}_enigma_dti_gwas.pheno
#The output file should be named: ${COHORT}_enigma_dti_gwas.covar