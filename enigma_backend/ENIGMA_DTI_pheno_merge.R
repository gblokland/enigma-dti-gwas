# ENIGMA-DTI-GWAS
# Authors: Gabriella Blokland, Nina Roth Mota
# Last edit 2025 August

####################################################################
cmdargs = commandArgs(trailingOnly=T)
csvFILE_FA=cmdargs[1]
csvFILE_MD=cmdargs[2]
csvFILE_AD=cmdargs[3]
csvFILE_RD=cmdargs[4]
localfamFILE=cmdargs[5]
outFolder=cmdargs[6]
#################

ROIS=c("ACR","ACR.L","ACR.R","ALIC","ALIC.L","ALIC.R","AverageFA","BCC","CC","CGC","CGC.L","CGC.R",
       "CGH","CGH.L","CGH.R","CR","CR.L","CR.R","CST","CST.L","CST.R","EC","EC.L","EC.R",
       "FX","FX.ST.L","FX.ST.R","FXST","GCC","IC","IC.L","IC.R","IFO","IFO.L","IFO.R","PCR",
       "PCR.L","PCR.R","PLIC","PLIC.L","PLIC.R","PTR","PTR.L","PTR.R","RLIC","RLIC.L","RLIC.R","SCC",
       "SCR","SCR.L","SCR.R","SFO","SFO.L","SFO.R","SLF","SLF.L","SLF.R","SS","SS.L","SS.R","UNC","UNC.L","UNC.R")

ROIS2=c("FA","MD","AD","RD")

#FA <- data.frame(read.csv(csvFILE_FA))
#MD <- data.frame(read.csv(csvFILE_MD))

FA <- data.frame(read.csv(csvFILE_FA,colClasses = "character")) # if you have column names with unstandard symbols ("-") use this
MD <- data.frame(read.csv(csvFILE_MD,colClasses = "character"))
AD <- data.frame(read.csv(csvFILE_AD,colClasses = "character")) # if you have column names with unstandard symbols ("-") use this
RD <- data.frame(read.csv(csvFILE_RD,colClasses = "character"))

VarNamesFA=names(FA)
VarNamesMD=names(MD)
VarNamesAD=names(AD)
VarNamesRD=names(RD)

NsubjFA=dim(FA)[1]
NsubjMD=dim(MD)[1]
NsubjAD=dim(AD)[1]
NsubjRD=dim(RD)[1]

subjectNames=FA[,1]

if (NsubjFA != NsubjMD) {
  stop("Number of Subjects with FA measures do not match those with MD. Please make sure your files are correct.")
}
if (NsubjFA != NsubjAD) {
  stop("Number of Subjects with FA measures do not match those with AD. Please make sure your files are correct.")
}
if (NsubjFA != NsubjRD) {
  stop("Number of Subjects with FA measures do not match those with RD. Please make sure your files are correct.")
}
if (NsubjMD != NsubjAD) {
  stop("Number of Subjects with MD measures do not match those with AD. Please make sure your files are correct.")
}
if (NsubjMD != NsubjRD) {
  stop("Number of Subjects with MD measures do not match those with RD. Please make sure your files are correct.")
}
if (NsubjAD != NsubjRD) {
  stop("Number of Subjects with AD measures do not match those with RD. Please make sure your files are correct.")
}

fam <- data.frame(read.table(localfamFILE,colClasses = "character")) # if you have column names with unstandard symbols ("-") use this
names(fam) <- c("FID","IID","PID","MID","SEX","PHENO")

NsubjFam=dim(fam)[1]

if (NsubjFam != NsubjFA) {
       print(NsubjFam); print(NsubjFA)
  stop("Number of Subjects with genotype data do not match those with FA measures. Please make sure your files are correct.")
}


##THIS DOESN'T APPLY BECAUSE WE ALREADY HAVE AVERAGE MEASURES ACROSS L&R HEMISPHERE:
# Avg_ALL=cbind(subjectNames)
# colnames(Avg_ALL)=c("SubjID")
# 
# for (i in 1:length(ROIS)) {
#   ALLcolnames = colnames(Avg_ALL);
#   L_roi=paste0("FA_",ROIS[i],".L")
#   R_roi=paste0("FA_",ROIS[i],".R")
#   AVG_roi=paste("FA_",ROIS[i],".Mean",sep="")
#   #tmp = 0.5*(FA[,which(VarNamesFA==L_roi)]+FA[,which(VarNamesFA==R_roi)])
#   tmp = 0.5*(as.numeric(FA[,which(VarNamesFA==L_roi)])+ as.numeric(FA[,which(VarNamesFA==R_roi)])) # if you have column names with unstandard symbols ("-") use this
#   Avg_ALL = cbind(Avg_ALL, tmp)
#   colnames(Avg_ALL)<-c(ALLcolnames,AVG_roi)
#   
#   ALLcolnames = colnames(Avg_ALL);
#   L_roi=paste0("FA_",ROIS[i],".L")
#   R_roi=paste0("FA_",ROIS[i],".R")
#   AVG_roi=paste("FA_",ROIS[i],".Mean",sep="")
#   #tmp = 0.5*(MD[,which(VarNamesMD==L_roi)]+MD[,which(VarNamesMD==R_roi)])
#   tmp = 0.5*(as.numeric(MD[,which(VarNamesMD==L_roi)])+ as.numeric(MD[,which(VarNamesMD==R_roi)])) # if you have column names with unstandard symbols ("-") use this
#   Avg_ALL = cbind(Avg_ALL, tmp)
#   colnames(Avg_ALL)<-c(ALLcolnames,AVG_roi)
# }
# 
# for (i in 1:length(ROIS2)) {
#   if(i == 1){
#     ALLcolnames = colnames(Avg_ALL);
#     L_roi=paste0("L",ROIS2[i])
#     R_roi=paste0("R",ROIS2[i])
#     AVG_roi=paste0("Mean_Full_",ROIS2[i])
#     #tmp = FA[,which(VarNamesFA==L_roi)]+FA[,which(VarNamesFA==R_roi)]
#     tmp = as.numeric(FA[,which(VarNamesFA==L_roi)]) + as.numeric(FA[,which(VarNamesFA==R_roi)])
#     Avg_ALL = cbind(Avg_ALL, tmp)
#     colnames(Avg_ALL)<-c(ALLcolnames,AVG_roi)
#   } else {
#     ALLcolnames = colnames(Avg_ALL);
#     L_roi=paste("L",ROIS2[i],sep="")
#     R_roi=paste("R",ROIS2[i],sep="")
#     AVG_roi=paste("Mean_Full_",ROIS2[i],sep="")
#     #tmp = 0.5 * (FA[,which(VarNamesFA==L_roi)]+FA[,which(VarNamesFA==R_roi)])
#     tmp = 0.5 * (as.numeric(FA[,which(VarNamesFA==L_roi)]) + as.numeric(FA[,which(VarNamesFA==R_roi)]))
#     Avg_ALL = cbind(Avg_ALL, tmp)
#     colnames(Avg_ALL)<-c(ALLcolnames,AVG_roi)
#   }
# }
# 
# SummaryAvg_ALL = Avg_ALL[,-1]
# SummaryAvg_ALL_dat = colMeans(matrix(as.numeric(unlist(SummaryAvg_ALL)),nrow=nrow(SummaryAvg_ALL)), na.rm=T)
# 
# write.csv(SummaryAvg_ALL_dat,file=paste0(outFolder,"/SummaryMeasures.csv"),quote=F,row.names=F)
# write.csv(Avg_ALL,paste0(outFolder,"/DTI_Measures_ENIGMA_ALL_Avg.csv"),quote=F,row.names=F)

print(colnames(FA)); print(colnames(AD)); print(colnames(MD)); print(colnames(RD)); 
colnames(FA)[colnames(FA) == "AverageFA"] ="GlobalAverage"
colnames(AD)[colnames(AD) == "AverageFA"] ="GlobalAverage"
colnames(MD)[colnames(MD) == "AverageFA"] ="GlobalAverage"
colnames(RD)[colnames(RD) == "AverageFA"] ="GlobalAverage"
colnames(FA)[2:length(colnames(FA))] <- paste0("FA_", colnames(FA)[2:length(colnames(FA))])
colnames(AD)[2:length(colnames(AD))] <- paste0("AD_", colnames(AD)[2:length(colnames(AD))])
colnames(MD)[2:length(colnames(MD))] <- paste0("MD_", colnames(MD)[2:length(colnames(MD))])
colnames(RD)[2:length(colnames(RD))] <- paste0("RD_", colnames(RD)[2:length(colnames(RD))])
print(colnames(FA)); print(colnames(AD)); print(colnames(MD)); print(colnames(RD)); 

#DTI_ALL <- merge(RD, merge(MD, merge(FA, AD, by = "SubjID"), by = "SubjID"), by = "SubjID")
DTI_ALL <- merge(RD, merge(MD, merge(FA, AD, by = "subjectID"), by = "subjectID"), by = "subjectID")

#Rename/copy subjectID to FID and IID:
names(DTI_ALL)[names(DTI_ALL) == "subjectID"] <- "IID"
# Copy the variable
DTI_ALL$FID <- DTI_ALL$IID
# Move it to the first column
DTI_ALL <- DTI_ALL[, c("FID", names(DTI_ALL)[names(DTI_ALL) != "FID"])]

#write.csv(DTI_ALL,paste0(outFolder,"/DTI_Measures_ENIGMA_ALL.csv"),quote=F,row.names=F)
write.csv(DTI_ALL,paste0(outFolder,"/Phenotypes.csv"),quote=F,row.names=F)
writeLines(colnames(DTI_ALL))
writeLines(paste(colnames(DTI_ALL), collapse=";"))
