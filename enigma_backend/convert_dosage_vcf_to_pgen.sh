#!/bin/bash
#Author: Gabriëlla Blokland for ENIGMA-DTI-GWAS working group
#Date: 2025-10-09

#To convert your (Michigan Imputation server) vcf files to pgen dosage files, run:
for i in {1..22}; do
plink2 --vcf ../chr"$i".dose.vcf.gz dosage=HDS --keep ${datafile}_QC4_keeplist.txt --maf 0.01 --extract-if-info "R2 >= 0.6" --min-ac 20 --make-pgen --out ${datafile}_chr"$i"
done 

#Use bcftools to extract PAR1, PAR2, and non-PAR regions of Chromosome X into separate vcf files:
bcftools index chrX.dose.vcf.gz

#PAR1 (GRCh37 positions):
bcftools view -r X:60001-2699520 chrX.dose.vcf.gz -Oz -o PAR1.dose.vcf.gz
#PAR2 (GRCh37 positions):
bcftools view -r X:154931044-155260560 chrX.dose.vcf.gz -Oz -o PAR2.dose.vcf.gz

#non-PAR (everything else on chrX):
# non-PAR before PAR1
bcftools view -r X:1-60000 chrX.dose.vcf.gz -Oz -o nonPAR_left.dose.vcf.gz
# non-PAR between PAR1 and PAR2
bcftools view -r X:2699521-154931043 chrX.dose.vcf.gz -Oz -o nonPAR_middle.dose.vcf.gz
# non-PAR after PAR2
bcftools view -r X:155260561-end chrX.dose.vcf.gz -Oz -o nonPAR_right.dose.vcf.gz

#You can merge them later with bcftools concat:

# Merge the PAR regions back together
bcftools concat -Oz -o chrX_PAR.dose.vcf.gz PAR1.dose.vcf.gz PAR2.dose.vcf.gz
# Index the PAR merged file
bcftools index chrX_PAR.dose.vcf.gz

# Merge the non-PAR regions back together
bcftools concat -Oz -o chrX_nonPAR.dose.vcf.gz nonPAR_left.dose.vcf.gz  nonPAR_middle.dose.vcf.gz nonPAR_right.dose.vcf.gz
# Index the non-PAR merged file
bcftools index chrX_nonPAR.dose.vcf.gz

#To convert your (Michigan Imputation server) vcf file for chrX_nonPAR to pgen dosage file, run:
plink2 --vcf ../chrX_nonPAR.dose.vcf.gz dosage=HDS --keep ${datafile}_QC4_keeplist.txt --maf 0.01 --extract-if-info "R2 >= 0.6" --min-ac 20 --make-pgen --out ${datafile}_chr23
#To convert your (Michigan Imputation server) vcf file for chrX_PAR to pgen dosage file, run:
plink2 --vcf ../chrX_PAR.dose.vcf.gz dosage=HDS --keep ${datafile}_QC4_keeplist.txt --maf 0.01 --extract-if-info "R2 >= 0.6" --min-ac 20 --make-pgen --out ${datafile}_chr25

#Calculate MAF from the filtered data for meta data file to accompany summary statistics:
for i in {1..22} 25; do
plink2 --pfile ${datafile}_chr"$i" --freq 'gz' --out ${datafile}_chr"$i" 
done
for i in 23; do
plink2 --pfile ${datafile}_chr"$i" --freqx 'gz' --out ${datafile}_chr"$i"
done

#Get INFO score and MAF from the VCF file to make meta data file to accompany summary statistics:
for i in {1..23} 25; do
#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\n' input.vcf.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/INFO\n' ${datafile}_chr"$i".vcf.gz | awk '{if ($5 >=0.01 || $5 <=0.99) print $1, $7 }' > info_scores_${datafile}_chr"$i".txt
done 


#OR IF YOU DID NOT USE MICHIGAN IMPUTATION SERVER, BUT IMPUTED LOCALLY:

#To convert your (Minimac4) dosage files to pgen dosage files, run:
for i in {1..23}; do
plink2 --import-dosage ../${datafile}_chr"$i".out.dosage.gz format=1 --keep ${datafile}_QC4_keeplist.txt --maf 0.01 --extract-col-con ${datafile}_chr"$i".info 7 --extract-col-cond-min 0.6 --make-pgen --out ${datafile}_chr"$i"
done 

for i in {1..23}; do
awk 'NR==1 { print $1, "INFO" } NR>1 {if ($5 >=0.01 || $5 <=0.99 ) print $1, $7 }' ../${datafile}_chr"$i".info > info_scores_${datafile}_chr"$i".txt
done 
