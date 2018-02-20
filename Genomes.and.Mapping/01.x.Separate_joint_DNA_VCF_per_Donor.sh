#!/usr/bin/env bash
'''
### 01.x.Separate_joint_DNA_VCF_per_Donor.sh
- DNA SNPs were jointly called for all samples.
- This script separates them into 5 Child-Mother pairs.
'''

cd /Users/abelvertesy/Google_Drive/X_react_Data/EXOME/perFamily/gVCF/raw_VCF

filename='gVCF_BGI_Leiden_exome_on_genome.raw.vcf'

cut -f 1-9,10-11 $filename > gVCF_BGI_Leiden.d1.vcf &
cut -f 1-9,12-13 $filename > gVCF_BGI_Leiden.d2.vcf &
cut -f 1-9,14-15 $filename > gVCF_BGI_Leiden.d3.vcf &
cut -f 1-9,16-17 $filename > gVCF_BGI_Leiden.d4.vcf &
cut -f 1-9,18-19 $filename > gVCF_BGI_Leiden.d5.vcf &

cut -f 7 gVCF_BGI_Leiden.d5.vcf > filters.d1.vec &

# select rows with a SNP: 'See Python script'
# tail -10000 gVCF_BGI_Leiden.d1.vcf > gVCF_BGI_Leiden.d1.play.vcf
# tail -20000 gVCF_BGI_Leiden_exome_on_genome.raw.vcf > gVCF_BGI_Leiden_exome_on_genome.play.vcf

echo "DONE"