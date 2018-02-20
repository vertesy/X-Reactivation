#!/usr/bin/python -u
'''
### 01.xb.DNA_SNP_Statistics_md_parser.py
- Calculate summary statistics on per-family DNA SNPs
- Parse the output into a markdown table
'''

import os, glob


os.chdir("~/Google_Drive/X_react_Data/EXOME/perFamily/gVCF/raw_VCF/PASS")
list_of_PASS_VCF = glob.glob("*.vcf"); print list_of_PASS_VCF


print '| Donor | Heterozygote SNPs | Het-Het SNPs | Heterozygote xSNPs | Het-Het xSNPs |'
print '|---|---|---|---|---|'; # GitHub requires 3 dashes!!!
d=0
for vcf in list_of_PASS_VCF:
	d+=1
	with open(vcf) as f:
		i=0; x=0; hethet =0; hethetX =0
		for line in f:
			i+=1
			chrom=line.split("\t")[0]
			if chrom == "chrX":
				x+=1
			matGT=line.split("\t")[10]
			if matGT.startswith('0/0') or matGT.startswith('1/1'):
				hethet+=1
				chrom=line.split("\t")[0]
				if chrom == "chrX":
					hethetX+=1


	print '| d'+str(d),'|',i, "\t| " ,hethet, "\t| ", x, "\t| ", hethetX, "|"