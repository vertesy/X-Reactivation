#!/usr/bin/python -u
'''
### 02.x.Filter_split_SNP_tables.py

- Filter VCF file by rows, matching all of:
	- 0/1 genotype, except for d3[male] chrX and chrY where 1/1 is selected
	- FILTER = PASS
	- situated on an ordinary chromosome
	- no indel len(REF) <2 and len(ALT) <2 [this causes a minimal % [>0.7%] of false negatives, e.g. ALT "A,CGTCC" ]
	- SCRIPT DOES NOT FILTER FOR TRUE HETEROZYGOTES [statistical modeling of the sampling process]
		- Apparently no strong biases on the RNA level
	- No filtering on Depth
		- There are only  ~150 SNPs are under 10. Previously done VQSR takes carte of it?
- Output is used in 06.x.modify_genome.py
####!/usr/bin/env python -u
'''
import os, glob


# setup ------------------------------------------------
InputDir = "~/Google_Drive/X_react_Data/EXOME/perFamily/gVCF/raw_VCF/"


# colnames = "CHROM" ,"POS" ,"ID" ,"REF" ,"ALT" ,"QUAL" ,"FILTER" ,"INFO" ,"FORMAT","child" ,"mother"
os.chdir(InputDir)
OutDir = InputDir+"PASS/"
if not os.path.exists(OutDir):
	os.makedirs(OutDir)

for infname in glob.glob("*.d[0-9].vcf"):
	donor = infname.split(".")[1]
	print donor
	outfname='PASS/pass.kidhet.'+infname
	outfile=open(outfname, 'w+')

	print "current outfile is: " + outfname
	with open(infname) as f:
		for line in f:
			if line.startswith('#C'):
				outfile.write(line)
			elif line.startswith('c'):	# print line
				splitted = line.split("\t")
				CHROM = splitted[0]
				REF = splitted[3]
				ALT = splitted[4]
				FILTER = splitted[6]
				child_splitted = splitted[9].split(":")
				GT_kid = child_splitted[0]
				# if (',' in ALT) and (GT_kid == "0/1"):												# Partial solution for: ALT "T,GAT"-problem > it is a false negative with len(ALT)<2 filter.
				# 	ALT = ALT.split(',')[0]															# split by comma, and select 1st entry [GT_kid = 0/1!]
				# 	splitted[4] = ALT 																# put back clean genotype
				# 	line = "\t".join(splitted)														# re-parse the line
				# elif (',' in ALT) and (GT_kid == "0/2"):
				# 	ALT = ALT.split(',')[1]
				# 	splitted[4] = ALT
				# 	child_splitted[0] = "0/1"														# reset kid's GT
				# 	splitted[9] = child_splitted													# reset kid's GT
				# 	line = "\t".join(splitted)
				if GT_kid == "0/1" and FILTER == "PASS" and len(CHROM)<6 and len(REF) <2 and len(ALT) <2:	# print CHROM; print GT_kid ; print FILTER
					outfile.write(line)
				if donor == 'd3': 													# allow X SNPs for donor 3 # print ' Special rules for the male donor 3'
					if CHROM == "chrY" and GT_kid == "1/1" and FILTER == "PASS" and len(REF) <2 and len(ALT) <2:
						outfile.write(line)
					elif CHROM == "chrX" and GT_kid == "1/1" and FILTER == "PASS" and len(REF) <2 and len(ALT) <2:
							outfile.write(line)

print "DONE"