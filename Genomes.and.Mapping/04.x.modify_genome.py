#!/usr/bin/python -u
'''
### 04.x.modify_genome.py

- Load genome into a dictionary [mm10.pkl]
- Read filtered VCF file line by line: ['pass.kidhet.mgp.v4.snps.dbSNP.vcf.columnFilt.head']
	- look up positions [0 indexing!]
	- check if match with REF and check upper- or lower-case.
	- ***check if kidGT == "0/1" and matGT == '1/1'***
	- replace with ALT converted to upper- or lower-case
	- write out to a reformat.fa [1 chr per line]
- Usage
	- replace outfname [pat or mat] in line 65
	- replace matGT == ['0/0' or 1/1] in line 95
	- replace glob (*.vcf) to restrict to a set of files in line 29
'''
print "04.x.modify_genome.py started Yuhuu"

import os, pickle, glob

# setup ------------------------------------------------
infname = 'hg19.fa.pkl'
# GenomeDir = "~/Dokumentumok/Tanulas/PhD/AvanO/00_Genomes/"
GenomeDir = "./genomes/pkl/"
InputDir = '~/Google_Drive/X_react_Data/EXOME/perFamily/gVCF/raw_VCF/PASS/'
os.chdir(InputDir)
VCF_Files = glob.glob("*.vcf"); print VCF_Files
# raw_input("Is the list of VCF Files correct? press ENTER / CTRL +C")

# open files ------------------------------------------------
infile = open(GenomeDir+infname, 'rb')
print 'Input genome .pkl: ', infile
per_chr_dict = pickle.load( infile )
print 'Chromosomes in Dictionary: ', per_chr_dict.keys()

# Count chromosome length a list
for key in sorted(per_chr_dict):
	print key, "seq. length: ", len(per_chr_dict[key])

# Create a list ------------------------------------------------
print "\n   Create a list from sequence of: "
for key in sorted(per_chr_dict):
	print key
	per_chr_dict[key] = list(per_chr_dict[key])

# per_chr_dict_original = copy.copy(per_chr_dict) # save original
per_chr_dict_original = dict(per_chr_dict.items())

# print (per_chr_dict_original['chr1'][1:10]) 	# the first ten bases
# print (per_chr_dict['chr1'][1:10])


# open vcf file ------------------------------------------------

for vcf in glob.glob("*.vcf"):
	print vcf
	donor = vcf.split(".")[3]

	# per_chr_dict = dict(per_chr_dict_original.items()) # overwrite with original [so SNPs are not inherited]
	per_chr_dict = {chrom: [base for base in seq] for chrom, seq in per_chr_dict_original.items()}
	print 'O-o-overwritten with original'

	total_lines = (sum(1 for line in open(vcf))-1)
	print 'VCF file: ', vcf, 'SNPs in file: ', total_lines

	# outfname = (infname.split(".")[0])+'.'+donor+'.mat.reformat.clean.fa'
	outfname = (infname.split(".")[0])+'.'+donor+'.HetHet.reformat.clean.fa'
	print 'Output genome: ', outfname

	old_CHROM = CHROM = ""
	with open(vcf) as f:
		i = HetHomSNPs = SNPs_changed = MultiALT = UnexpectedReferenceBase = 0
		for line in f:
			i += 1
			if not line.startswith('#'):
				splitted = line.split("\t")
				old_CHROM = CHROM
				CHROM = splitted[0]
				if old_CHROM != CHROM: 								# report the progress
					print CHROM
				POS = int(splitted[1])-1
				REF = splitted[3]
				if len(REF) > 1:
					print 'Deletion! len(REF) > 1'; break
				ALT = splitted[4]
				kidGT = (splitted[9].split(":")[0])
				matGT = (splitted[10].split(":")[0])
				# if len(REF) == 1 and len(ALT) == 1 and kidGT == "0/1" and matGT == '1/1':
				if len(REF) == 1 and len(ALT) == 1 and kidGT == "0/1" and matGT == '0/1':
					HetHomSNPs += 1
					toReplace = per_chr_dict[CHROM][POS]
					if toReplace.upper() == REF: # print "OK:   ", REF,  "in line: ", i  # print per_chr_dict.get(CHROM)[POS:POS+1]
						SNPs_changed += 1
						if toReplace.isupper():
							per_chr_dict[CHROM][POS] = ALT
						elif toReplace.islower():
							per_chr_dict[CHROM][POS] = ALT.lower() # print per_chr_dict.get(CHROM)[POS], ALT, "\n"
					else:
						UnexpectedReferenceBase += 1
						print "    Mismatch, in genome: ", per_chr_dict[CHROM][POS], " at: ", CHROM, POS, "in line: ", i, "REF in .vcf:", REF
				elif len(ALT) > 1:									# Redundant error catch, as operates on .PASS.vcf that fulfills all conditions in if() statement above.
					MultiALT += 1
					print '      UNEXPECTED LINE (Multiple ALT bases), ALT:', ALT, 'genotypes:', kidGT, matGT, 'line nr: ', i
	print "          SNPs_changed: ", SNPs_changed, " in ", outfname
	print "          UnexpectedReferenceBase: ", UnexpectedReferenceBase, " in ", outfname
	print "          MultiALT: ", MultiALT, " in ", outfname

	lf=open(GenomeDir + outfname + '.logs.tsv','w')
	lf.write("\nHetHomSNPs: \t" + str(HetHomSNPs))
	lf.write("\nSNPs_changed: \t" + str(SNPs_changed))
	lf.write("\nMultiALT: \t" + str(MultiALT))
	lf.write("\nUnexpectedReferenceBase: \t" + str(UnexpectedReferenceBase))
	lf.close()



	# Join a list for each chr ---------------------------------------------------------------
	print "\n   Joining sequence of: "
	for key in sorted(per_chr_dict):
		print key
		per_chr_dict[key] = "".join(per_chr_dict[key])

	print "\n   Writing out the file: ", outfname
	f=open(GenomeDir + outfname,'w')
	for key in sorted(per_chr_dict):
		f.write('>'+key+"\n")
		f.write(per_chr_dict[key]+'\n')
		print key
	f.close()

print 'DONE'