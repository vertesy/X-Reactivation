#!/usr/bin/python -u
'''
### 03.genome_fa_reformatter_1chr_per_line.py
- read a multifasta into a dictionary object
- write out .reformatted.fa:
	- with 1chr per line
	- weird chromosomes removed [by name string > 6, eg; ">chr6_cox_hap2"]
'''
import os

# Setup ------------------------------------------------------
print "\n 03.genome_fa_reformatter_1chr_per_line.py", 1111
InputDir = "~/Downloads/00_zacc/"
inGenomeName = "hg19.fa"

# Init ------------------------------------------------------
os.chdir(InputDir)
outfname = inGenomeName.rsplit('.', 1)[0] + ".reformat.fa"; # print outfname
if os.path.isfile(outfname):
	raw_input("Output file exist, do you want to overwrite? press ENTER / CTRL +C")


# Read the genome into a dic ------------------------------------------------------
per_chr_dict = {}										# initialize variables
chrom_count=0
with open(inGenomeName) as f:								# Read the sequence lines up to the blank line.
	for line in f:
		line = line.strip()
		if line.startswith('>'):						# HEADER / CHROM name > Key
			print line
			if len(line) >6:							# weird chromosome headers do not get written ">chr22" vs ">chr6_mcf_hap5"
				append_it = False
			elif len(line) <=6:							# normal chromosomes [by name]
				append_it = True
				chrom_count += 1
				if chrom_count >= 2:					# put elements in the per_chr_dict just after it encounters the 2nd Header
					per_chr_dict[key] = "".join(sequence_lines)
				key =line.replace(">", ""); print "   Key:", key
				sequence_lines = []
		elif not line.startswith('>') and (append_it == True) and line !="":	# append non-header lines
				sequence_lines.append(line)
per_chr_dict[key] = "".join(sequence_lines)				# Join the last chromosome
print "\n Genome importing is finished"

for key in sorted(per_chr_dict):						# Count chromosome's length [bases]
	print key, "seq. length: ", len (per_chr_dict[key])

# Write the genome out to reformatted .fa file ------------------------------------------------------
print "\n   Witing out the file: ", outfname
f=open(outfname,'w')
for key in sorted(per_chr_dict):
	f.write('>'+key+"\n")
	f.write(per_chr_dict[key]+'\n')
	print key
f.close()

print per_chr_dict.keys()

print "DONE Reformatted genome is written out to:", outfname
