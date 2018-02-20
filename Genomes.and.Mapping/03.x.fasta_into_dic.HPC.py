#!/usr/bin/env python
'''
### 03.x.fasta_into_dic.HPC.py
- read a multifasta into a dictionary object
- write out the disctionary as pkl object
'''

import pickle, os

# setup ------------------------------------------------
infname = "hg19.fa"; print infname
InputDir = "/hpc/hub_oudenaarden/gene_models/human_gene_models/"
os.chdir(InputDir)

if os.path.isfile(infname+'.pkl'):
	raw_input("Output file exist, do you want to overwrite? press ENTER / CTRL +C")

# initialize variables
per_chr_dict = {}

# Read the sequence lines up to the blank line.
i=0
with open(infname) as f:
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if len(line) <9:			# get rid of weird chromosomes
				i += 1
				if i > 1:				# put elements in the per_chr_dict just after it encounters the 2nd chr
					sequence_of_chr = "".join(sequence_lines)
					per_chr_dict[key] = sequence_of_chr
				key =line.replace(">", "")
				print key
			sequence_lines = []
		if not line.startswith('>') and line !="":
				sequence_lines.append(line)
		# if line == "" or ENDOFFILE:
		if line == "":
			print "\n"
if line != "": 			# if the file did not contain a last empty line
	sequence_of_chr = "".join(sequence_lines)
	per_chr_dict[key] = sequence_of_chr

# print "      LAST LINE:    "+line

# Pickle per_chr_dict using protocol 0.
output = open(infname+'.pkl', 'wb')
pickle.dump(per_chr_dict, output)


print 'Genome is written out to: ', InputDir, infname,'.pkl'

# print per_chr_dict.values()
print per_chr_dict.keys()


print "DONE"