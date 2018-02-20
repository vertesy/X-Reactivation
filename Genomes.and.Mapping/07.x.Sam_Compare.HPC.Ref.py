#!/usr/bin/python -u
'''
### 07.x.Sam_Compare.HPC.Ref.py

- Split .sam files into allele specific outputs according to the smaller edit distance
- Read the 2 differentially mapped input sam files and compare line-by-line
- Discard reads that
	- are not mapped in both genomes [ alignment FLAGs each either ==16 or ==0]
	- have more than one *optimal* alignment [optimal by X0 tag = 1, allowing more suboptimal alignments]
- Line-by-line compare edit distances [by NM tag]
- Output the better mapping read into a
	- .mapped.Ref1.sam OR
	- .mapped.Ref2.sam file.
- If mapped equally good, output read from ref_name1.sam to:
	- .mapped.uninformative.sam
- Usage
	- Provide: SamDir, ref_names [also a subfolder under samdir]
'''

import os, sys, glob, itertools


# Setup ----------------------
# ref_name1 = "mat"
# ref_name2 = "pat"
ref_name1 = "HetHet"
ref_name2 = "Ref"

SamDir = "./X_inact/SAM_ERCC/RefSeq/Ref/"

# Prepare----------------------
OutDir = SamDir+"Allelic_HetHet_Ref_SAMs/";
LogDir = OutDir+"Logs/"; print OutDir, LogDir
if not os.path.exists(OutDir):	os.makedirs(OutDir)
if not os.path.exists(LogDir):	os.makedirs(LogDir)


print (SamDir+ref_name1)
os.chdir(SamDir+ref_name1+'/')
ls_Sam1 = sorted(glob.glob("*.sam"));	print ls_Sam1
os.chdir(SamDir+ref_name2+'/')
ls_Sam2 = sorted(glob.glob("*.sam"));	print ls_Sam2

ls_cellIDs = [i.split('.')[0] for i in ls_Sam1]; print ls_cellIDs

cmdargs = list(sys.argv)
i = int(cmdargs[1])-1; print i
cellID = ls_cellIDs[i]

SAM1_fnp = SamDir+ref_name1+'/'+ls_Sam1[i]; print SAM1_fnp
SAM2_fnp = SamDir+ref_name2+'/'+ls_Sam2[i]; print SAM2_fnp

outSAM1_fnp = OutDir+ls_cellIDs[i]+'.mapped.'+ref_name1+'.allelic_exp.sam'; print outSAM1_fnp
outSAM2_fnp = OutDir+ls_cellIDs[i]+'.mapped.'+ref_name2+'.allelic_exp.sam'; print outSAM2_fnp
outSAM3_fnp = OutDir+ls_cellIDs[i]+'.mapped.uninformative'+'.allelic_exp.sam'; print outSAM3_fnp

if any([os.path.isfile(j) for j in outSAM1_fnp, outSAM2_fnp, outSAM3_fnp]):		# raw_input("Output file exist, do you want to overwrite? press ENTER / CTRL +C")
	os.remove(outSAM1_fnp); os.remove(outSAM2_fnp); os.remove(outSAM3_fnp)

reads = readsA = readsB = readsU = 0
samFiles = [open(f) for f in [SAM1_fnp, SAM2_fnp]]

os.chdir(LogDir)
print "    ", cellID
for line_tuple in itertools.izip(*samFiles):
	A_line = line_tuple[0]
	B_line = line_tuple[1]
	# print A_line, B_line
	if A_line.startswith("@") and B_line.startswith("@"):
		with open(outSAM1_fnp, 'a') as OutSamFile1:
			OutSamFile1.write(A_line)
		with open(outSAM2_fnp, 'a') as OutSamFile2:
			OutSamFile2.write(B_line)
		with open(outSAM3_fnp, 'a') as OutSamFile3:
			OutSamFile3.write(A_line)
	elif (A_line.startswith("@")) ^ (B_line.startswith("@")):	# if one of the samFiles end the header earlier
		print "ERROR: Header lenght is unequal"; break
	if not A_line.startswith("@"):
		reads += 1
		if not (reads % 10000):	print reads
		FLAG_A=A_line.split('\t')[1]							# 2nd field is the ALIGNMENT_FLAG:
		FLAG_B=B_line.split('\t')[1]
		if (FLAG_A in ['0','16']) and (FLAG_B in ['0','16']): # if not unmapped
			ReadID_A = A_line.split('\t')[0]					# 1st field is the read-ID
			ReadID_B = B_line.split('\t')[0]
			NM_tag_A = (A_line.split('\t')[12]).split(":")[2]	# 13rd field is the NM_tag, encoding the edit distance.
			NM_tag_B = (B_line.split('\t')[12]).split(":")[2]	# tidy up to a single integer in the 3rd field of "NM:i:0"	# NM_tag_B = NM_tag_B.split(":")[2]
			X0_tag_A = A_line.split('\t')[13]					# X0_tag reports the "Number of best hits"
			X0_tag_B = B_line.split('\t')[13]
			if (X0_tag_A == "X0:i:1") and (X0_tag_B == "X0:i:1") and (ReadID_A == ReadID_B):	# if mapped uniquely
				if NM_tag_A < NM_tag_B:							# if maps better to ref_A >>> .ref_A.sam
					readsA += 1;
					with open(outSAM1_fnp, 'a') as OutSamFile1:
						OutSamFile1.write(A_line)
				elif NM_tag_A > NM_tag_B:							# if maps better to ref_B >> .ref_A.sam
					readsB += 1;
					with open(outSAM2_fnp, 'a') as OutSamFile2:
						OutSamFile2.write(B_line)
				elif NM_tag_A == NM_tag_B:						# if equally maps to both references >> .uninformative.sam
					readsU += 1
					with open(outSAM3_fnp, 'a') as OutSamFile3:
						OutSamFile3.write(A_line)					# it does not matter if A_line or B_line is written out
			elif ReadID_A != ReadID_B:
				print "          NON MATCHING LINES: DIFFERENT READ ID IN THE SAME LINE of the 2 samFiles", reads
for f in samFiles:
	f.close()

total_uq_mapped = (readsA+readsB+readsU)
print reads, "reads processed,", readsA, "are from Ref_A", readsB, "are from Ref_B, and in total", total_uq_mapped,"or", 100*total_uq_mapped/reads, "% were uniquely mapped [optimally and suboptimally]\n"
print "DONE"

## Check output by
# grep "DONE" * | wc -l;			# Nr of files finished in order
# grep -L "DONE" *.o* | wc -l;		# Nr of files not finished on time
# grep -L "DONE" *.o* ;				# The files not finished on time

