#!/usr/bin/python -u
'''
### 06.x.Invite_BWA.HPC.3.py
/hpc/hub_oudenaarden/Abel/bin/x_reactivation/DevCell_analysis/06.x.Invite_BWA.HPC.3.py

- Align each .fastq file towards N different genomes [fq extension incompatible!]
	- Genomes can be e.g.g: REF ALT PAT MAT kidHET
- The script parses and passes on the bash command for the UMC server
- INDEXING
	- generate .sai index file for each reaf
- ALIGNMENT
	- [BWA aln aligning parameters](http://bio-bwa.sourceforge.net/bwa.shtml) are:
		-q INT 	Parameter for read trimming. BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length. [0]
		-n NUM 	Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
		-k INT 	Maximum edit distance in the seed [2]
		-l INT 	Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for "-k 2". [inf]
		-t INT 	Number of threads (multi-threading mode) [1]
	- Parameters
		- bwa aln -q 30 -n 0.04 -k 2 -l 200 -t 1
- OUTPUT CONVERSION
	- generate pe /se .sam files
- Usage: Replace the following
	- FastqDir1, FastqDir2, GenomeDir, SamDir, SAIDir, LogDir
	- paired_end
	- SELECT 1: indexing, aligning, sam
'''

print 11
import os, glob
from subprocess import call

# Setup ------------------------------------------------------
indexing = 0
aligning = 0
sam = 1
paired_end = 0

JOBname = "RS_aln"
REFS = "RefSeq"

# FastqDir1 = "/hpc/hub_oudenaarden/Abel/X_inact/FASTQ2/"; 		# FastqDir1
FastqDir1 = "/hpc/hub_oudenaarden/Abel/X_inact/FASTQ/"; 		# FastqDir1
GenomeDir = "/hpc/hub_oudenaarden/Abel/genomes/X_React_used/"
SamDir = "/hpc/hub_oudenaarden/Abel/X_inact/SAM3/"+REFS+"/"
SAIDir = "/hpc/hub_oudenaarden/Abel/X_inact/SAI3/"+REFS+"/"


if indexing:
	LogDir = GenomeDir+"Logs"; print LogDir;
elif aligning:
	LogDir = SAIDir+"Logs"; print LogDir;
elif sam:
	LogDir = SamDir+"Logs"; print LogDir;
if not os.path.exists(LogDir):	os.makedirs(LogDir)

BWA_location = '/hpc/hub_oudenaarden/bin/software/bwa-0.6.2/bwa'
PARAMETERS='aln -q 30 -n 0.04 -k 2 -l 200 -t 4'
qsubme='qsub -q long -cwd -V -M a.vertesy@hubrecht.eu -m beas -pe threaded 4 -N '+JOBname

os.chdir(FastqDir1)
list_of_fastq = glob.glob("*.fastq"); print list_of_fastq

# if paired_end:
# 	os.chdir(FastqDir2)
# 	list_of_fastq2 = glob.glob("*.fastq"); print list_of_fastq2

list_of_fname_trunks =[i.split('.fastq')[0] for i in list_of_fastq]

print "\n", list_of_fname_trunks

os.chdir(GenomeDir)
list_of_genomes = glob.glob("*.fa"); print list_of_genomes

os.chdir(LogDir)
for GENOME in list_of_genomes:
	i= -1
	print "Aligned against: ", GENOME
	sex=GENOME.split(".")[2]; # print sex
	donor=GENOME.split(".")[1]; print donor
	if indexing: 	# reference indexing ------------------------------------------------------------------------------------------
		# qsub_call = 'echo "bwa index -a bwtsw '+GenomeDir+GENOME+'" | '+qsubme
		qsub_call ='echo "',BWA_location+' index -a bwtsw '+GenomeDir+GENOME+'" | '+qsubme
		qsub_call="".join(qsub_call)
		print qsub_call
		# call(qsub_call, shell=True)
	elif aligning: 	# Find SA coordinates for each read ------------------------------------------------------------------------------------------
		for FQ1 in list_of_fastq:
			i = i+1; 			# print i, FQ1
			if FQ1.split('_', 1)[0] == donor:
				# if paired_end:
				# 	qsub_call='echo "',BWA_location+' '+PARAMETERS+' '+GenomeDir+GENOME+' '+FastqDir1+FQ1+' '+FastqDir2+list_of_fastq2[i]+' > '+SAIDir+list_of_fname_trunks[i]+'.'+sex+'.sai" | ',qsubme
				if not paired_end:
					qsub_call='echo "',BWA_location+' '+PARAMETERS+' '+GenomeDir+GENOME+' '+FastqDir1+FQ1+' > '+SAIDir+list_of_fname_trunks[i]+'.'+sex+'.sai" | ',qsubme
				qsub_call="".join(qsub_call)
				print i, qsub_call
				# call(qsub_call, shell=True)
	elif sam: 	# Create .sam output from .sai ------------------------------------------------------------------------------------------
		if paired_end:
			PARAMETERS = ' sampe -n 100 -N 100 '
		elif not paired_end:
			PARAMETERS = ' samse -n 100 '
		for trunk in list_of_fname_trunks:
			i = i+1; 			# print i, FQ1
			if trunk.split('_', 1)[0] == donor:
				qsub_call='echo "',BWA_location+PARAMETERS+GenomeDir+GENOME+' '+SAIDir+trunk+'.'+sex+'.sai '+FastqDir1+trunk+'.fastq > '+SamDir+trunk+'.'+sex+'.sam" | ',qsubme
				qsub_call="".join(qsub_call)
				print trunk, i, qsub_call
				call(qsub_call, shell=True)

print "Walhalla 2.0"
