#!/usr/bin/python -u
'''
### 08.x.Generate_cout_table.py
- call process_sam_cel_seq for each .allelic.sam file

### process_sam_cel_seq.pl
- usage: usage: -in=INPUTFILE.sam
	- s=1 if single end 0 is paired end -cmp=LIST_OF_READS.csv
	- s_flag= 1 or 0 ( if 1 then produce separate files for sense and antisense strand )
	- u= 1 or 0 ( if 1 then only map reads that map to only one strand, optional)
	- gff=TRANSCRIPTOME.gff
	- uniq=1 (optional, only unique reads)
- Check Output:
	-Nonempty s-, r- and cout files in the input (Allelic SAM) folder
'''

import os, glob
from subprocess import call

# Setup ----------------------
print 1111111
JOBname='cout_RS2_'
gene_set = "RefSeq"

options = "-s=1 -u=0 -uniq=1 -s_flag=0"
qsubme='qsub -q veryshort -cwd -V -pe threaded 1 -N '

# Prepare ----------------------
# AllelicSamDir ="/hpc/hub_oudenaarden/Abel/X_inact/SAM4/RefSeq/mat/"
AllelicSamDir ="/hpc/hub_oudenaarden/Abel/X_inact/SAM_ERCC/RefSeq/ERCC/"

LogDir = AllelicSamDir+"Logs/"
if not os.path.exists(LogDir):
	os.makedirs(LogDir)
# for dire in [RoutDir, SoutDir, LogDir]:
# 	print dire
# 	if not os.path.exists(dire):
# 		os.makedirs(dire)

print AllelicSamDir

os.chdir(AllelicSamDir)
ls_AllSam_fp =sorted(glob.glob("*.sam")); print ls_AllSam_fp

os.chdir(LogDir)
for i in range(0, len(ls_AllSam_fp)):
	AllSam_fp = AllelicSamDir+ls_AllSam_fp[i]
	JOBname_i = JOBname+str(i)
	qsub_call='echo "process_sam_strand.pl -in='+AllSam_fp, options, '" | ', qsubme, JOBname_i
	qsub_call=" ".join(qsub_call)
	print i, qsub_call
	call(qsub_call, shell=True)

print "Madonna 3.0"