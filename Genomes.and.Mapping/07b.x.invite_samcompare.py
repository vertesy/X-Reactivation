#!/usr/bin/python -u
'''
### 07b.x.invite_samcompare.py
- Usage: change:
	- JOBname
	- LogDir
	- -q "queue" in qsubme
	- -pe threaded "4" in qsubme
- Change:
	- NrCells & Script: 07.x.Sam_Compare.HPC.Ref.py
'''

from subprocess import call
import os

NrCells = 63

JOBname = 'sCmp_RS'
LogDir = "/RefSeq/Allelic_SAMs/Logs"

if not os.path.exists(LogDir):	os.makedirs(LogDir)
os.chdir(LogDir)
qsubme='qsub -q short -cwd -V -m beas -pe threaded 4 -N '

for i in range(1,NrCells+1):
	# qsub_call='echo "', "07.x.Sam_Compare.HPC.Ref.py ", str(i),   '" | ', qsubme, JOBname, "_", str(i)
	qsub_call='echo "', "07.x.Sam_Compare.HPC.py ", str(i),   '" | ', qsubme, JOBname, "_", str(i)
	qsub_call="".join(qsub_call)
	print i, qsub_call
	call(qsub_call, shell=True)

print "Mari Carmen"

# os.chdir(SamDir)
# os.system("mv *.pat.sam pat/")
# os.system("mv *.mat.sam mat/")
