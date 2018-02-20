#!/bin/bash
'''
### 09.x.Global_merger.sh
- Copy to separate folders
	- .cout.csv to X_inact/COUNT
	- .rout.csv to /COUNT/ROUT
	- .sout to /COUNT/SOUT
- Check if all files contain the same number of genes
	- Less genes mean unfinished previous script
- Merge all "finished" *.cout.csv
	- add rownames and header
- Usage: replace
	- gene_set
	- expected_nr_entries
'''

# Setup ------------------------------------------------------------------
UseUninf=0; # 0 if ususal pat mat mapping!!!
gene_set="RefSeq"
# expected_nr_entries=26362
expected_nr_entries=92

ref_name1="mat"
ref_name2="pat"
# ref_name1="Ref"
# ref_name2="HetHet"
merged="global."$gene_set".cout.tsv"
echo $merged

cd
. .bash_profile


# Prepare ------------------------------------------------------------------
# AllelicSamDir="/hpc/hub_oudenaarden/Abel/X_inact/SAM4/RefSeq/mat/"
AllelicSamDir="/hpc/hub_oudenaarden/Abel/X_inact/SAM_ERCC/RefSeq/ERCC/"

# CoutDir="/hpc/hub_oudenaarden/Abel/X_inact/COUNT5/"$gene_set"/"; mkdir $CoutDir
# CoutDir="/hpc/hub_oudenaarden/Abel/X_inact/COUNT5_total/"; mkdir $CoutDir
CoutDir="/hpc/hub_oudenaarden/Abel/X_inact/COUNT_ERCC/"; mkdir $CoutDir
# CoutDir="/hpc/hub_oudenaarden/Abel/X_inact/COUNT2/"$gene_set"/Uninformative/"
SoutDir=$CoutDir"SOUT/"; mkdir $SoutDir
UninfDir=$CoutDir"UninfDir/"; mkdir $UninfDir

cd $AllelicSamDir

mv *.cout.csv $CoutDir
rm -rf *.rout.csv
mv *.sout $SoutDir

cd $CoutDir
mv  *.uninformative.allelic_exp.cout.csv $UninfDir

echo 111

if [ $UseUninf -eq 1 ] ; then
	cd $UninfDir
fi


FILES=*.cout.csv
for f in $FILES
do
	mv $f ${f%%mapped*}${f##*mapped.};		# Rename files
done


touch $merged
FILES=*.cout.csv
for f in $FILES
do
	echo "Processing $f file..."
	echo "             "$(wc -l $f)
	if [ $(wc -l $f | cut -d ' ' -f 1) -eq $expected_nr_entries ]
	then
		echo ${f%%.allelic*} > $f'.read_col';			# Sample ID
		cut $f -f 2 >> $f'.read_col';					# read counts
		echo "OK"
	else
		echo $(wc -l $f | cut -d ' ' -f 1) "ENTRIES IN UNFINISHED FILE" $f
			# echo ${f%%.allelic*} > $f'.read_col';			# Sample ID
			# cut $f -f 2 >> $f'.read_col';					# read counts
			# "This is a quick fix above!!!"
	fi

done
ls *.read_col

# add rownames
echo "\t" >  'RowNames.vec'
cut $f -f 1 >>  'RowNames.vec'

# # merge them
paste RowNames.vec *'.read_col' > $merged

head $merged

echo \n\n"Eurekaaaaa"
