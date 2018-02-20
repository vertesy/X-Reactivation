#!/usr/bin/env bash
'''
### 05.x.Reference_generation_from_genome_pipeline.sh
- convert UCSC to GTF by *01.ucsc2gtf.pl*
	- http://genome.ucsc.edu/cgi-bin/hgTables
- Merge of all exons in all isofrom by *02.merge_isoforms_gtf.pl*
- Convert the genome from normal [multifasta] genome.fa > genome.refomat.fa with: 1 chr / 1 line
- Extract sequences from genome.fa by *03.gtf2fa.pl*
	- Repeat for each genome
- Add ERCC SpikeIn sequences
- Replace poly-A stretches by poly N-s by *04.mask_polyA.pl*
'''


cd '/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg19/RefSeq/'

scripts_path='/Users/abelvertesy/x_reactivation/analysis/DevCell_analysis/05.5.Reference_generation/'

# Transform ucsc it into gtf format
$scripts_path'01.ucsc2gtf.pl' -in='01.hg19.RefSeq.ucsc' -out='02.hg19.Refseq.gtf' -m='02.hg19.Refseq.tsv'
'maybe new line after last'
# gtf 1 line per exon, multiple isoforms for the same gene follow each other: each marked by 'transcript' in col3

#  Merge of all exons in all isofrom
$scripts_path'02.merge_isoforms_gtf.pl' -in='02.hg19.Refseq.gtf' -cl='02.hg19.Refseq.tsv' -out='03.hg19.Refseq_genes.gtf'
# -cl= key value: gene name >> all isoforms
# -out= no more different isoforms




#  genome reformatter : 1line/ chr
cd '/Users/abelvertesy/Downloads'
more hg19.fa |perl -ane 'chomp; $x = $_; $x = "\n".$_."\n" if /^\>/; print $x;' >jk
more jk |perl -ane '$i++; print if $i > 1;' > jk2
echo 1|perl -ane 'print "\n";' > jk2
cat jk jk2 >  'hg19.reformat.fa'
rm jk
rm jk2


# extract sequences from genome.fa
### genome has to be in 1 line / chr format
$scripts_path'03.gtf2fa.pl' -in='03.hg19.Refseq_genes.gtf' -ref='/Users/abelvertesy/Downloads/hg19.reformat.fa' > 'hg19.REF.Refseq_genes.fa'


# Add ERCC SpikeIn sequences
ERCC_fnp='/Users/abelvertesy/Downloads/ERCC/ERCC92_reformat.fa'
cat 'hg19.REF.Refseq_genes.fa' $ERCC_fnp > 'hg19.REF.Refseq_genes_ERCC.fa'


# Replace poly-A stretches by poly N-s [BWA counts as mismatch]
$scripts_path'04.mask_polyA.pl' -in='hg19.REF.Refseq_genes_ERCC.fa' -n=15 -out='hg19.REF.Refseq_genes_ERCC_polyN.fa'

# 'Make an index with bwa in the next script'



# 'Reference GENCODE -----------'
# $scripts_path'03.gtf2fa.pl' -in='/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg19/Gencode/03.hg19_GencodeV19_genes.gtf' -ref='/Users/abelvertesy/Downloads/00_zacc/hg19_reformat.fa' > '/Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg19/Gencode/hg19.REF.GencodeV19.fa'
# cd /Users/abelvertesy/Dokumentumok/Tanulas/PhD/AvanO/Data_analysis/UCSC_GTF_tables/hg19/Gencode/
# cat 'hg19.REF.GencodeV19.fa' $ERCC_fnp > 'hg19.REF.GencodeV19_genes_ERCC.fa'
# $scripts_path'04.mask_polyA.pl' -in='hg19.REF.GencodeV19_genes_ERCC.fa' -n=15 -out='hg19.REF.GencodeV19_genes_ERCC_polyN.fa'