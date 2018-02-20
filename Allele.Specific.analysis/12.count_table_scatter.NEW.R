######################################################################################################
# Count table analysis for Abels Tables
######################################################################################################
# source("~/analysis/DevCell_analysis/12.count_table_scatter.R")
rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# '''
# ### 12.count_table_scatter.R
# - Read in merged read count table, with alternating pat / mat columns.
# - Filter genes not present in the dataset in sufficient count
# 	- Thr: never above 5
# - ...
# '''

# Functions ------------------------
try (source ('~/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F)
library(stringr)

# Setup ------------------------
InputDir = "~/Google_Drive/X_react_Data/Abelz"
infile = '~/Google_Drive/X_react_Data/Abelz/zz.trimmed.RefSeq.cout.processing/AllCells.tsv'
setwd(InputDir)

setup_MarkdownReports(OutDir = "~/Google_Drive/X_react_Data/Abelz/X_react_2016_PGC_clustering", scriptname = "12.count_table_scatter.NEW.R", append = F)
OutDirOrig = OutDir

# Parameters ------------------------
SecondPart = F
nr_cells = 108
nr_cols = nr_cells*2

# minReadCount = 5000		# Minimum nr. of Reads that a cell has to have
# minSNP = 400			# Minimum nr. of Expressed Genes with SNP
# DP_thr = 25 			# Used only to demonstrate that higher min RC only makes good and bad cells less distighuisable.

# Meta data ------------------------
# CellTypeColz = read.delim("~/Google_Drive/X_react_Data/RNA/Expression_analysis/TranscriptomeCountTables_Sw_Hubr/ColorTheme/ColorThemeSimple.tsv", row.names=1)
# a = as.character(unlist(CellTypeColz)); names(a) =rownames(CellTypeColz); CellTypeColz = a
Impr_MAT = sort(read.simple.vec("~/Dropbox/X_reactivation/Metadata_d/Imprinted_genes/HUGO_updated_and_ICR_updated/Imprinted_confirmed_union_MAT_HUGO.vec"))
Impr_PAT = sort(read.simple.vec("~/Dropbox/X_reactivation/Metadata_d/Imprinted_genes/HUGO_updated_and_ICR_updated/Imprinted_confirmed_union_PAT_HUGO.vec"))
IC_members_iPat = sort(read.simple.vec ("~/Dropbox/X_reactivation/Metadata_d/Imprinted_genes/ICs/ICR_members_iPat.vec"))
IC_members_iMat = sort(read.simple.vec ("~/Dropbox/X_reactivation/Metadata_d/Imprinted_genes/ICs/ICR_members_iMat.vec"))
	IC_members_iPat = intersect(IC_members_iPat,Impr_PAT)
	IC_members_iMat = intersect(IC_members_iMat,Impr_MAT)

escapees = read.simple.vec ('~/Dropbox/X_reactivation/Metadata_d/X_escapees/X_Escape_intersect_always_esc.vec')
escapees_het = read.simple.vec ('~/Dropbox/X_reactivation/Metadata_d/X_escapees/X_Escape_intersect_HeterogenousData.vec')

# Manual annotation file Sample_Name	Sample_Name_simple	Sample_Name_old	Cell_type_annot	Cell_type_annot_simple ------------
# metadata_manual_annot = read.simple.tsv("~/Google_Drive/X_react_Data/Abelz/Metadata2/Manual_Metadata/metadata_manual_annot_AllCells.RSEM_FPKM.tsv")
metadata_manual_annot = read.simple.tsv("~/Google_Drive/X_react_Data/Abelz/Metadata2/Manual_Metadata/metadata_manual_annot_AllCells.PGC_clustering.tsv")
attach_w_rownames(metadata_manual_annot);
cell_IDs = rownames(metadata_manual_annot)

# Colors ---
# CellTypeColz = read.simple.tsv.named.vector("~/Google_Drive/X_react_Data/Abelz/Metadata2/ColorTheme/ColorTheme_old/ColorTheme.tsv")
CellTypeColz = read.simple.tsv.named.vector("~/Google_Drive/X_react_Data/Abelz/Metadata2/ColorTheme/ColorThemeNew/ColorTheme.final.tsv")
cell_type_col		 = CellTypeColz[cell_type_annot[cell_IDs]]; names (cell_type_col) = cell_IDs

# CellTypeColz_Simple = read.simple.tsv.named.vector("~/Google_Drive/X_react_Data/Abelz/Metadata2/ColorTheme/ColorTheme_old/ColorThemeSimple.tsv")
CellTypeColz_Simple = read.simple.tsv.named.vector("~/Google_Drive/X_react_Data/Abelz/Metadata2/ColorTheme/ColorThemeNew/ColorThemeSimple.final.tsv")
cell_type_col_simple = CellTypeColz_Simple[cell_type_annot_simple[cell_IDs]]; names (cell_type_col_simple) = cell_IDs

CellTypeColz_Simplest = CellTypeColz_Simple[c("SOMATIC", "EGC", "L+MGC")]
cell_type_col_simplest = CellTypeColz_Simplest[cell_type_annot_simplest[cell_IDs]]; names (cell_type_col_simplest) = cell_IDs

names(CellTypeColz_Simplest)[3] = 'L+MGC'
cell_type_annot_simplest = cell_type_annot_simple
cell_type_annot_simplest[cell_type_annot_simplest ==  "MGC"] = "L+MGC"
cell_type_annot_simplest[cell_type_annot_simplest ==  "LGC"] = "L+MGC"

# Read in Data ------------------------
df_ReadCount = read.table(infile, row.names=1, header =T); dim(df_ReadCount); # dimnames(df_ReadCount)

# Gene_Expr_in_Germ_Cells = na.omit(unlist(df_ReadCount))
# Gene_Expr_in_Germ_Cells=Gene_Expr_in_Germ_Cells[Gene_Expr_in_Germ_Cells>0 & Gene_Expr_in_Germ_Cells<1000]
# whist (Gene_Expr_in_Germ_Cells, breaks =100)

# - Filter genes not present in the dataset in sufficient count -----------
GenesCovered = (rowSums (df_ReadCount)>0);
llprint(sum(GenesCovered), "transcripts out of", dim(df_ReadCount)[1], "annotated transcripts are detected. (", percentage_formatter(sum(GenesCovered)/dim(df_ReadCount)[1]), ")")

df_ReadCount = df_ReadCount[GenesCovered,]
dim(df_ReadCount); nr_rows = dim(df_ReadCount)[1]; nr_cols == dim(df_ReadCount)[2]

# sort cell alphanumerically
index = sort(colnames(df_ReadCount))
df_ReadCount = df_ReadCount[ ,index]

# - Rename rows -----------
name_split = str_split_fixed (dimnames(df_ReadCount)[[1]],"__",2)
stopifnot ( sum(duplicated (name_split[,1]))==0 ); # There are non-unique rownames!
gene_names  = name_split[,1]
rownames (df_ReadCount) = gene_names
stopif(l(c(intersect(setdiff(IC_members_iPat,Impr_PAT), gene_names), intersect(setdiff(IC_members_iMat,Impr_MAT), gene_names))), message = "There are IC genes that are not part of the iPat nor iMat vectors, but are found in our dataset")

chrz = name_split[,2]; names (chrz) = gene_names
# Determine	SNP_type gene_category ---------------------------
gene_category = rep('autosomal',nr_rows); names (gene_category) = gene_names

llprint ("## Check how many of annotated gene categories we see in our dataset")
llogit ("#### Chr X")
gene_category[chrz =='chrX'] = 'X'
hits = lookup (escapees_het,rownames (df_ReadCount))
gene_category[hits$hit_poz] = 'escapees_het'
hits = lookup (escapees,rownames (df_ReadCount))
gene_category[hits$hit_poz] = 'escapees'
llogit ("#### Imprinted genes")
hits = lookup (Impr_MAT,rownames (df_ReadCount))
gene_category[hits$hit_poz] = 'Impr_MAT'
hits = lookup (Impr_PAT,rownames (df_ReadCount))
gene_category[hits$hit_poz] = 'Impr_PAT'

GeneAnnotationStatistic = table (gene_category);GeneAnnotationStatistic
MarkDown_Table_writer_NamedVector (GeneAnnotationStatistic)
xist_here = which(rownames(df_ReadCount) =='XIST')

gene_category_used = gene_category

# General Sample statistics ----------------------------------------------------------------------
llprint("## General Sample statistics")
NR_of_reads_per_cell = colSums(df_ReadCount);
wbarplot(NR_of_reads_per_cell, col=c(2,4), mdlink = T)
llprint ("There is still a clear difference between batch 1 and 2 in terms of raw read counts.
			 It is also visible on a histogram, where the 2nd peak consists the cells of the 1st batch:")
whist(NR_of_reads_per_cell, breaks=20 )

llprint("## Chromosome-X statistics")
NR_X_Reads_per_cell = colSums(df_ReadCount[gene_category =='X',])
wbarplot(NR_X_Reads_per_cell, col=c(2,4), mdlink = T)
llprint ("The dip in the middle is because you have a sole copy of chr X in D3, the male individual.
			 Otherwise the same trend is observed as above.")

# Sample annotation vectors -----------------------------------------------------------------------------
mat_colz = seq(1, nr_cols, 2)
pat_colz = seq(2, nr_cols, 2)
cell_IDs = substr(colnames(df_ReadCount)[pat_colz],1,8)

cell_NR = substr (cell_IDs,7,7); names(cell_NR) = cell_IDs
celltype = substr (cell_IDs,4,5); names(celltype) = cell_IDs
donor =	paste("D",substr (cell_IDs,2,2),sep=""); names(donor) = cell_IDs
# trimester =	as.numeric(donor == 'D3' | donor == 'D5') +1; names(trimester) =	cell_IDs
gender = as.numeric(donor != 'D3'); names(gender) = cell_IDs
# batch = replicate (nr_cells,2); batch [which(cell_NR[1:34] < 27)] = 1; names(batch) = cell_IDs ; # all batch 1 cells have NR<27 except d3_27_m
# batch[34]=1; 'Fix d3_27_am2'; sum(batch==1) # control: there should be 24
weeks_5 = c(D1 ='11.1', D2 ='10', D3 ='20', D4 ='12', D5 ='16.4'); weeks_5
weeks = paste (weeks_5[donor],"wk", sep=""); names(weeks) = cell_IDs

# Control =   (cell_type_annot == "GO" | cell_type_annot == "AD"  | cell_type_annot == "MALE" )
# Control =   (cell_type_annot_simplest == "LQ" | cell_type_annot_simplest == "SOMATIC")
Control =   (cell_type_annot_simplest == "SOMATIC")


OneLetterAnnot = substr(cell_type_annot,1,1); OneLetterAnnot[cell_type_annot == "MALE"] ="Y"
OneLetterDonorID = substr(donor,2,2); OneLetterDonorID
"Write out metadata later this script with HQ_samples"

# Split main table into smaller ones per mapping reference 'd1_01_fg_a	d1_01_fg_m	d1_01_fg_p	d1_01_fg_r' ----------
DP_mat = df_ReadCount[,mat_colz];  colnames (DP_mat) = cell_IDs;  colnames(DP_mat) =  cell_IDs
DP_pat = df_ReadCount[,pat_colz];  colnames (DP_pat) = cell_IDs;  colnames(DP_pat) =  cell_IDs;
DP_SUM_parental =  (DP_mat + DP_pat);  colnames(DP_SUM_parental) =  cell_IDs
write.simple.tsv (DP_pat); write.simple.tsv (DP_mat); write.simple.tsv (DP_SUM_parental)

# Correct for Xist ------------------------------------------
tmp = DP_pat[xist_here,]
DP_pat[xist_here,] = DP_mat[xist_here,]
DP_mat[xist_here,] = tmp
gene_category[xist_here] = "X"

PP_matrix = DP_pat / DP_SUM_parental
PP = cbind(PP_matrix, chrz); # write.simple.tsv (PP)
colnames(PP_matrix) = cell_IDs
write.simple.tsv(PP_matrix)


# Quality Filtering ------------------------------------------------------------------------------
llprint ("## Filtering HQ cells")
# llprint("Filtering is based on the number of reads useful for this analysis [SNP covering reads] and the total number of SNPs seen in a cell")
# FilterCriteria = c ( "minReadCount" = minReadCount, "minSNP" = minSNP )
# MarkDown_Table_writer_NamedVector (FilterCriteria)

# Number of SNPs containing_genes ------------------------------------------------------------------------
Sum_depth = (DP_SUM_parental);
NR_of_SNP_Containing_Genes = colSums(Sum_depth > 0)
NR_of_SNP_Containing_Genes_log2 = sort(log2(NR_of_SNP_Containing_Genes), decreasing = T)

# Number of reads per cell ------------------------------------------------------------------------
NR_of_reads_per_cell = colSums(DP_SUM_parental) # THIS VARIABLE IS REDIFNED HERE!!!
NR_of_reads_per_cell_log2 = log2(sort(colSums(DP_SUM_parental), decreasing = T))

llprint("We remove the lower quintile of the data that is heterogenously low quality.")
DeriveThresholds = T
if (DeriveThresholds) {
	minSNP = ceiling(2^percentile2value(NR_of_SNP_Containing_Genes_log2, percentile = .2)) # We remove the lower quintile, see below why
	minReadCount= ceiling(2^percentile2value(NR_of_reads_per_cell_log2, percentile = .2))
}

# Rationale for the removing the lower quintile -------------------------------------------------
llprint("### Removing the lowest quintile of the data is optimal to have the maximal number samples with homogenous quality.")
Cumulative_CV_of_ReadCounts = Cumulative_CV_of_gene_count_w_SNPs = NULL
for (i in 1:l(NR_of_SNP_Containing_Genes_log2)) {
	Cumulative_CV_of_gene_count_w_SNPs[i] = cv(NR_of_SNP_Containing_Genes_log2[1:i])
	Cumulative_CV_of_ReadCounts[i] = cv(NR_of_reads_per_cell_log2[1:i])
}
Cumulative_CV_of_ReadCounts[is.na(Cumulative_CV_of_ReadCounts)] = 0
Cumulative_CV_of_gene_count_w_SNPs[is.na(Cumulative_CV_of_gene_count_w_SNPs)] = 0
names(Cumulative_CV_of_ReadCounts) = names (NR_of_reads_per_cell_log2)
names(Cumulative_CV_of_gene_count_w_SNPs) = names(NR_of_SNP_Containing_Genes_log2)

quantiles = c(tercile = 0.333, quartile = 0.25, quintile = 0.2, sextile = 0.166, septile = 0.143)
boundaries = floor(nr_cells*(1-quantiles));boundaries

n=3
ccc = c(rep(3,boundaries[n]), rep(2, (nr_cells-boundaries[n])))

llprint("Cumulative CV shows, how the vairance of the data increases by including the next, lower quality sample:")
wbarplot(Cumulative_CV_of_gene_count_w_SNPs, vline =boundaries, col=ccc, lcol = c(1,1,2,1,1), lty = c(3,3,1,3,3), lwd = c(.75, .75, 1, .75, .75), mdlink = T)

wbarplot(Cumulative_CV_of_ReadCounts, vline = boundaries, col = ccc, lcol = c(1,1,2,1,1), lty = c(3,3,1,3,3), lwd = c(.75, .75, 1, .75, .75), mdlink = T)


#  Number of SNPs seen >0 ------------------------------
llprint ("### Filter 1: min SNP count")
wbarplot (NR_of_SNP_Containing_Genes_log2, sub = paste0("We exclude samples having less than ",minSNP," SNPs"), hline = log2(minSNP), mdlink = T)
llprint ("We eliminate ", sum(NR_of_SNP_Containing_Genes<minSNP), "cells having less than", minSNP, "SNP containing genes [expression != 0]: ",
			 names(which(NR_of_SNP_Containing_Genes<minSNP)))

NR_X_SNPs_per_Sample =		colSums ( 	DP_SUM_parental[chrz=="chrX",] > 0 , na.rm=T	)
NR_X_proper_SNPs_per_Sample = colSums ( 	DP_SUM_parental[gene_category=="X",] > 0, na.rm=T	)
NR_X_escapee_per_Sample =	colSums ( 	DP_SUM_parental[(gene_category=="escapees" | gene_category=="escapees_het"),] > 0 , na.rm=T	)
NR_auto_per_Sample =		colSums ( 	DP_SUM_parental[(gene_category=="autosomal"),] > 0 , na.rm=T	)
NR_Impr_SNPs_per_Sample = 	colSums ( 	DP_SUM_parental[gene_category=="Impr_PAT"  | gene_category=="Impr_MAT",] !=0 , na.rm=T	)

#  Number of SNPs seen above DP_thr ------------------------------
# llprint ("Looking at genes above a minimal, eg.",DP_thr,"reads gives less distinctive power for the low qual end, therefore we go for genes at any expression level. See supfig in folder.")
# Sum_depth[Sum_depth < DP_thr] = NA
# NR_SNP_genes_above_thr_log2 = log2(sort(colSums(Sum_depth != 0, na.rm =T), decreasing = F))
# NR_SNP_genes_above_thr_log2[NR_SNP_genes_above_thr_log2<0] =0
# wbarplot (NR_SNP_genes_above_thr_log2, mdlink = F)
# abline (h=log2(minSNP))

llprint ("### Filter 2: min SNP containing read count")
wbarplot(NR_of_reads_per_cell_log2, hline = log2(minReadCount), sub = paste0("We exclude samples having less than ",minReadCount," reads"), mdlink = T)

Filters= cbind(
	"SNPs" = NR_of_SNP_Containing_Genes_log2[cell_IDs],
	"Reads" = NR_of_reads_per_cell_log2[cell_IDs]
)
HQ_by_Transcriptome = as.named.vector(metadata_manual_annot[, "HQ_by_Transcriptome", drop=F])
HQ_samples = ( HQ_by_Transcriptome &  (NR_of_reads_per_cell[cell_IDs] > minReadCount)  & (NR_of_SNP_Containing_Genes[cell_IDs] > minSNP)); names(HQ_samples) = cell_IDs
nr_HQ_cells = sum(HQ_samples)
write.simple.vec(input_vec = HQ_samples, ManualName = "~/Google_Drive/X_react_Data/Abelz/Metadata2/HQ_samples.vec")
HQ_Cellnames = names(which(HQ_samples))
DiscardedCells = names(which(!HQ_samples))
llprint ("After both steps, we discard", l(DiscardedCells),"cells from the dataset:", DiscardedCells)

NR_of_SNP_Containing_Genes["d1_31_af"] > as.numeric(minSNP)+1

wplot(Filters, col =(HQ_samples+2), sub = "We discard the lower quintile (20%) of the data in both dimensions")
abline(h=log2(minReadCount), lty =3, col =2)
abline(v=log2(minSNP), lty =3, col =2)
plotname = kollapse("FigS.2_", plotnameLastPlot)
wplot_save_this(plotname = plotname, mdlink = T)

text(x = Filters, labels = rownames(Filters), cex =.5, srt =-45, pos =3)
wplot_save_this(plotname = kollapse(plotname, "_names"))

# table(substr(which_names(HQ_samples),7,8))
# table(substr(ControlNamez_HQ,7,8))

# which_names(HQ_by_Transcriptome)
ControlNamez_HQ = which_names(Control & HQ_samples)
GermNamez_HQ = which_names( !Control & HQ_samples)
# which_names(HQ_samples)
# metadata_manual_annot$

# -------------------------------------------------------------------------------------------------------------------------------------
llprint ("## The filtered dataset")
llprint (sum(HQ_samples), "samples remain for analysis. Further filtering is applied at each analysis step separately.
			 The cells are annotated as the following:")
CellTypes = table(cell_type_annot[HQ_samples]);
wpie (CellTypes, percentage = F, mdlink = T)

DP_SUM_parental_HQ = DP_SUM_parental[, HQ_samples]

llprint ("### Read Count statistics for the filtered dataset")
ReadCounts_HQ = colSums(DP_SUM_parental_HQ)
subb = kollapse ("Median Read Count: ", median (ReadCounts_HQ), xlab ="Read Count")
whist (ReadCounts_HQ, breaks = 25, sub = subb, mdlink = T)

properX_ReadCounts_HQ = colSums(DP_SUM_parental_HQ[gene_category == "X", ])
subb = kollapse ("Median Read Count: ", median (properX_ReadCounts_HQ))
whist (properX_ReadCounts_HQ, breaks = 25, sub = subb, mdlink = T)

table(gene_category)
Impr_ReadCounts_HQ = colSums(DP_SUM_parental_HQ[ (gene_category == "Impr_PAT" | gene_category == "Impr_MAT") , ])
subb = kollapse ("Median Read Count: ", median (Impr_ReadCounts_HQ))
whist (Impr_ReadCounts_HQ, breaks = 25, sub = subb, mdlink = T)

# Write out Metadata ----------------------------------------------------------------------------------
create_set_OutDir("~/Google_Drive/X_react_Data/Abelz/Metadata2")

metadata_cells = cbind(  cell_IDs, cell_NR, HQ_samples, Control, gender,
	weeks,  donor, metadata_manual_annot, OneLetterAnnot, cell_type_col, cell_type_col_simple )

write.simple.tsv(metadata_cells)

metadata_genes = cbind (chrz, gene_category);
write.simple.tsv (metadata_genes)

OutDir = OutDirOrig

# Derive Moniallelic Expression ----------------------------------------------------------------------------------
MonAllExp_pat = PP_matrix
MonAllExp_pat[PP_matrix > 0.98 & !is.na(PP_matrix)] = 'Pat'
MonAllExp_pat[PP_matrix < 0.02 & !is.na(PP_matrix)] = 'Mat'
MonAllExp_pat[ PP_matrix < 0.98 & PP_matrix > 0.02 &  !is.na(PP_matrix)] = 'Both'
BiallGenes = colSums(MonAllExp_pat == 'Both', na.rm = T)
MonAllExp_pat[MonAllExp_pat == NaN] =  NA
ExprGenes =	colSums(MonAllExp_pat != 'StringXXX', na.rm = T) # Colsums only works if you remove, NA-s. Then you want to count every value, hence the toy string.
ABR = BiallGenes / ExprGenes

clean = which_names(gene_category =="X" | gene_category == "autosomal"); l(clean)

DP_pat_clean =DP_pat[clean,]
DP_mat_clean =DP_mat[clean,]; dim(DP_pat_clean)
chrz_clean = chrz[clean]
chromosomes = sort(names(table (chrz_clean)))
nr_chr = l(chromosomes)

if (SecondPart) {

	# plot PNR and PP -----------------------------------------------
	# source("~/analysis/DevCell_analysis/13a.PNR_PP_plots.R")

	# 13.a.TestOfReactivation_AgainstSomaticX.R -----------------------------------------------
	# source ("~/analysis/DevCell_analysis/13.a.TestOfReactivation_AgainstSomaticX.R")

	# 13b.Parental_scatter_plots -------------------------------------------------------------------------
	source("~/analysis/DevCell_analysis/13.b.Parental_scatter_plots.R")

	# 13.d.Proportion_of_Biallelic_Genes_per_Chr ---------------------------------------------------------
	# source ("~/analysis/DevCell_analysis/13.d.Proportion_of_Biallelic_Genes_per_Chr.R")

	# Chromosome Sums ------------------------------------------

	ChrSums_Mat = aggregate(DP_mat_clean, by=list(chrz_clean),  FUN=sum, na.rm=TRUE);
	ChrSums_Mat=	data.matrix(ChrSums_Mat[,2:(nr_cells+1)])
	ChrSums_Mat= ChrSums_Mat[,HQ_samples]							# HQ samples only
	rownames (ChrSums_Mat) = chromosomes

	ChrSums_Pat = aggregate(DP_pat_clean, by=list(chrz_clean),  FUN=sum, na.rm=TRUE);
	ChrSums_Pat=	data.matrix(ChrSums_Pat[,2:(nr_cells+1)])
	ChrSums_Pat= ChrSums_Pat[,HQ_samples]							# HQ samples only
	rownames (ChrSums_Pat) = chromosomes

	ChrSums = ChrSums_Pat + ChrSums_Mat

	write.simple.tsv (ChrSums_Pat, ManualName = kollapse(OutDirOrig,"/ChrSums_Pat.tsv"))
	write.simple.tsv (ChrSums_Mat, ManualName = kollapse(OutDirOrig,"/ChrSums_Mat.tsv"))


	# 13.e.Chr_wide_scatters ------------------------------------------------------------------------------
	# source ("~/analysis/DevCell_analysis/13.e.Chr_wide_scatters.R")

	# 13.f How many genes determine the chrom sum? Is it only a handful? ----------------------------------
	# source ("~/analysis/DevCell_analysis/13.f.Gene_Composition.R")

	# 13.f2 ----------------------------------
	# source ("~/analysis/DevCell_analysis/13.f2.GeneCompositionNew.R")

	# 13.h chr X Reactivation testing ---------------------------------------------------------------------
	source ("~/analysis/DevCell_analysis/13.h.chrX_Reactivation.R")

	# 13.h.5 GEOM MEAN chr X Reactivation testing ---------------------------------------------------------------------
	# source ("~/analysis/DevCell_analysis/13.h5.ReactivationTest_GeomMean.R")
	# 13.h2 chr X log2 Reactivation testing ---------------------------------------------------------------------
	# source ("~/analysis/DevCell_analysis/13.h2.ReactivationTest_log2.R")
	# 13.h3 MEDIAN based chr X Reactivation testing ---------------------------------------------------------------------
	# source ("~/analysis/DevCell_analysis/13.h3.chrX_Reactivation_MEDIAN.R")

	# 13.k.chrX_chrBoxplots.R -----------------------------------------------------------------------------
	source ("~/analysis/DevCell_analysis/13.k.chrX_chrBoxplots.R")

	# 13.b2.PlotAllChrX_on1plot.R -----------------------------------------------------------------------------
	source ("~/analysis/DevCell_analysis/13.b2.PlotAllChrX_on1plot.R")

	# 13.m.EscapeeGenePrediction: Detect genes that are expressed from the inactive, or from  both alleles in Somatic cells ---------------------
	# source ("~/analysis/DevCell_analysis/13.m.EscapeeGenePrediction.R")

	### Is the expression of late cells universally lower? -----------------------------------------------------------------------------------------------------------------------------------------
	# source("~/analysis/DevCell_analysis/19.DownRegulationInLGC.R")

	# ChrRead Count for SNP containing reads ------------------------------------------------------------------------------------
	AutoBialleleicRate = AutoMatRate = AutoPatRate = numeric(nr_cells)
	XBialleleicRate  = XMatRate = XPatRate = numeric(nr_cells)
	XescBialleleicRate  = XescMatRate = XescPatRate = numeric(nr_cells)
	ImprBialleleicRate  = ImprMatRate = ImprPatRate = numeric(nr_cells)
	index_auto = 	which_names (gene_category == "autosomal")
	indexImpr = 	which_names (gene_category == "Impr_MAT" | gene_category == "Impr_PAT")
	indeXX = 	 	which_names (gene_category == "X")
	indeXesc = 	 	which_names (gene_category == "escapees" | gene_category == "escapees_het")
	MonAllExp_pat[MonAllExp_pat == NaN] =  NA

	for (j in 1:nr_cells) {
		AutoBialleleicRate[j] =	percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ index_auto,j ]),  'Both'), digitz=6)
		AutoPatRate[j] =		percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ index_auto,j ]),  'Pat'), digitz=6)
		AutoMatRate[j] = 		percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ index_auto,j ]),  'Mat'), digitz=6)
		XBialleleicRate[j] = 	percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indeXX ,j]),  'Both'), digitz=3)
		XMatRate[j] = 			percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indeXX,j ]), 'Mat'), digitz=3)
		XPatRate[j] = 			percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indeXX,j ]), 'Pat'), digitz=3)
		XescBialleleicRate[j] = percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indeXesc ,j]),  'Both'), digitz=3)
		XescMatRate[j] = 		percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indeXesc,j ]), 'Mat'), digitz=3)
		XescPatRate[j] = 		percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indeXesc,j ]), 'Pat'), digitz=3)
		ImprBialleleicRate[j] = percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indexImpr,j ]),  'Both'), digitz=3)
		ImprMatRate[j] = 		percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indexImpr,j ]), 'Mat'), digitz=3)
		ImprPatRate[j] = 		percentage_formatter(pc_in_total_of_match(table(MonAllExp_pat[ indexImpr,j ]), 'Pat'), digitz=3)
	}

	ReadCount_per_Chr_w_SNPs_box = remove_outliers(t(ChrSums_Mat + ChrSums_Pat))
	wboxplot (ReadCount_per_Chr_w_SNPs_box, mdlink = T)

	OutDir = OutDirOrig
	ST2_TranscriptStats = cbind (
		"Sample_Name"	= Sample_Name,
		"HQ_samples"	= HQ_samples,
		"SNPContGenes"		= colSums(DP_SUM_parental != 0, na.rm =T),
		"NR_X_SNPs_per_Sample"		= NR_X_SNPs_per_Sample,
		"NR_Impr_SNPs_per_Sample"		= NR_Impr_SNPs_per_Sample,
		"AutoBialleleicRate" 	= AutoBialleleicRate,
		"ImprBialleleicRate" 	= ImprBialleleicRate,
		"XBialleleicRate" 	= XBialleleicRate
		) # View(ST2_TranscriptStats)
	write.simple.tsv(ST2_TranscriptStats)

	xmat = DP_mat[ gene_category=="X", ]
	xpat = DP_pat[ gene_category=="X", ]
	MatXisActive = colSums(xmat) > colSums(xpat)

	ST2a_MonoAllGeneStats = cbind (
		"Annotation" = cell_type_annot,
		"HQ_samples" = HQ_samples,
		"Total_SNPs" = colSums(DP_SUM_parental != 0, na.rm =T),
		"MaternalActiv" = MatXisActive,
		"NR_X_SNPs_per_Sample"		= NR_X_SNPs_per_Sample,
		"NR_auto_per_Sample"		= NR_auto_per_Sample,
		"AutBi" = AutoBialleleicRate,
		"AutMat" = AutoMatRate,
		"AutPat" = AutoPatRate,
		"NR_X_proper_SNPs_per_Sample" = NR_X_proper_SNPs_per_Sample,
		"X_Bi" = XBialleleicRate,
		"X_Mat" = XMatRate,
		"X_Pat" = XPatRate,
		"NR_X_escapee_per_Sample" = NR_X_escapee_per_Sample,
		"Xesc_Bi" = XescBialleleicRate,
		"Xesc_Mat" = XescMatRate,
		"Xesc_Pat" = XescPatRate
	)
	write.simple.tsv(ST2a_MonoAllGeneStats)

	#  ------------------------------------------------------------------------------------------------------------------------------------------------------
	source("~/analysis/DevCell_analysis/13.ST2b_EachXGene.R")
	# source("~/analysis/DevCell_analysis/13.ST3b_EachImprGene.R")

} # if (SecondPart)

# ICR annotation - only once in a lifetime!!! ------------------------------------------
# source("~/Dropbox/X_reactivation/Metadata_d/Imprinted_genes/ICs/zz.ICR_gene_annotation")
memspace_file = kollapse(OutDir, "/memspace12.rd")
save.image(file=memspace_file)

MiniStats = T
if (MiniStats) {
  X_geneIndex = which(rowSums(DP_SUM_parental_HQ[ which(chrz == "chrX" ), ] > 0)>0)
  NrXgens = l(X_geneIndex);NrXgens
  Nr_X_proper = l(setdiff(setdiff(names(X_geneIndex), escapees), escapees_het)); Nr_X_proper
  MiniStats = c(NrXgens = NrXgens, Nr_X_proper = Nr_X_proper, Nr_Esc = NrXgens-Nr_X_proper)
  MarkDown_Table_writer_NamedVector(MiniStats)
}
