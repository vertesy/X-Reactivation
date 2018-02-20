######################################################################################################
# Box Plots without Exome seq
######################################################################################################
# source ("~/analysis/DevCell_analysis/13.l.chrX_chrBoxplots.BiallelicRate.REBUTTAL.R")
try(dev.off())
assign("last.warning", NULL, envir = baseenv())

# Setup -------------------------------------------------------------------------------------
OutDir = OutDirOrig ; OutDir
scriptname = "13.l.chrX_chrBoxplots.BiallelicRate.REBUTTAL.R"
setup_MarkdownReports(OutDir = kollapse(OutDirOrig,"/13.l.chrX_chrBoxplots.BiallelicRate.REBUTTAL.R"), scriptname = scriptname, append = F)

# Percent_of_Biallelic_Genes -------------------------------------------------------------------------------------


  ST2_TranscriptStats = read.simple.tsv("~/Google_Drive/X_react_Data/Abelz/X_react_Paper/ST2_TranscriptStats.Fraction.tsv")
  attach_w_rownames(ST2_TranscriptStats)

  rm(list(colnames(ST2_TranscriptStats)))

  Celltypes_= (cell_type_annot_simplest[rownames(ST2_TranscriptStats)])
  Percent_of_Biallelic_Genes = split(100*ST2_TranscriptStats$XBialleleicRate, f = Celltypes_)

  display_order =  c("SOMATIC", "EGC", "L+MGC")
  Percent_of_Biallelic_Genes = Percent_of_Biallelic_Genes[display_order]
  wstripchart(Percent_of_Biallelic_Genes, col=CellTypeColz_Simplest[display_order], colorbyColumn = T, jitter = .35, ylim =c(0,100), main="Without Haplotypes Technical Bias Masks Reactivation",
              ylab="Percentage of Biallelically Expressed Genes on chrX", sub = "Biallelic expression defined as less than 95% bias towards either allele")

  oo()

# Same Graph in simulated Reference-Alternative dimension -------------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------------------------------------------------

  InputDir = "~/Google_Drive/X_react_Data/RNA/Intersections/gVCF/per_rna_sample/chrx/"
  ls_F = list.files(InputDir, pattern = ".txt$")
  cellNamezz = substr(ls_F, start = 0, stop = 8)

  nrf = l(ls_F)

  f=1

  Coverage = BiallelicRate = Ls_cells= NULL
  for (f in 1:nrf) {
    print(f)
    Ls_cells[[f]] = read.simple.table(InputDir,ls_F[f])
    BiallelicRate[f] = sum(Ls_cells[[f]]$"alt_dp", na.rm = T) / sum(Ls_cells[[f]]$"total_dp", na.rm = T)
  }

  names(BiallelicRate) = cellNamezz

  symdiff(names(cell_type_annot_simplest), cellNamezz)

  foundd = intersect(names(cell_type_annot_simplest), cellNamezz)
  l(cell_type_annot_simplest[foundd])

  BiallelicRate_ls = split(100*BiallelicRate, f=cell_type_annot_simplest[foundd])
  BiallelicRate_ls = BiallelicRate_ls[display_order]


  wstripchart(BiallelicRate_ls, col=CellTypeColz_Simplest[display_order], colorbyColumn = T, jitter = .35, ylim =c(0,100),
              main="Pooling Information from Individual SNPs is not possible without ",
              ylab="%ALT reads on chrX (average)", sub = "Pooling information is problematic without haplotype")


  # ------------------------------------------------------------------------------------------------------------------------



# Setup -------------------------------------------------------------------------------------

GeoMeanTesting = F
# It corresponds to certain settings in 13.h
if (GeoMeanTesting) {
	log_settings_MarkDown(GeoMeanTesting, minXGMean, minGenes, SignLevel, NrBins)
	PP_of_metric = PP_ChrGeomMean_Total
	Test = "GeoMean-based test"
} else { # if Sum based testing
	PP_of_metric = PP_chr_level
	log_settings_MarkDown(GeoMeanTesting, minGenes, SignLevel, minXReads, rc_thr)
	Test = "SUM-based test"
}



# Calculate New metric -------------------------------------------------------------------------------------
dim(PP_matrix)

isMonoallelic = (PP_matrix >= .95 | PP_matrix <= .95)
xx = aggregate(isMonoallelic, by=list(chrz),  FUN=mean, na.rm=TRUE)[ ,-1]


llprint ("# Boxplots of Allelic bias in chromosome X")
llprint ("## ", Test)
llprint ("There are four kind of boxplots can be made.")
llprint ("As data we can use the **percentile position** of chr-X in the autosomal distribution or we can use the  **% of paternal reads**.")
llprint ('We can plot these numbers **"Raw"** where you have to comapre the differences in variation, which is bit unusual
			 and hard to test for significance, as this is not normally distributed data. We chose for *Mood\'s Test*, but we also calculate the
			 significance at the *Ansari-Bradley Test*, which yields typically very similar results.')
llprint ("Alternatively we can reduce the dimensions, and bring it to **Biases-to-Unbiased scale**, where you compare means ( & test with *Mann-Whitney-Wilcoxon*),
			 but you lose the paternal-maternal information. The transformation is essentially folding the raw distribution into half. )")

# Setup ------------------------------------------------------------------------------------------------

bxpOrder = c("SOMATIC","EGC", "L+MGC")
# cell_type_annot_simplest = cell_type_annot_simple
# cell_type_annot_simplest[cell_type_annot_simplest=="MGC" | cell_type_annot_simplest=="LGC" ] = 'L+MGC'
# metadata_manual_annot = read.simple.tsv("~/Google_Drive/X_react_Data/Abelz/Metadata2/Manual_Metadata/metadata_manual_annot_AllCells.TPM.TW_clustering.tsv")
metadata_manual_annot = read.simple.tsv("~/Google_Drive/X_react_Data/Abelz/Metadata2/Manual_Metadata/metadata_manual_annot_AllCells.PGC_clustering.tsv")
attach_w_rownames(metadata_manual_annot);

CellTypeColz_Simplest = CellTypeColz_Simple[-4]
# names(CellTypeColz_Simplest)[3] = 'L+MGC'

# The 4 plots in this script:
PP_PatMat_plot = T
PP_BiasUnbias_plot = T
Percentile_BiasUnbias_plot = T
Percentile_PatMat_plot = T

WithMean = T
ExcludeMales  = F
ControlForOutliers = F # Do not check TCEAL genes effect

# Prepartaion -------------------------------------------------------------------------------------
# PP_of_metric["chrX","d1_03_gf"]
QCsum = F
if (QCsum) {
	cxx = rep ("grey33", ncol(PP_of_metric)); l(cxx)
	colnames(PP_of_metric) ==HQ_Cellnames

	cxx[HQ_Cellnames %in% hitzzz] = 2; table(cxx)
	ppChrSum = PP_of_metric["chrX",]
	idx = order(ppChrSum)
	wbarplot(ppChrSum[idx], col=cxx[idx])
}

ChrX_Allelic_Expression_Bias = PP_of_metric[23, which_names(HQ_samples) ]
percentiles_of_chrX_Bias = percentiles_of_chrX_ALL

percentiles_of_chrX_Bias[which_names(fail_X_QC)] =NaN
percentiles_of_chrX_Bias = na.omit.strip(percentiles_of_chrX_Bias[pass_X_QC_names])

ChrX_Allelic_Expression_Bias = 	na.omit.strip(ChrX_Allelic_Expression_Bias[pass_X_QC_names]) 				# Remove NA-s from PP
names_tested = names (ChrX_Allelic_Expression_Bias)

Categories = cell_type_annot_simplest[names_tested]
# Categories = cell_type_annot_simplest[names(percentiles_of_chrX_Bias)]

LinTransform = F
if(LinTransform) {	RelativeAllelicBias = rescale	(ChrX_Allelic_Expression_Bias, from = 0, upto = 100); "this makes the mean shif to 60  "
} else { RelativeAllelicBias = ChrX_Allelic_Expression_Bias}

PPforBXP_ls =split(RelativeAllelicBias, Categories)[bxpOrder] 				# Split into a list per annotation category
	PPforBXP_ls_NR = 	split(RelativeAllelicBias[hitzzz], 	Categories[hitzzz])[bxpOrder] 				# Split into a list per annotation category
	PPforBXP_ls_R = 	split(RelativeAllelicBias[non_hitzzz], Categories[non_hitzzz])[bxpOrder] 				# Split into a list per annotation category

PercentileforBXP_ls =	 split(percentiles_of_chrX_Bias, Categories)[bxpOrder]
	PercentileforBXP_ls_ls_NR =	 split(percentiles_of_chrX_Bias[hitzzz], 	 Categories[hitzzz])[bxpOrder] 				# Split into a list per annotation category
	PercentileforBXP_ls_ls_R =	 split(percentiles_of_chrX_Bias[non_hitzzz], Categories[non_hitzzz])[bxpOrder] 				# Split into a list per annotation category

Control_HQ = which_names(Control[names_tested])
Germ_HQ =    which_names(!Control[names_tested])
# if (ExcludeMales) {  }
PP_mean_of_categories =	unlist(lapply(PPforBXP_ls,  mean)); PP_mean_of_categories
Percentile_mean_of_categories =unlist(lapply(PercentileforBXP_ls,  mean, na.rm = T)); Percentile_mean_of_categories

mean(PercentileforBXP_ls$SOMATIC, na.rm = T)
# ---------------------------------------------------------------------------------------------------------------------------------------
# Boxplot1: Chr-level % of paternal reads with mean -------------------------------------------------------------------------------------
llogit("-----------")
llprint("## Percentage of Paternal Reads")
llprint(l(ChrX_Allelic_Expression_Bias), "cells have a valid PP ratio, these are plotted (Males NOT excluded).")

fname = "1. Percent of Paternal reads on chr-X"
llprint ("### ", fname)
llprint("Plotted with mean")

# http://en.wikipedia.org/wiki/Variance#Tests_of_equality_of_variances
MOODtest = mood.test(percentiles_of_chrX_Bias[Control_HQ], percentiles_of_chrX_Bias[Germ_HQ])
ANStest = suppressWarnings(ansari.test(percentiles_of_chrX_Bias[Control_HQ], percentiles_of_chrX_Bias[Germ_HQ]))
AnsP = iround(	as.numeric(as.vector(ANStest[2]))	)
MoodP = iround(	as.numeric(as.vector(MOODtest[2]))	)
llprint("The difference is significant on both applicable test known to me. P-value at Mood's test:", MoodP ,"and at Ansari's test:", AnsP)

pMain = "Paternal-Maternal Bias in chr X"
pSub = kollapse('Germ vs Somatic p-value @ Mood\'s Test: ', MoodP )
pYlim = range(PPforBXP_ls)
xlb = "Cell types"
ylb = "% reads from the paternal chromosome"
# ylb = " gm(PatRC)/(gm(PatRC) +gm(MatRC))"
# coll = CellTypeColz_Simplest[names(PPforBXP_ls)]
coll = CellTypeColz_Simplest[names(PPforBXP_ls)]

a =boxplot(PPforBXP_ls, plot=F)
if (WithMean) { a$stats[3,] = PP_mean_of_categories }						# Replace mean with median
bxp(a, outpch = NA, ylim = pYlim, main=pMain, sub=pSub, xlab = xlb, ylab= ylb)
stripchart(PPforBXP_ls_R, vertical = TRUE, # method = "stack",
		   method = "jitter", jitter = .35,
		   pch = 16, add = TRUE, cex=1.5, col= "grey33" )
stripchart(PPforBXP_ls_NR, vertical = TRUE, # method = "stack",
		   method = "jitter", jitter = .35,
		   pch = 16, add = TRUE, cex=1.5, col= coll )

wplot_save_this (fname, mdlink = T)

# Boxplot2: Chr-level % of Active reads with mean -------------------------------------------------------------------------------------
fname = "2. Percent of reads from Active allele on chr-X"
llprint ("### ", fname)
llprint("Plotted with median")

ChrX_Allelic_Expression_Bias[ChrX_Allelic_Expression_Bias < 50] =100-ChrX_Allelic_Expression_Bias [ChrX_Allelic_Expression_Bias < 50] 		# Transforming Pat-Mat to Unbiased-Biased
PPforBXP_ls =split(ChrX_Allelic_Expression_Bias, Categories)[bxpOrder] 				# Split into a list per annotation category

# https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
MWW = suppressWarnings(wilcox.test( ChrX_Allelic_Expression_Bias[Control_HQ], ChrX_Allelic_Expression_Bias[Germ_HQ]))
any_print ( 'p-value @ MWW: ', iround(as.numeric(MWW[3]),2) )

pMain = "Allelic Bias in chr X"
pSub = kollapse('Germ vs Somatic p-value @ MWW: ', iround(as.numeric(MWW[3]),2) )
pYlim = range(PPforBXP_ls);pYlim
coll = CellTypeColz_Simplest[names(PPforBXP_ls)];coll
xlb = "Cell types"

if (GeoMeanTesting) {
	ylb = " gm(ActiveRC)/(gm(ActiveRC) +gm(SilentRC))"
} else { # if Sum based testing
	ylb = "% of reads from the active allele of chr-X"
}

a =boxplot(PPforBXP_ls, plot=F)
# if (WithMean) { a$stats[3,] = PP_mean_of_categories }						# Replace mean with median
bxp(a, outpch = NA, ylim = pYlim, main=pMain, sub=pSub, xlab = xlb, ylab= ylb)
stripchart(PPforBXP_ls, vertical = TRUE, method = "jitter", jitter = .3,
		   pch = 16, col = coll, add = TRUE, cex=1.5 )
wplot_save_this (fname, mdlink = T)

llprint ("I was thinking about median normalization, but then for control, you would have values around 140%, so I abandoned the idea.")

# ---------------------------------------------------------------------------------------------------------------------------------------
# Boxplot3: Boxplot with percentile rank -------------------------------------------------------------------------------------
llogit("---------------")
llprint("## Percentile in within the null-distribution")
llprint(l(percentiles_of_chrX_Bias), "cells are plotted for percentile position (Males excluded).")

fname = "3. Percentile position of chr-X in the autosomal distribution"
llprint ("### ", fname)
llprint("Plotted with mean")

# http://en.wikipedia.org/wiki/Variance#Tests_of_equality_of_variances
MOODtest = mood.test(percentiles_of_chrX_Bias[Control_HQ], percentiles_of_chrX_Bias[Germ_HQ])
ANStest = suppressWarnings(ansari.test(percentiles_of_chrX_Bias[Control_HQ], percentiles_of_chrX_Bias[Germ_HQ]))
AnsP = iround(	as.numeric(as.vector(ANStest[2]))	)
MoodP = iround(	as.numeric(as.vector(MOODtest[2]))	)
llprint("The difference is significant on both applicable test known to me. P-value at Mood's test:", MoodP ,"and at Ansari's test:", AnsP)

pMain = "Paternal-Maternal percentile of chr X"
pSub = kollapse('Germ vs Somatic p-value @ Mood\'s Test: ', MoodP )
pYlim = range(PercentileforBXP_ls, na.rm = T);pYlim
coll = CellTypeColz_Simplest[names(PPforBXP_ls)]
xlb = "Cell types"
ylb = "Percentile in the total distributiion"

a =boxplot(PercentileforBXP_ls, plot=F)
if (WithMean) { a$stats[3,] = Percentile_mean_of_categories }						# Replace mean with median
bxp(a, outpch = NA, ylim = pYlim, main=pMain, sub=pSub, xlab = xlb, ylab= ylb)
stripchart(PercentileforBXP_ls_ls_R, vertical = TRUE, # method = "stack",
		   method = "jitter", jitter = .3,
		   pch = 16, add = TRUE, cex=1.5, col= "grey33" )
stripchart(PercentileforBXP_ls_ls_NR, vertical = TRUE, # method = "stack",
		   method = "jitter", jitter = .3,
		   pch = 16, add = TRUE, cex=1.5, col= coll )
# stripchart(PercentileforBXP_ls, vertical = TRUE, method = "jitter",
# 		   pch = 16, col =coll, add = TRUE, cex=1.5 )
wplot_save_this (fname, mdlink = T)
# dev.off()


# Boxplot4: Boxplot with percentile rank -------------------------------------------------------------------------------------
fname = "4. Symmetric Percentile position of chr-X in the autosomal distribution"
llprint ("### ", fname)
llprint("Plotted with median")

percentiles_of_chrX_Bias = percentiles_of_chrX_Bias[!is.na(percentiles_of_chrX_Bias)]
percentiles_of_chrX_Bias[percentiles_of_chrX_Bias < 50] =100-percentiles_of_chrX_Bias [percentiles_of_chrX_Bias < 50] 		# Transforming Pat-Mat to Unbiased-Biased
percentiles_of_chrX_Bias = rescale(percentiles_of_chrX_Bias, from = 0, upto = 100)

PercentileforBXP_ls = 	split(percentiles_of_chrX_Bias, Categories[names(percentiles_of_chrX_Bias)])[bxpOrder] 															# Split into a list per annotation category


# https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
MWW = suppressWarnings(wilcox.test( percentiles_of_chrX_Bias[Control_HQ], percentiles_of_chrX_Bias[Germ_HQ]))
any_print ( 'p-value @ MWW: ', iround(as.numeric(MWW[3]),2) )

pMain = "Allelic Bias of chr X"
pSub = kollapse('Germ vs Somatic p-value @ MWW: ', iround(as.numeric(MWW[3]),2) )
pYlim = range(PercentileforBXP_ls, na.rm = T);pYlim
coll = CellTypeColz_Simplest[names(PPforBXP_ls)]
xlb = "Cell types"
ylb = "Symmetric % in the total distributiion"

a =boxplot(PercentileforBXP_ls, plot=F)
# if (WithMean) { a$stats[3,] = Percentile_mean_of_categories }						# Replace mean with median
bxp(a, outpch = NA, ylim = pYlim, main=pMain, sub=pSub, xlab = xlb, ylab= ylb)
stripchart(PercentileforBXP_ls, vertical = TRUE, method = "jitter", jitter = .3,
		   pch = 16, col = coll, add = TRUE, cex=1.5 )
wplot_save_this (fname, mdlink = T)

##############################################################################################################
if (ControlForOutliers) {
	WeirdControls_Perc = head(sort(PercentileforBXP_ls$Control), n=2)
	llprint("The 2 outlier controls are:", names(WeirdControls_Perc), "The same as in the PP plots")

	llogit("---------------")
	llprint("## Control for outliers")
	WeirdControls_Perc

	llprint("d2_18_af indeed looks terrible. See the linked image below:")
	llogit(	MarkDown_ImgLink_formatter(OutDirOrig, "/Scatter_PAT_log2/SNP_scatter_log2_d2_18_af_Control.plot.pdf") 	)

	llprint("In d2_35_gf there is on highly expressed, beautifully biallelic gene is causign the effect. (And possibly XIST - to be solved)")
	llogit(	MarkDown_ImgLink_formatter(OutDirOrig, "/Scatter_PAT_log2/SNP_scatter_log2_d2_35_gf_Control.plot.pdf") 	)


	OutlierzGeenz = DP_SUM_parental_HQ[, "d2_35_gf"]
	OutlierzGeenzNamez = which_names(gene_category[OutlierzGeenz>0] == "X")

	AllelicExprin_d2_35_gf = rbind (
				"DP_mat" = DP_mat_HQ[OutlierzGeenzNamez, "d2_35_gf"],
				"DP_pat" = DP_pat_HQ[OutlierzGeenzNamez, "d2_35_gf"],
				"DP_sum" = DP_SUM_parental_HQ[OutlierzGeenzNamez, "d2_35_gf"]
	)
	colnames(AllelicExprin_d2_35_gf) = OutlierzGeenzNamez

	MarkDown_Table_writer_DF_RowColNames(AllelicExprin_d2_35_gf)

	llprint("Well, TCEAL4 is a predicted escapee, but it is not listed in our annotation.")
	'TCEAL4' %in% escapees
	'TCEAL4' %in% escapees_het # it should be here

	tceal = as.vector(PP_matrix_HQ['TCEAL4',]);tceal

	TCEAL4_AllelicExpr = t(rbind(
		"DP_pat" = as.vector(DP_pat_HQ['TCEAL4',]),
		"DP_mat" = as.vector(DP_mat_HQ['TCEAL4',])
	))
	dim(TCEAL4_AllelicExpr)

	llprint("If I look at all the cells, is sometimes bi- sometimes mono-allelic, but it seems that all genes look like that (bursting/sampling):")
	cc = CellTypeColz[cell_type_annot_simple[cell_IDs[HQ_samples]]]
	wplot(TCEAL4_AllelicExpr, col = cc, mdlink = T, pch =16)
	text(TCEAL4_AllelicExpr, labels = cell_IDs[HQ_samples], cex = 0.5, pos = 3)
	wplot_save_this("TCEAL4_AllelicExpr.plot", ManualName = T)
} # if (ControlForOutliers)


# Autosomal Bias per cell category --------------------------------------------------------------------------------------------------------
llprint("# Autosomal Bias per cell category")
llprint("## Smaller autosomal bias in LGC-s explain the smaller PA bias in chr-X in the same cells.")

Categz = names(table(cell_type_annot_simplest));Categz

CellByCateg = list(
	"SOMATIC" = na.omit.strip(as.vector(PP_of_metric[1:22, which_names(HQ_samples & cell_type_annot_simplest == "SOMATIC")])),
	"EGC" = 	na.omit.strip(as.vector(PP_of_metric[1:22, which_names(HQ_samples & cell_type_annot_simplest == "EGC")])),
	"L+MGC" = 	na.omit.strip(as.vector(PP_of_metric[1:22, which_names(HQ_samples & cell_type_annot_simplest == "L+MGC")]))
)

library("vioplot"); # install.packages("vioplot")


pname = "Percentage of reads from the Paternal Allele in Autosomes"
par(las=1,bty="l")  ## my preferred setting
## set up empty plot
ylb = "Percentage of Paternal Reads [%]"
plot(0,0,type="n",xlim=c(0.5,3.5), ylim=c(0,100),  xaxt = 'n', xlab ="", ylab = ylb,  main =pname)
abline(h=c(10, 30, 50, 70, 90), lty=2)
for (i in 1:3) { vioplot(100*CellByCateg[[i]], at = i, add = T, col = CellTypeColz_Simplest[i]) }
axis(side=1,at=1:3,labels=Categz[3:1])

wplot_save_this(pname, mdlink = T)

# --------------------------------------------------------------------------------------------------------
llprint("## The percentile plots")
llprint("The percentile plot is not really meaningful: position in its own distribution. ")
llprint("Yet there is some trend between cell categories")


par(mfrow = c(3,2))  # 3 rows and 2 columns
Percentile_all_autsomes = 100*ecdf(as.vector(PP_of_metric[1:22, ])) (as.vector(PP_of_metric[1:22, ]))
whist(Percentile_all_autsomes)
AutosomalPercentiles = list (
	Percentile_SOMA_autsomes = ecdf(as.vector(PP_of_metric[1:22, ])) (as.vector(PP_of_metric[1:22, which_names(HQ_samples & cell_type_annot_simplest == "SOMATIC")])),
	Percentile_EGC_autsomes = ecdf(as.vector(PP_of_metric[1:22, ])) (as.vector(PP_of_metric[1:22, which_names(HQ_samples & cell_type_annot_simplest == "EGC")])),
	Percentile_LGC_autsomes = ecdf(as.vector(PP_of_metric[1:22, ])) (as.vector(PP_of_metric[1:22, which_names(HQ_samples & cell_type_annot_simplest == "L+MGC")]))
)
for (i in 1:3) { hist(100*AutosomalPercentiles[[i]], col = CellTypeColz_Simplest[i], main = names(AutosomalPercentiles)[i],
					  breaks=20, xlab = "Percentile in the distribution") }
abline(h=c(6,12))

llprint("### The LGC in not significanlty different from random expectation for that many autosomes")

NrOfAutsomes_in_LGC = l(AutosomalPercentiles[[i]])
PercentileExpectationRandom = sample(Percentile_all_autsomes, NrOfAutsomes_in_LGC)
whist(PercentileExpectationRandom)

NrTrials = 1000
SimulatedLGC = matrix(data = NaN, nrow = NrTrials, ncol = 20)
for (j in 1:NrTrials) {
	SimulatedLGC[j,] = hist(sample(Percentile_all_autsomes, NrOfAutsomes_in_LGC), breaks=20, plot=F)$counts
}

xx = hist(sample(Percentile_all_autsomes, NrOfAutsomes_in_LGC), breaks=20, plot=F)$mids
SimulatedMedian = colMedians (SimulatedLGC)
wbarplot(SimulatedMedian, errorbar = T, upper = apply(SimulatedLGC, MARGIN = 2, sd))
# aaa = barplot(SimulatedMedian, add=F, ylim =c(0,15), plot=F)
# error.bar( x = aaa, y =SimulatedMedian,upper =  apply(SimulatedLGC, MARGIN = 2, sd))

par(mfrow = c(1,1))
wplot_save_this(plotname = "Percentile_all_autsomes.hist", w = 8.2, h =11.69, ManualName = T, mdlink = T)

# Write out ranking ---------------------------------------------------------------
OutDir = OutDirOrig; OutDir
OutDir = create_set_OutDir(OutDir,"/Rankings"); OutDir

nr_HQ_Germs = l(percentiles_of_chrX_Bias[Germ_HQ])
ChrX_Reactivation = sort(percentiles_of_chrX_Bias[Germ_HQ], decreasing = T)
ChrX_Reactivation =  2*(100 - ChrX_Reactivation) 										# Scale up to 1-100 from early to late

# write out --------------------------------------------------
OutDir = OutDirOrig; OutDir

xRanking_PA = RelativeAllelicBias[cell_type_annot_simplest[names(RelativeAllelicBias)] != "SOMATIC" & !is.na(RelativeAllelicBias)]
xRanking_PA[xRanking_PA<.5] = 1-xRanking_PA[xRanking_PA<.5]

xRanking_PA = sort (xRanking_PA, decreasing = T)
xRanking_PA = cbind (xRanking_PA, Sample_Name[names(xRanking_PA)])
write.simple.tsv(xRanking_PA)

xRanking_Percentile = percentiles_of_chrX_Bias[cell_type_annot_simplest[names(percentiles_of_chrX_Bias)] != "SOMATIC" & !is.na(percentiles_of_chrX_Bias)]
xRanking_Percentile[xRanking_Percentile<50] = 100-xRanking_Percentile[xRanking_Percentile<50]
xRanking_Percentile = sort (xRanking_Percentile, decreasing = T)
xRanking_Percentile = cbind (xRanking_Percentile, "Sample_Name" = Sample_Name[names(xRanking_Percentile)])
write.simple.tsv(xRanking_Percentile)

PlotRanking = T
if (PlotRanking) {
  Nr = dim(xRanking_Percentile)[[1]]
  ccc = as.numeric(as.factor(donor[rownames(xRanking_Percentile)]))
  plot(rep(1,Nr), pch=19, cex=1, col=ccc, ylim = c(.9, 1.1), xlab = "", ylab = "", axes = F, main="X-ranking order")
  ccc = cell_type_col_simple[rownames(xRanking_Percentile)]
  points(rep(.95,Nr), pch=15, cex=1, col=ccc)
  wplot_save_this(plotname = "Fig2_f.X_Ranking")
}

