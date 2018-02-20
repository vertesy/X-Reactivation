######################################################################################################
# Test if chrX has a different allelic bias than autosomes
######################################################################################################
# source ("~/analysis/DevCell_analysis/13.h.chrX_Reactivation.R")
try(dev.off(), silent = T)

# Setup ----------------------------
OutDir = OutDirOrig; OutDir
OutDir = create_set_OutDir(OutDir,"/13.h.X-Reactivation"); OutDir
setup_logging_markdown ("13.h.chrX_Reactivation.R", append = F)
log2plot =T

# Parameters ----------------------------
Autosomes = (ChrSums[ -23,])
minXReads = min(Autosomes[Autosomes>0])
BinSize = 200
rc_thr = 10
minGenes = 3
SignLevel = 0.95
log_settings_MarkDown(minXReads, BinSize, minGenes, SignLevel)

ChrGeneCount = aggregate((DP_SUM_parental_HQ>rc_thr), by=list(chrz),  FUN=sum, na.rm=TRUE)[ ,-1]

ChrSums_rawChrAnnot = aggregate(DP_SUM_parental_HQ, by=list(chrz),  FUN=sum, na.rm=TRUE)[ ,-1]
ChrSums_rawChrAnnot = unlist(ChrSums_rawChrAnnot)
ChrGeneCount_ = unlist(ChrGeneCount)

# Testing of reactivation ----------------------------------------------------------------------
llogit (" ## Finding the best threshold of reliability for Chromosome X-reactivation test")
ChrSums = ChrSums_Mat[ ,HQ_Cellnames] + ChrSums_Pat[ ,HQ_Cellnames]
PP_chr_level = ChrSums_Pat[ ,HQ_Cellnames] / ChrSums
write.simple.tsv(PP_chr_level, ManualName = kollapse(OutDirOrig,"/PP_chr_level.tsv"))

tested = names(which(ChrSums[23,]>minXReads))
nr_chrX_ExprAboveThr = l(tested)

# Reliability ---------------------------------------------------------------------------------------------------------------------------------------------------------
llprint ("### False discovery in control cells is low, and independent of chr total expression level in control cells")
chrxex = ChrSums[23, ]

conn = intersect(which_names (chrxex > 0), which_names(Control))
# conn = intersect(conn , pass_X_QC_names);

xex = sort(ChrSums[23, conn]); xex

PP_chrX_sorted_by_X_expression = PP_chr_level[23, names(xex)]
PP_chrX_sorted_by_X_expression[PP_chrX_sorted_by_X_expression<.5] = 1-PP_chrX_sorted_by_X_expression[PP_chrX_sorted_by_X_expression<.5]
llprint ("Surprsingly for chrX there is no clear correlation between relaibility and expression level, **unlike the clear correlation observed in imprinted genes!**
			 Disreagrd male cells (d3). : ")
pname= "ChrX Allelic Expr. bias in Controls vs ReadCount"
wbarplot(PP_chrX_sorted_by_X_expression, ylab="Relative Allelic Expression bias", col="olivedrab3",
		 plotname = pname, sub = "Read Counts are labeled within the bars")
a=barplot(PP_chrX_sorted_by_X_expression, plot=F)
text (x=a, y=0.1, labels = xex, srt=90, cex=1)
wplot_save_this (plotname = plotnameLastPlot, mdlink = T)
llprint("Of course, X-inactivation has one more degree of freedom as compared to the e.g. Paternally imprited genes, where reads must come from specifically the maternal allele.
			Yet, I believe, this cannot explain the difference between X-inact and Imprinting.")

# Reliability in gene counts ---------------------------------------------------------------------------------------------------------------------------------------------------------
llprint ("### False discovery does is not affected by low genecounts either")

cc = unlist(ChrGeneCount[23,])
conn = intersect(which_names ( cc> 1), which_names(Control)); l(conn)
xex = sort(ChrGeneCount[23, conn]); xex

PP_chrX_sorted_by_X_expression = PP_chr_level[23, names(xex)]
PP_chrX_sorted_by_X_expression[PP_chrX_sorted_by_X_expression<.5] = 1-PP_chrX_sorted_by_X_expression[PP_chrX_sorted_by_X_expression<.5]
pname= "ChrX Allelic Expr. bias in Controls vs GeneCount"
wbarplot(PP_chrX_sorted_by_X_expression, ylab="Relative Allelic Expression bias", col="olivedrab3",
		 plotname = pname, sub = "Observed Gene Counts are labeled within the bars")
a=barplot(PP_chrX_sorted_by_X_expression, plot=F)
text (x=a, y=0.1, labels = xex, srt=90, cex=1)
wplot_save_this (plotname = plotnameLastPlot, mdlink = T)
llprint ("Nevertheless we must exclude cells having a low number of genes. Most genes are monoallelically detected,
			 so a chromosome with e.g. 2 monoallelic genes have 50% chance to show consistent monoallelic expression for the whole chromosome.")

# Filtering  ---------------------------------------------------------------------------------------------------------------------------------------------------------
GenesOnChrX = unlist(ChrGeneCount[23, ])
fail_X_QC = ( GenesOnChrX < minGenes)
NrPass = sum(fail_X_QC)
fail_X_QC_names = which_names(fail_X_QC)
pass_X_QC_names = which_names(!fail_X_QC)

llprint(sum(fail_X_QC), "samples fail the X-chromosome QC, they are removed from the analysis. Cells have to have more than:", minGenes ,"genes expressed on chromosome X.")
CellsPassing_xQC = table(cell_type_annot_simple[pass_X_QC_names])
MarkDown_Table_writer_NamedVector(CellsPassing_xQC)

llprint("### The removed cells are:")
CellsFail_xQC = table(cell_type_annot_simple[fail_X_QC_names])
MarkDown_Table_writer_NamedVector(CellsFail_xQC)

fail_X_QC_names_germ = intersect(fail_X_QC_names, GermNamez_HQ)
pass_X_QC_names_germ = intersect(pass_X_QC_names, GermNamez_HQ)
llprint("HQ germ cells of these:", fail_X_QC_names_germ)


# llogit (" ## Filtering on minimal number of genes")
# llprint ("Filter out cells having less than", minGenes, "genes expressed above",rc_thr,"reads. For chromosomes above this,
# 	 the false discovery [false concordant monoallelic rate] due to integer effect is fairly low: 0.5 ^ 3 = 12.5%.
# 	 This is assuming that all detected genes are coming from a single transcript, which is the worst possible scenario,
# 	 and certainly unlikely for highly expressed genes.")
# llprint ("We discard", NrPass, "cells from the dataset. **Cells passing:**", pass_X_QC_names)
# llprint ("The discarded cells have max.", max(ChrSums[23, names(which(fail_X_QC))]), "reads on ChrX")
# whist(GenesOnChrX, breaks = max(GenesOnChrX), vline = minGenes-1, filtercol = 1)

# minXReads = min(chrxex[pass_X_QC_names]); minXReads
# llprint("**This filtering as a side effect removes cell with less than",minXReads, "reads.**")
#
#
# # Tricky correlation with the minimal gene Filtering  ---------------------------------------------------------------------------------------------------------------------------------------------------------
# llprint("One could filter based on how many genes determine the transcriptome of each give chromosome, hence exclude based on actual transcripts that determine the chr allelicity,
# *and not on genes above some arbitraty threshold.*")
#
#
# # Gene and Read ount correlations  ---------------------------------------------------------------------------------------------------------------------------------------------------------
# # Why can we exclude lowly expressed chromosomes? > Their read counts often come from a very few genes only > They have a high chance to sit on the axes----------------------------------------------------------------------
# "What is the correlation between the number of reads per Chr and the nr of observed genes ?"
# GeneCountVSReadCount = T
# if (GeneCountVSReadCount) {
# 	llprint ("### Correlation between read count and gene count underpins to remove LE chromosomes")
# 	llprint ("1. Best correlation is observed fore genes expressed 25+ reads.", rc_thr, "reads are chosen for analysis of the LE regime" )
# 	llprint ("- Not much correlation with samples [not depicted], some noticable correlation with certain chromosomes. This is probably explainable by gene density.")
#
# 	chr_colz = matrix(data = rainbow(nr_chr), nrow = nr_chr, ncol = nr_HQ_cells) # 	cell_colz = t(matrix(data = rainbow(nr_HQ_cells), nrow = nr_HQ_cells, ncol =nr_chr ))
# 	fname = "Gene-and_Read-Count_correlation_per_Chr_total"
# 	plot(1, type="n", main=fname,
# 		 sub = "Each color marks one ot the 23 chromosomes",
# 		 xlab="Nr of genes expressed 10+ reads", ylab="Total Read Count", xlim=c(0, 75), ylim=c(0, 12000))
# 	points(ChrGeneCount_, ChrSums_rawChrAnnot, pch =c(3,4), cex=0.5, col = chr_colz)
# 	wplot_save_this(plotname = fname, mdlink = T)
#
# 	fname = "Gene-and_Read-Count_correlation_per_Chr_LE"
# 	plot(1, type="n", main=fname,
# 		 sub = "Each color marks one ot the 23 chromosomes",
# 		 ylab="Nr of genes expressed 25+ reads", xlab="Total Read Count", ylim=c(0, 10), xlim=c(0, 1000))
# 	points(ChrSums_rawChrAnnot, ChrGeneCount_, pch =c(3,4), cex=0.5, col = chr_colz)
# 	abline (v=minXReads)
# 	wplot_save_this(plotname = fname, mdlink = F)
# 	llprint (" - Below" ,minXReads, "reads, chromosomes often have 5 or less genes.
# 				 So few genes can lead to high 'artificial' chromosome-level monoallelic observation.
# 				 With 5 genes, if all corresponding reads are only PCR-duplicates of single original molecules,
# 				 the chance of observing consistent parental bias is: 0.5^5 =  3,25%. See ij supfig in folder." )
# } # if


# ## Chromosome X-reactivation ---------------------------------------------------------------------------------------------------------------------------------------------------------
llogit (" ## Chromosome X-reactivation ")

llprint (""); llprint ("The alleic expression of chrX is determined by adding up all paternal, and all maternal reads respectively in each cell.")
llprint ("The allelic expression bias is summarized in the xPP-ratio: PatReads_chrX / ( PatReads_chrX + MatReads_chrX )")
llprint ("The significance of allelic bias is determined by comparing each cell's xPP-ratio to the distribution autosomal-PP-ratios. If they fall extreme of the range within 90, 95 or 99% of the autosomes fall, the bias is unlikely to be caused by chance, hence called siginificant.")
llprint ("The full distribution autosomal-PP-ratios is binned by expression level, so that chrX can be compared to autosomes with the same expression level")

llprint  (" We apply the following statistical test to account for expression dependent biases:")
llprint ("1. We add up all reads coming for the paternal and the maternal allele, to gain statistical power.")
llprint ("- We remove low quality samples, leaving us with", nr_HQ_cells, "cells")
llprint ("- We disregard chromosomes that has less than ",minXReads," reads in total.", nr_chrX_ExprAboveThr, "cells express chrX sufficiently highly.")
llprint ("- We slice up the remaining distribution in bins of 100 data points, from low to high expression")
llprint ("- The autosomes's total expression yield a distribution in the paternal-maternal dimension.")
llprint ("- The percentile [position] in these distributions mark the chance of a the chromosomal bias being random.")
llprint ("- If chrX is falling in the  trailing 2.5%-s of the distribution, the observaiton has <5% probability that it happened by chance, so we call the bias significant.")

# dim(PP_chr_level)

# Bin /slice by total expression level --------------------------------------------------------------------------------
autosomez_exp = sort(unlist(ChrSums[-23,]))
autosomez_exp = autosomez_exp [autosomez_exp > minXReads]

borderz = seq (0,l(autosomez_exp), by=BinSize); borderz						# Get bin sizes with 100 entries, by sorting the data!
slices = na.omit(unique(autosomez_exp[borderz]))						# Read count values of the 100th, 200th etc entry
slices[l(slices)] = max(autosomez_exp)									# Extend the upper limit to the maximal data point
slices = c(minXReads, slices); slices									# Add the lower limit

ChrSums_backup = ChrSums
ChrSums[23,fail_X_QC] =  0 										# Remove chr-X entries for genes that did not make it through the filter
sum((ChrSums[23,]==0)) # control

percentiles_of_chrX_ALL =  numeric(0)
hitzzz =  non_hitzzz = NullDistr = TestChr = percentiles_of_chrX_ = h = nh = LevelOfSiginif = list(NA)
OneSidedPC = (1-SignLevel)/2
for (s in 1:(l(slices)-1)) {
	# s=1
	bin_diag = c(slices[s],slices[s+1]);	print (bin_diag)
	# InDaBin= (ChrSums >= bin_diag[1])& (ChrSums < bin_diag[2]); sum(InDaBin)
	InDaBin= (ChrSums[ , HQ_Cellnames] >= bin_diag[1])& (ChrSums[ , HQ_Cellnames] < bin_diag[2]); sum(InDaBin)

	any_print(sum(InDaBin, na.rm = T), "data points in da bin.")
	xIDB 	= InDaBin[23, ]
	autoIDB = InDaBin[1:22,]

	NullDistr[[s]] = PP_chr_level[1:22,][autoIDB]
	LevelOfSiginif[[s]] = quantile(NullDistr[[s]], c(OneSidedPC, 1-OneSidedPC), na.rm = T)      			# Trailing quantiles given sign. level.
	TestChr[[s]]   = PP_chr_level[23,][xIDB]
	nr_tested = l(TestChr[[s]])
	if (nr_tested){
		# FullDistr = unlist(c(NullDistr[[s]], TestChr[[s]] ) )
		# percentiles_of_chrX_[[s]] = 100*ecdf(FullDistr) (TestChr[[s]]);	names(percentiles_of_chrX_[[s]]) = names (TestChr[[s]]);percentiles_of_chrX_[[s]] 		# calculate the concrete percentiles belonging to chrX-data point
		percentiles_of_chrX_[[s]] = 100*ecdf(NullDistr[[s]]) (TestChr[[s]]);	names(percentiles_of_chrX_[[s]]) = names (TestChr[[s]]);percentiles_of_chrX_[[s]] 		# calculate the concrete percentiles belonging to chrX-data point
		# ECDF is exaclty doing the following: 100*match(x= xIDB, table = names(sort(FullDistr)))/l(FullDistr)
		SignifSkewedAllelicExp = TestChr[[s]] < LevelOfSiginif[[s]][1] | TestChr[[s]] > LevelOfSiginif[[s]][2]; SignifSkewedAllelicExp
		h[[s]] =	names(SignifSkewedAllelicExp[SignifSkewedAllelicExp==T])
		nh[[s]] = names(SignifSkewedAllelicExp[SignifSkewedAllelicExp==F])
	} else { h[[s]] = nh[[s]] = percentiles_of_chrX_[[s]] = NA } # if positive length
	any_print("Nr of chrX in the plots: ", nr_tested, "; Significant:", sum(SignifSkewedAllelicExp), "; Germ:", sum(!Control[h[[s]]]))
	any_print("SignifSkewedAllelicExp: ", h[[s]]); 	# any_print("NOT SignifSkewedAllelicExp: ", nh[[s]])
} # for
hitzzz = sort(na.omit.strip(unlist (h))); hitzzz =  setdiff(hitzzz, fail_X_QC_names_germ); hitzzz
non_hitzzz =  sort(na.omit.strip(unlist (nh))); non_hitzzz =  setdiff(non_hitzzz, fail_X_QC_names_germ); non_hitzzz

GermCellNonReactivated = which_names(!Control[hitzzz])
GermCells_Reactivated = which_names(!Control[non_hitzzz])

nr_germ = sum(cell_type_annot_simplest!="SOMATIC")
nr_germ_HQ = l(GermNamez_HQ)
llprint("Out of",nr_germ,"germ cells,", nr_germ_HQ, "are HQ.",l(pass_X_QC_names_germ), "of that is passing X-QC." )
stopifnot(l(non_hitzzz) + l(hitzzz) + l(fail_X_QC_names_germ) == l(unique(c((non_hitzzz), (hitzzz), (fail_X_QC_names_germ)))))


percentiles_of_chrX_ALL = rep (NaN, l(cell_IDs)); names(percentiles_of_chrX_ALL) = cell_IDs; SignifSkewedExp_chrX = percentiles_of_chrX_ALL;# make  a full, named vectors
a = sort(unlist (percentiles_of_chrX_))
percentiles_of_chrX_ALL[names(a)] = a

SignifSkewedExp_chrX[hitzzz] = T
SignifSkewedExp_chrX[non_hitzzz] = F

ChrSums = ChrSums_backup

# "look up in proper annotation" -------------------------------------------------------------------------------------
llogit("### Results")

names(cell_type_annot) = cell_IDs

a=lookup (c("GO", "AD","MALE"), cell_type_annot[hitzzz], exact =T, report=F)
llprint ("**- ", a$ln_hits, " Somatic cells have a significantly biased X-expression**, which are: ", names(a$hits), " annotated as: " ,a$hits )
llprint ("**- ", l(a$nonhits), "Germ cells have a significantly biased X-expression**, which are: ", names(a$nonhits), " annotated as: " ,a$nonhits )
chrX_bias = names(a$nonhits);
write.simple.vec(chrX_bias, ManualName = "~/Google_Drive/X_react_Data/Abelz/Metadata2/chrX_bias.vec")

a=lookup (c("GO", "AD","MALE"), cell_type_annot[non_hitzzz], exact =T, report=F)
llprint ("**- ", a$ln_hits, " Somatic cells have an insignificantly biased X-expression**, which are: ", names(a$hits), " annotated as: " ,a$hits)
chrX_WeirdoControls = names(a$hits)
llprint("**- ", l(a$nonhits), "Germ cells have an insignificantly biased X-expression**, which are: ", names(a$nonhits), " annotated as: " ,a$nonhits )
write.simple.tsv(chrX_WeirdoControls)
# chrX_WeirdoControls = c("d2_35_gf", "d4_02_gf", "d1_31_af", "d2_18_af")

# Barplot summary -------------------------------------------------------------------------------------
tbl_all_samples = table(cell_type_annot_simple[tested])
kat = names (tbl_all_samples)
tbl_hitz = table_fixed_categories(cell_type_annot_simple[hitzzz], categories_vec = kat)
tbl_non_hitz = table_fixed_categories(cell_type_annot_simple[non_hitzzz], categories_vec = kat)

# Barplot for celltype and ratios -----
tbl_hitz / tbl_all_samples
Pc_of_Cells_w_sign_Xchr_expr_bias = 100*tbl_hitz / tbl_all_samples
wbarplot(Pc_of_Cells_w_sign_Xchr_expr_bias, mdlink = T)

Sensitivity = percentage_formatter(tbl_hitz / tbl_all_samples); names (Sensitivity) = names (tbl_all_samples)
llprint (Sensitivity["SOMATIC"], " of ", tbl_all_samples["SOMATIC"], " SOMATIC cells;\n ",
		 Sensitivity["EGC"]," of ", tbl_all_samples["EGC"],  " EGC;\n",
		 Sensitivity["LGC"], " of ", tbl_all_samples["LGC"],  "the LGC  have a significantly parentally biasef chr-X." )
# Bootstrap for barplot -------------------------------------------------------------------------------------
# "re-sample the autosomal distribution with replacement, 	do the same for X, 	use the curtoff calculation"
# http://stackoverflow.com/questions/10264025/r-sample-a-vector-with-replacement-multiple-times.

# write out metadata -------------------------------------------------------------------------------------
TestResults_X = cell_NR[HQ_samples]
TestResults_X[fail_X_QC] = "Fail QC"
TestResults_X[hitzzz] = "Sign. Different from Autosomes"
TestResults_X[non_hitzzz] = "Biallelic"
Rank_X = (2*abs((percentiles_of_chrX_ALL[HQ_Cellnames])-50))

ST2c_XchrTest =  t(rbind(
	"Annotation" = cell_type_annot[HQ_Cellnames],
	# "HQ_samples" = HQ_samples[HQ_Cellnames],
	"HQ_samples" = HQ_Cellnames,
	"DP_chrX" = ChrSums[23, HQ_Cellnames],
	"DP_pat_chrX" = ChrSums_Pat[23, HQ_Cellnames],
	"DP_mat_chrX" = ChrSums_Mat[23, HQ_Cellnames],
	"Rank_X" = Rank_X,
	"fail_X_QC" 			= fail_X_QC,
	"percentiles_of_chrX" = percentiles_of_chrX_ALL[HQ_Cellnames],
	"SignifSkewedExp_chrX" = SignifSkewedExp_chrX[HQ_Cellnames],
	"TestResults_X" 		= TestResults_X
)) # metadata_Xreact

write.simple.tsv(ST2c_XchrTest, ManualName = "~/Google_Drive/X_react_Data/Supplementary_tables/ST2c_XchrTest.tsv")
write.simple.tsv(ST2c_XchrTest, ManualName = "~/Google_Drive/X_react_Data/Abelz/Metadata2/metadata_Xreact.tsv")

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the outliers the ------------------------------------------------------------------------------------------------------------------------------------

thr = 1500
not_on_plot= sum(ChrSums_Mat[23,]>thr);not_on_plot

fname = "X_bias"
plot(1, type="n", main=fname,
	 sub = paste ("not_on_plot: ",not_on_plot),
	 xlab="Maternal Chromosomal Expression (read count)", ylab="Paternal Chromosomal Expression (read count)", xlim=c(0, thr), ylim=c(0, thr))
for (sl in slices){ 															# Add Slices
	abline (a=sl, b=-1, col ="darkgrey")
}
bins = T
if (bins) {
  slices_ = slices[-l(slices)] # remove last element
  slices_ = slices
  RC_Floors = slices_[-l(slices_)]
	RC_Ceilings =  slices_[-1]
	NrOfSlices = l(slices_)
	Limitz = sort(c(RC_Floors, RC_Ceilings) )
	xx = as.data.frame(LevelOfSiginif) #
	colnames(xx) = slices_[-1]; xx

	Xorder = c(1,2,2,1) 											# to draw the rectangle, you need the corners anti-clockwise

	for (r in 1:NrOfSlices) { print (r);
		q =r+1
		Ypz = xx[ Xorder, r]					# Get the Paternal read-ratios in the order of the rectangle
		Xpz = (1-xx)[ Xorder, r]				# Maternal = 1 -Paternal
		mplyr = slices_[c(r,r,q,q)]; mplyr;	# Multiply read-ratios with actual read counts corresponding to each '-1'-diagonals
		x= mplyr * Xpz
		y= mplyr * Ypz; x+y
		polygon (x,y, density = NULL, angle = 45, border = NA, col = rgb(1, 1, 0.5, 0.3) )
	} # for
} # if (bins)

points (ChrSums_Mat[-23,], ChrSums_Pat[-23,], pch =20, cex =0.5 )				# Add null distribution

points (ChrSums_Mat[23,ControlNamez_HQ], ChrSums_Pat[23,ControlNamez_HQ], pch =16, col=cell_type_col_simplest[ControlNamez_HQ])
# text (ChrSums_Mat[23,ControlNamez_HQ], ChrSums_Pat[23,ControlNamez_HQ], labels = ControlNamez_HQ, pos =4, cex = 0.2)
points (ChrSums_Mat[23,GermNamez_HQ], ChrSums_Pat[23,GermNamez_HQ], pch =16, col=cell_type_col_simplest[GermNamez_HQ])
# points (ChrSums_Mat[23,][GermCellNonReactivated], ChrSums_Pat[23,][GermCellNonReactivated], pch =16, col = "red")
# which(!Control[names(ChrSums_Mat[23,])])
text (ChrSums_Mat[23,][GermNamez_HQ], ChrSums_Pat[23,][GermNamez_HQ], label =GermNamez_HQ, col="darkblue",pos = 3, cex=0.2)
# text (ChrSums_Mat[23,][GermCellNonReactivated], ChrSums_Pat[23,][GermCellNonReactivated], label =cell_NR[GermCellNonReactivated], col="darkblue",pos = 3, cex=1)
x = c(0, minXReads, 0); y = c(0, 0, minXReads) 			# Draw triangle in unevaluated region
polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))

lx = c( "Autosomes", names(CellTypeColz_Simplest))
lc =  c(1, CellTypeColz_Simplest)

# lpc = c(20,16,16,16)
lcex = c(0.5,1,1,1)
legend("topright", lx, col =lc, inset = .02, pch =16, cex =0.8)

wplot_save_this(plotname = fname, mdlink = T)

I_believe = c( "d1_04_gf", "d2_12_gf", "d2_15_gf", "d2_16_gf")


# Look up gene-wise scatter plots -------------------------------------------------------------------------------------
llprint("# Look up what actual genes do in the candidate cells")

llprint("## Log2 plots look very unconvincing for the",l(GermCellNonReactivated),"candidates")
prefix =kollapse( OutDirOrig, "/Scatter_PAT_log2/PNG/SNP_scatter_log2_", print =F )
for (g in GermCellNonReactivated) {
	llogit (MarkDown_ImgLink_formatter( prefix, g, "_",cell_type_annot_simple[g],".plot.png"))
}
llogit("Except the d2_15 and d2_16 the rest is not convincing. To understand why they pop up as non-reactivated, we should look at the linear scale plots.")


llprint("## Linear plots show that single, highly expressed genes determine, which cell are called non-reactivated")
llogit("*A single genes effect is especially apparent in d1_03_gf, note the X-gene at around y=300.*")
prefix =kollapse( OutDirOrig, "/Scatter_PAT/PNG/SNP_scatter_", print =F )
for (g in GermCellNonReactivated) {
	llogit (MarkDown_ImgLink_formatter( prefix, g,".plot.png"))
}

llprint("## Look up what actual genes do in already reactivated cells")
prefix =kollapse( OutDirOrig, "/Scatter_PAT_log2/PNG/SNP_scatter_log2_", print =F )
for (g in GermCells_Reactivated) {
	llogit (MarkDown_ImgLink_formatter( prefix, g,"_",cell_type_annot_simple[g],".plot.png"))
}

llogit("----------------")
llprint("## The ", l(chrX_WeirdoControls) ,"Insignificant Controls look all  but 1 indeed very noisy")
for ( c in chrX_WeirdoControls) {
	llogit(MarkDown_ImgLink_formatter(OutDirOrig, "/Scatter_PAT_log2/PNG/SNP_scatter_log2_",c, "_SOMATIC.plot.png"))
}

ChrSums_backup = ChrSums # restore without NA-s
# LogPlot -------------------------------------------------------------------------------------

if (log2plot) {
	ChrSums_Mat_log2 = log2(ChrSums_Mat+1)
	ChrSums_Pat_log2 = log2(ChrSums_Pat+1)

	Autsomes_log2 = cbind(
		Mat = as.numeric(unlist(ChrSums_Mat_log2[-23, ])),
		Pat = as.numeric(unlist(ChrSums_Pat_log2[-23, ]))
	)
	range(Autsomes_log2)
	X_log2 = cbind(
		Mat = as.numeric(unlist(ChrSums_Mat_log2[23, ])),
		Pat = as.numeric(unlist(ChrSums_Pat_log2[23, ]))
	)
	rownames(X_log2) = colnames(ChrSums_Mat_log2)
	pname = "X-bias_in_log-space"
	ll =c(0,15)

	plot(Autsomes_log2, xlim = ll, ylim = ll, main = pname, type ="n")
	nrz = nrz_rev = NULL
	for (i in 1:l(slices)) {
		nrz[[i]] = (seq(from = 1, to= slices[i], by = ceiling(slices[i]/100)))
		nrz_rev[[i]] = rev(nrz[[i]])
	}

	for (i in 1:l(slices)) {
		lines(spline(x=log2(nrz_rev[[i]]), y =log2(nrz[[i]])), lty =3, lwd =.75)
	}
	points(Autsomes_log2, pch=20, col = rgb(0,0,0,.2))
	ccc= cell_type_col_simple[rownames(X_log2)]
	points(X_log2, col=ccc, pch =18)
	points(X_log2[hitzzz, ], col=2)
	wplot_save_this(pname, mdlink = T)
}

# DeriveThr =F
# if (DeriveThr) {	source("~/analysis/DevCell_analysis/13.h2.chrX_Reactivation_CV_Based_thr.R")}
