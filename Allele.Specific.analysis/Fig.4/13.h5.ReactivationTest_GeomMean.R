######################################################################################################
# Geom mean based test of reactivaiton
######################################################################################################
# source ("~/analysis/DevCell_analysis/13.h5.ReactivationTest_GeomMean.R")

# Setup ----------------------------
OutDir = OutDirOrig; OutDir
OutDir = create_set_OutDir(OutDir,"/13.h5.ReactivationTest_GeomMean.R"); OutDir
setup_logging_markdown ("13.h5.ReactivationTest_GeomMean.R", append = F)

# Parameters ----------------------------
minXGMean = 2  # determined on the fly, put back here
minGenes = 3
SignLevel = 0.99
# SignLevel = 0.95
	SelectedCandidates = paste0("GM_",100*SignLevel)
NrBins =  5
log_settings_MarkDown(minXGMean, minGenes, SignLevel, NrBins)

# ----------------------------------------------------------------------
DP_m = DP_mat_clean[,which_names(HQ_samples)]							# HQ samples only
DP_p = DP_pat_clean[,which_names(HQ_samples)]
DP_all = DP_p + DP_m


ChrGeomMean_Mat = aggregate(DP_m, by=list(chrz_clean),  FUN=gm_mean, na.rm=TRUE);
ChrGeomMean_Mat=	data.matrix(ChrGeomMean_Mat[,-1])
ChrGeomMean_Pat = aggregate(DP_p, by=list(chrz_clean),  FUN=gm_mean, na.rm=TRUE);
ChrGeomMean_Pat=	data.matrix(ChrGeomMean_Pat[,-1])
dim(ChrGeomMean_Pat)
"ChrGeomMean_Total cannot be calculated as gm_mean(sum), it must be calculated as sum(allelic gm_means). See script: zz.WhyNotBinByTotal.R"
# ChrGeomMean_Total = aggregate(DP_all, by=list(chrz_clean),  FUN=gm_mean, na.rm=TRUE);
# ChrGeomMean_Total=	data.matrix(ChrGeomMean_Total[,2:59])
ChrGeomMean_Total = ChrGeomMean_Mat + ChrGeomMean_Pat
rownames (ChrGeomMean_Total) = rownames (ChrGeomMean_Mat) = rownames (ChrGeomMean_Pat) = chromosomes

MaxGM = max(c(ChrGeomMean_Mat[23,], ChrGeomMean_Pat[23,]), na.rm = T)
write.simple.tsv(ChrGeomMean_Mat)
write.simple.tsv(ChrGeomMean_Pat)
write.simple.tsv(ChrGeomMean_Total)

# ------------------------------------------------------------------------------------------------------------------------------------------------
# QC:  Read Count based ordering ----------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
llogit (" ## Finding the best threshold of reliability for Chromosome X-reactivation test")
ChrSumX = colSums(DP_SUM_parental_HQ[gene_category =="X",ControlNamez_HQ[(donor[ControlNamez_HQ] != "D3")]])
FemaleSoma = names(ChrSumX)
ChrSumX = ChrSumX[ChrSumX>0]
PP_GeomMeanChr =  ((ChrGeomMean_Pat) / (ChrGeomMean_Mat + ChrGeomMean_Pat))

old = F
# if (old) {
# 	# PP_GeomMeanChrX = ((ChrGeomMean_Pat-0.9999) / (ChrGeomMean_Mat + ChrGeomMean_Pat -1.9999))[23,]
# 	# PP_GeomMeanChrX = PP_GeomMeanChrX *200
# 	# PP_GeomMeanChrX[PP_GeomMeanChrX<.5] = 1-PP_GeomMeanChrX[PP_GeomMeanChrX<.5];
# } else if (!old) {
	PP_GeomMeanChrX = PP_GeomMeanChr [23, names(ChrSumX) ];PP_GeomMeanChrX;PP_GeomMeanChrX
	PP_GeomMeanChrX[PP_GeomMeanChrX<.5] = 1-PP_GeomMeanChrX[PP_GeomMeanChrX<.5];
	PP_GeomMeanChrX = PP_GeomMeanChrX - min(PP_GeomMeanChrX) 						# Scale to 0-100
	Scaler = 100/max(PP_GeomMeanChrX)
	PP_GeomMeanChrX = Scaler * PP_GeomMeanChrX
# }

hist(PP_GeomMeanChrX, breaks = 20)

ord = names(sort(ChrSumX))
PP_GeomMeanChrX = PP_GeomMeanChrX[ord]

llprint ("### Surprisingly for chrX there is no clear correlation between reliability and expression level, **unlike the clear correlation observed in imprinted genes!**: ")
pname= "ChrX Allelic Expr. bias in Controls vs ReadCount GMean"
wbarplot(PP_GeomMeanChrX-min(PP_GeomMeanChrX), ylab="Relative Allelic Expression bias", col="olivedrab3",
		 plotname = pname, main = "Somatic Allelic Bias of GMean, ordered by RC ", sub = "Read Counts are labeled within the bars")
a=barplot(PP_GeomMeanChrX, plot=F)
text (x=a, y=10, labels = sort(ChrSumX), srt=90, cex=1)
pname = kollapse(pname, ".barplot")
wplot_save_this (plotname = pname, mdlink = T, ManualName = T)

# ------------------------------------------------------------------------------------------------------------------------------------------------
# Geom Mean based ordering ----------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
TotalGmeanX = ChrGeomMean_Pat[23,] + ChrGeomMean_Mat[23,]

ord= names(sort(TotalGmeanX[FemaleSoma]))
PP_GeomMeanChrX = PP_GeomMeanChrX[ord]

llprint ("## Not much correlation between reliability and the gm_mean of the total read count: ")
pname= "ChrX Allelic Expr. bias in Controls vs Total Geom Mean"
wbarplot(PP_GeomMeanChrX-min(PP_GeomMeanChrX), ylab="Relative Allelic Expression bias", col="olivedrab3",
		 plotname = pname, main = "Somatic Allelic Bias of GMean, ordered by gm_mean ", sub = "Read Counts are labeled within the bars")
a=barplot(PP_GeomMeanChrX, plot=F)
text (x=a, y=10, labels = sort(iround(TotalGmeanX[ord])), srt=90, cex=1)
pname = kollapse(pname, ".barplot")
wplot_save_this (plotname = pname, mdlink = T, ManualName = T)

# ------------------------------------------------------------------------------------------------------------------------------------------------
# GeneCountX based ordering ----------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------
GeneCountX = colSums(DP_SUM_parental_HQ[gene_category =="X",ControlNamez_HQ[(donor[ControlNamez_HQ] != "D3")]] >0)
GeneCountX = GeneCountX[GeneCountX>0]

ord = names(sort(GeneCountX[FemaleSoma]))
PP_GeomMeanChrX = PP_GeomMeanChrX[ord]


llprint ("## The minimum of 3 genes can serve as a reasonable threshold")
pname= "ChrX Allelic Expr. bias in Controls vs GENE COUNT GMean"
wbarplot(PP_GeomMeanChrX-min(PP_GeomMeanChrX), ylab="Relative Allelic Expression bias", col="olivedrab3",
		 plotname = pname, main = "Somatic Allelic Bias of GMean, ordered by #Genes ", sub = "Read Counts are labeled within the bars")
a=barplot(PP_GeomMeanChrX, plot=F)
text (x=a, y=10, labels = GeneCountX[ord], srt=90, cex=1)
pname = kollapse(pname, ".barplot")
wplot_save_this (plotname = pname, mdlink = T, ManualName = T)



# Minimum Gene Count Filter

GeneCountX = colSums(DP_SUM_parental_HQ[gene_category =="X",]>0)
fail_X_QC = (GeneCountX <= minGenes)
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

llogit("The minimum gene count filter is not applied on the autosomes, but it is irrelevant for them: only 5 times happen that an autosome expresses <=", minGenes,"genes out of 1276 cases.")
# for (c in chromosomes) { any_print(c,"-", sum(colSums(DP_SUM_parental_HQ[chrz ==c,]>0)<=3))} # chr 21 is LE

# ------------------------------------------------------------------------------------------------------------------------
llogit (" ## Chromosome X-reactivation ")

# Bin /slice by total expression level --------------------------------------------------------------------------------
autosomez_exp = sort(unlist(ChrGeomMean_Total[-23,]))
autosomez_exp = autosomez_exp [autosomez_exp > minXGMean];


# binsize = 200
binsize = floor(l(autosomez_exp)/NrBins); binsize
borderz = seq (0,l(autosomez_exp), by=binsize); borderz						# Get bin sizes with 100 entries, by sorting the data!
slices = na.omit(unique(autosomez_exp[borderz]))						# Read count values of the 100th, 200th etc entry
slices[l(slices)] = max(autosomez_exp)									# Extend the upper limit to the maximal data point
slices = c(minXGMean, slices); slices							# Add the lower limit


percentiles_of_chrX_ALL = hitzzz =  non_hitzzz =  numeric(0)
hitzzz =  non_hitzzz = NullDistr = TestChr = percentiles_of_chrX_ = h = nh = LevelOfSiginif = LevelOfSiginif_plot = list(NA)
OneSidedPC = (1-SignLevel)/2
for (s in 1:(l(slices)-1)) {
	# s=1
	bin_diag = c(slices[s],slices[s+1]);	any_print ("Slice",s,", range:",iround(bin_diag))
	InDaBin= (ChrGeomMean_Total[ ,HQ_Cellnames] >= bin_diag[1])& (ChrGeomMean_Total[ ,HQ_Cellnames] < bin_diag[2])
	xIDB 	= which_names(InDaBin[23,]);xIDB
	autoIDB = InDaBin[1:22,]
	NrDataPoints = sum(autoIDB, na.rm = T)
	any_print(NrDataPoints, "data points in da bin.")

	NullDistr[[s]] = PP_GeomMeanChr[1:22,][autoIDB]; range(NullDistr[[s]])
	LevelOfSiginif[[s]] = quantile(NullDistr[[s]], c(OneSidedPC, 1-OneSidedPC), na.rm = T)      			# Highly Siginif;
	# FirstDataPoint = round(SignLevel* NrDataPoints /2); FirstDataPoint
	# LastDataPoint = NrDataPoints- FirstDataPoint ; LastDataPoint
	# LevelOfSiginif[[s]] = sort(NullDistr[[s]], na.last = NA)[c(FirstDataPoint,LastDataPoint)]; LevelOfSiginif[[s]]
	TestChr[[s]]   = PP_GeomMeanChr[23,][xIDB]
	nr_tested = l(TestChr[[s]])
	if (nr_tested){
# 		FullDistr = unlist(c(NullDistr[[s]], TestChr[[s]] ) )
# 		percentiles_of_chrX_[[s]] = 100*ecdf(FullDistr) (TestChr[[s]]);	names(percentiles_of_chrX_[[s]]) = names (TestChr[[s]]);percentiles_of_chrX_[[s]] 		# calculate the concrete percentiles belonging to chrX-data point
		percentiles_of_chrX_[[s]] = 100*ecdf(NullDistr[[s]]) (TestChr[[s]]);	names(percentiles_of_chrX_[[s]]) = names (TestChr[[s]]);percentiles_of_chrX_[[s]] 		# calculate the concrete percentiles belonging to chrX-data point
		# ECDF is exaclty doing the following: 100*match(x= xIDB, table = names(sort(FullDistr)))/l(FullDistr)
		SignifSkewedAllelicExp = TestChr[[s]] < LevelOfSiginif[[s]][1] | TestChr[[s]] > LevelOfSiginif[[s]][2]; SignifSkewedAllelicExp
		h[[s]] =	names(SignifSkewedAllelicExp[SignifSkewedAllelicExp==T])
		nh[[s]] = names(SignifSkewedAllelicExp[SignifSkewedAllelicExp==F])
	} else { h[[s]] = nh[[s]] = percentiles_of_chrX_[[s]] = NA } # if positive length
	any_print("Nr of chrX in the plots: ", nr_tested, "; Significant:", sum(SignifSkewedAllelicExp), "; Germ:", sum(!Control[h[[s]]]))
	any_print("SignifSkewedAllelicExp: ", h[[s]]); 	 any_print("NOT SignifSkewedAllelicExp: ", nh[[s]])
} # for

hitzzz = sort(na.omit.strip(unlist (h))); hitzzz =  setdiff(hitzzz, fail_X_QC_names_germ); hitzzz
non_hitzzz =  sort(na.omit.strip(unlist (nh))); non_hitzzz =  setdiff(non_hitzzz, fail_X_QC_names_germ); non_hitzzz

nr_germ = sum(cell_type_annot_simplest!="SOMATIC")
nr_germ_HQ = l(GermNamez_HQ)
llprint("Out of",nr_germ,"germ cells,", nr_germ_HQ, "are HQ.",l(pass_X_QC_names_germ), "of that is passing X-QC." )
	stopifnot(l(non_hitzzz) + l(hitzzz) + l(fail_X_QC_names_germ) == l(unique(c((non_hitzzz), (hitzzz), (fail_X_QC_names_germ)))))

percentiles_of_chrX_ALL = rep (NaN, l(cell_IDs)); names(percentiles_of_chrX_ALL) = cell_IDs; SignifSkewedExp_chrX = percentiles_of_chrX_ALL;# make  a full, named vectors
a = sort(unlist (percentiles_of_chrX_))
percentiles_of_chrX_ALL[names(a)] =a

hist(percentiles_of_chrX_ALL, breaks = 25)

SignifSkewedExp_chrX[hitzzz] = T
SignifSkewedExp_chrX[non_hitzzz] = F
# ChrGeomMean_Total == ChrGeomMean_bac

llogit ("### Results")
a=lookup (needle = c("GO", "AD","MALE"), haystack = cell_type_annot[hitzzz], exact =T, report=F)
llprint ("**- ", a$ln_hits, " Somatic cells have a significantly biased X-expression**, which are: ", names(a$hits), " annotated as: " ,a$hits )
llprint ("**- ", l(a$nonhits), "Germ cells have a significantly biased X-expression**, which are: ", names(a$nonhits), " annotated as: " ,a$nonhits )
chrX_bias = names(a$nonhits);

a=lookup (c("GO", "AD","MALE"), cell_type_annot[non_hitzzz], exact =T, report=F)
llprint ("**- ", a$ln_hits, " Somatic cells have an insignificantly biased X-expression**, which are: ", names(a$hits), " annotated as: " ,a$hits)
chrX_WeirdoControls = names(a$hits)
llprint("**- ", l(a$nonhits), "Germ cells have an insignificantly biased X-expression**, which are: ", names(a$nonhits), " annotated as: " ,a$nonhits )

GermCellNonReactivated = which_names(!Control[hitzzz])
llprint(l(GermCellNonReactivated), "Germ cells inactive at", SignLevel, "siginificance level:", sort(GermCellNonReactivated))

# -------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the outliers the ------------------------------------------------------------------------------------------------------------------------------------
ChrGeomMean_Mat_m1 = ChrGeomMean_Mat
ChrGeomMean_Pat_m1 = ChrGeomMean_Pat

LowThr = 1
thr  = (max(c(ChrGeomMean_Pat_m1[23, ],ChrGeomMean_Mat_m1[23, ]), na.rm = T))
thr = max(cbind(ChrGeomMean_Pat_m1,ChrGeomMean_Mat_m1), na.rm = T)
not_on_plot= sum(ChrGeomMean_Mat_m1[23,]>thr);not_on_plot

fname = "Parental_Bias_of_ChrX_GM_Means"
plot(1, type="n", main=fname,
	 sub = paste ("not_on_plot: ",not_on_plot),
	 xlab="G.Mean X-Read Count on Maternal Chromosome", ylab="G.Mean X-Read Count on Paternal Chromosome", xlim=c(LowThr, thr), ylim=c(LowThr, thr))

for (sl in slices){ 															# Add Slices
	abline (a=sl, b=-1, col ="darkgrey")
}
bins = T

if (bins) {
	slices_ = slices;
	slices_[1]=2
	Ratio = t(as.data.frame(LevelOfSiginif)); 	rownames(Ratio) = paste("UpTo",iround(slices_[-1]), sep=""); Ratio
	# to draw the rectangle, you need the corners anti-clockwise
	RatioOrder = c(1,1,2,2) 			# The lower and upper signif level, chosen along the polygon in an anti-clockwise order
	SumOrder = c(1,2,2,1) 				# grey -1 diagonal line is the sum
	for (r in 1:l(slices)-1) { print (r);
		SumOrder_x = SumOrder + (r-1)
		Ratt = Ratio[ r, RatioOrder]; Ratt
		Summ = slices_[SumOrder_x]; Summ
		Ypz = Ratt * Summ;
		Xpz = Summ - Ypz
		polygon (Xpz,Ypz, density = NULL, angle = 45, border = NA, col = rgb(1, 1, 0.5, 0.3) )
	} # for
} # if (bins)


points (ChrGeomMean_Mat_m1[-23,], ChrGeomMean_Pat_m1[-23,], pch =20, cex =0.5 )				# Add null distribution
points (ChrGeomMean_Mat_m1[23,ControlNamez_HQ], ChrGeomMean_Pat_m1[23,ControlNamez_HQ], pch =16, col=cell_type_col_simple[ControlNamez_HQ])
# text (ChrGeomMean_Mat_m1[23,ControlNamez_HQ], ChrGeomMean_Pat_m1[23,ControlNamez_HQ], labels = ControlNamez_HQ, pos =4, cex = 0.2)
points (ChrGeomMean_Mat_m1[23,GermNamez_HQ], ChrGeomMean_Pat_m1[23,GermNamez_HQ], pch =16, col=cell_type_col_simple[GermNamez_HQ])
# points (ChrGeomMean_Mat_m1[23,][GermCellNonReactivated], ChrGeomMean_Pat_m1[23,][GermCellNonReactivated], pch =16, col = "red")
# which(!Control[names(ChrGeomMean_Mat_m1[23,])])
text (ChrGeomMean_Mat_m1[23,][GermNamez_HQ], ChrGeomMean_Pat_m1[23,][GermNamez_HQ], label =GermNamez_HQ, col="darkblue",pos = 3, cex=0.2)
text (ChrGeomMean_Mat_m1[23,][GermCellNonReactivated], ChrGeomMean_Pat_m1[23,][GermCellNonReactivated], label =cell_NR[GermCellNonReactivated], col="darkblue",pos = 3, cex=1)

x = c(1, minXGMean-1, 1); y = c(1,1, minXGMean-1) 			# Draw triangle in unevaluated region
polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))

lx = c( "Autosomes", names(CellTypeColz_Simple))				# LEGEND
lc =  c(1, CellTypeColz_Simple)
legend("topleft", lx, col =lc, inset = .02, pch =16, cex =0.8, bty = "n")

fn = kollapse(fname,"_Sign_",SignLevel,"_BinBy_",binsize)
wplot_save_this(plotname = fn, mdlink = T)


# Look up gene-wise scatter plots -------------------------------------------------------------------------------------
llprint("## Look up what actual genes do in the candidate cells")
GermCellNonReactivated = sort(which_names(!Control[hitzzz]))

prefix =kollapse( OutDirOrig, "/Scatter_PAT_log2/PNG/SNP_scatter_log2_", print=F)
for (g in GermCellNonReactivated) {
	llogit (MarkDown_ImgLink_formatter( prefix, g,"_",cell_type_annot_simple[g],".plot.png"))
}

llprint("### Linear plots support roughly the same conclusion")
llogit("*d1_06_gf might look like reactivated, but the actual nr of genes on the maternal axis is 6, vs only 3 expressed above 2 reads.
	   As the genes have the same nr of reads, it is hard to see, but their distinctness still visible if you zoom in (from the distorted shape of the square)*")

prefix =kollapse( OutDirOrig, "/Scatter_PAT/PNG/SNP_scatter_", print=F)
for (g in GermCellNonReactivated) {
	llogit (MarkDown_ImgLink_formatter( prefix, g,".plot.png"))
}

llprint("## Look up what actual genes do in Reactivated cells")
GermCells_Reactivated = sort(which_names(!Control[non_hitzzz]))

prefix =kollapse( OutDirOrig, "/Scatter_PAT_log2/PNG/SNP_scatter_log2_", print=F)
for (g in GermCells_Reactivated) {
	llogit (MarkDown_ImgLink_formatter( prefix, g,"_",cell_type_annot_simple[g],".plot.png"))
}
llogit("d4_08_gf is the only cell that I think might be wrongly classified as a reactivated ")

sort(PP_matrix[gene_category == "X","d1_06_gf"])

p95=F
if (p95){
	Suspicous = c("d1_01_gf","d1_06_gf", "d1_27_gf")
	llprint ("From the plots above, it seems to me that there are",l(Suspicous),"cells which might be falsly identified as non-reactivated **at 95% significance level**.
				 These are:", Suspicous, "No Suspicous cell is observed at 99% significance level.")

	llprint("## Linear plots of allelic gene expression explain some of the suspicous cells")
	prefix =kollapse( OutDirOrig, "/Scatter_PAT/PNG/SNP_scatter_", print=F)
	for (g in Suspicous) {
		llogit (MarkDown_ImgLink_formatter( prefix, g,".plot.png"))
	}
	llogit ("d1_01_gf seems clealry biased on the linear scale (~6 mat-MA VS ~2 pat-MA genes), also d1_06_gf (~5 mat-MA VS ~2 pat-MA genes)")
}


# Write out test results -------------------------------------------------------------------------------------
ChrGeomMean_Total[23,fail_X_QC] = ChrGeomMean_Mat[23,fail_X_QC] = ChrGeomMean_Pat[23,fail_X_QC] =  NaN 										# Remove chr-X entries for genes that did not make it through the filter

PP_ChrGeomMean_Total = (ChrGeomMean_Pat /(ChrGeomMean_Mat + ChrGeomMean_Pat)); dim(PP_ChrGeomMean_Total)
sum(is.nan(PP_ChrGeomMean_Total))


# TestScore_SumOfLog2 = ChrGeomMean_Total[23, ]
# TestScore_SumOfLog2 = abs(TestScore_SumOfLog2-50)*2
# TestScore_SumOfLog2 = as.vector(TestScore_SumOfLog2); names (TestScore_SumOfLog2) = names(ChrGeomMean_Total[23, ])
# TestScore_SumOfLog2 = TestScore_SumOfLog2[!fail_X_QC]
# write.simple.vec(TestScore_SumOfLog2)

# write out metadata -------------------------------------------------------------------------------------
#
TestResults_X = cell_NR[HQ_Cellnames]
TestResults_X[fail_X_QC] = "Fail QC"
TestResults_X[hitzzz] = "Sign. Different from Autosomes"
TestResults_X[non_hitzzz] = "Biallelic"
table(TestResults_X)

a = ChrGeomMean_Pat[23, HQ_Cellnames]
b = ChrGeomMean_Mat[23, HQ_Cellnames]
PP_GMean_minus1_chr_X = (a-1) / (a+b-2)

ST2c_XchrTest =  t(rbind(
	"Annotation" = cell_type_annot[HQ_Cellnames],
	"Cellnames" = HQ_Cellnames,
	"ChrGeomMean_Total" =   ChrGeomMean_Total[23, HQ_Cellnames],
	"ChrGeomMean_Pat" = ChrGeomMean_Pat[23, HQ_Cellnames],
	"ChrGeomMean_Mat" = ChrGeomMean_Mat[23, HQ_Cellnames],
	"PP_GMean_minus1_chr_X" = PP_GMean_minus1_chr_X[HQ_Cellnames],
	"fail_X_QC" 			= fail_X_QC[HQ_Cellnames],
	"percentiles_of_chrX" = percentiles_of_chrX_ALL[HQ_Cellnames],
	"SignifSkewedExp_chrX" = SignifSkewedExp_chrX[HQ_Cellnames],
	"TestResults_X" 		= TestResults_X
)) # metadata_Xreact


fname = kollapse("~/Google_Drive/X_react_Data/Supplementary_tables/ST2c_XchrTest_",SelectedCandidates, ".tsv")
write.simple.tsv(ST2c_XchrTest, ManualName = fname)

fname = kollapse("~/Google_Drive/X_react_Data/Abelz/Metadata2/metadata_Xreact_",SelectedCandidates, ".tsv")
write.simple.tsv(ST2c_XchrTest, ManualName = fname)


# chrX_bias = hitzzz[hitzzz %in% GermNamez_HQ]
fname = kollapse("~/Google_Drive/X_react_Data/Abelz/Metadata2/chrX_bias_",SelectedCandidates,".vec")
write.simple.vec(chrX_bias, ManualName =fname )


