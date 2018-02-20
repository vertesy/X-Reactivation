######################################################################################################
# Plot all Chr X in Active - Inactive Dimension
######################################################################################################
# source ("~/analysis/DevCell_analysis/13.b2.PlotAllChrX_on1plot.R")

# Functions ----------------------------
binner <- function (vec, binsize) {
	ranks = seq (0,l(vec), by=binsize)
	borders = sort (vec)[ranks]
	borders = c(min(vec), borders)
	names (borders) = NULL
	return (borders)
}

# function to obtain the mean
library(boot)
Bmean <- function(data, indices) {
	d <- data[indices] # allows boot to select sample
	return(mean(d))
}

# Setup ----------------------------
OutDir = OutDirOrig; OutDir
OutDir = create_set_OutDir(OutDir,"/13.b2.PlotAllChrX_on1plot"); OutDir
setup_logging_markdown ("13.b2.PlotAllChrX_on1plot.R", append = F)

# Parameters ----------------------------
# NrBins = 2
NrBins = 3
GeneClass = "escapees"
# GeneClass = "X"
UseConfInterval = FALSE
UseBoxplot = F
SecondPart =  T
Plot_LowThr_log2 = 0

CategOrdered = c("SOMATIC","EGC", "L+MGC")
cell_type_annot_simplest = cell_type_annot_simple
cell_type_annot_simplest[cell_type_annot_simplest=="MGC" | cell_type_annot_simplest=="LGC" ] = 'L+MGC'
CellTypeColz_Simplest = CellTypeColz_Simple[-4]
names(CellTypeColz_Simplest)[3] = 'L+MGC'

########################################################################################################################
gene_category_simpleESC = gene_category
gene_category_simpleESC[gene_category_simpleESC == "escapees_het"]


# Infer Active Inactive ---------------------------------------------------------------------------------------------------------------------------------------------------------
xmat = DP_mat[ gene_category_simpleESC==GeneClass, HQ_Cellnames]
xpat = DP_pat[ gene_category_simpleESC==GeneClass, HQ_Cellnames]
MatXisActive = colSums(xmat) > colSums(xpat)
	colSums(xmat) + colSums(xpat) == sum(xpat) + sum(xmat) 	# QC

# convert P-M to A-InA dimension---------------------------------------------------------------------------------------------------------------------------------------------------------
DP_active = xpat; DP_inactive = xmat
DP_active[, MatXisActive] = 	xmat[, MatXisActive]
DP_inactive[, MatXisActive] = 	xpat[ ,MatXisActive]
	sum(colSums(DP_active)) +sum(colSums(DP_inactive)) == sum(xpat) + sum(xmat) 	# QC
Expressed = (DP_active > 0) | (DP_inactive > 0)

xact= (log2(unlist(DP_active)[which(Expressed)]+1))
xinact= (log2(unlist(DP_inactive)[which(Expressed)]+1))
# xact[xact == -Inf] = -1
# xinact[xinact == -Inf] = -1
	sum(xact); sum(xinact)


########################################################################################################################
dim(Expressed)
sum(Expressed)
NrGenes = colSums(Expressed)
GeneInCell = rep (names(NrGenes), times =NrGenes); table(cell_type_annot_simplest[GeneInCell])

IndexSoma = which(cell_type_annot_simple[GeneInCell] == "SOMATIC"); 	l(IndexSoma)
IndexEGC = which(cell_type_annot_simple[GeneInCell] == "EGC"); 	l(IndexEGC)
IndexLGC = which(cell_type_annot_simple[GeneInCell] == "LGC" | cell_type_annot_simple[GeneInCell] == "MGC"); 	l(IndexLGC)

xi_LGC = 	xinact[IndexLGC]
xi_EGC = 	xinact[IndexEGC]
xi_Soma = 	xinact[IndexSoma]
xa_LGC = 	xact[IndexLGC]
xa_EGC = 	xact[IndexEGC]
xa_Soma = 	xact[IndexSoma]

Exp_LGC = xi_LGC + xa_LGC
Exp_EGC = xi_EGC + xa_EGC
Exp_Soma = xi_Soma + xa_Soma


slices_LGC = 	binner(Exp_LGC, binsize = round(l(Exp_LGC)/NrBins)); 	slices_LGC
# slices_EGC = 	binner(Exp_EGC, binsize = round(l(Exp_EGC)/NrBins)); 	slices_EGC
# slices_Soma = 	binner(Exp_Soma, binsize = round(l(Exp_Soma)/NrBins)); 	slices_Soma
slicess = slices_LGC
slicess[1] = -1.1
slicess[NrBins+1] = ceiling(max(xact+xinact, na.rm = T)); names(slicess) = NULL; slicess

BinIndex_LGC = slices_LGC
# BinIndex_EGC = slices_EGC
# BinIndex_Soma = slices_Soma

hist(Exp_LGC, breaks = 20)
# hist(Exp_EGC, breaks = 20)
# hist(Exp_Soma, breaks = 20)
rowwn = paste("UpTo_", iround(slicess[-1]),sep="")
ActStdErr = InactStdErr = ActMeans = InactMeans = ActMedians = InactMedians = PointsInBins = matrix(data = NA, nrow = NrBins, ncol = 3, dimnames = list(rowwn, CategOrdered))

# BinnedData =  split(CategOrdered, x = 1:3)
BinnedData = NULL

for (i in 1:NrBins) { # i=1
	HP = slicess[i]
	LP = slicess[i+1]

	BinnedData[["SOMATIC Act"]][[i]] = IDB_xa_Soma = xa_Soma[Exp_Soma > HP & Exp_Soma <= LP];
	BinnedData[["EGC Act"]][[i]] 	= IDB_xa_EGC = xa_EGC[Exp_EGC > HP & Exp_EGC <= LP];
	BinnedData[["L+MGC Act"]][[i]] 	= IDB_xa_LGC = xa_LGC[Exp_LGC > HP & Exp_LGC <= LP];
	BinnedData[["SOMATIC Inact"]][[i]] = IDB_xi_Soma = xi_Soma[Exp_Soma > HP & Exp_Soma <= LP];
	BinnedData[["EGC Inact"]][[i]] 		= IDB_xi_EGC = xi_EGC[Exp_EGC > HP & Exp_EGC <= LP];
	BinnedData[["L+MGC Inact"]][[i]] 	= IDB_xi_LGC = xi_LGC[Exp_LGC > HP & Exp_LGC <= LP];
	}

# calculate statisitcs
for(i in 1:3) { print(names(BinnedData)[i])
	ActStdErr[, i] 	= unlist(lapply(BinnedData[[i]], sem))
	ActMeans[,i ] 	= unlist(lapply(BinnedData[[i]], mean))
	ActMedians[,i ] 	= unlist(lapply(BinnedData[[i]], median))
	PointsInBins[ ,i ]	= unlist(lapply(BinnedData[[i]], length))
}

for(i in 4:6) { print(names(BinnedData)[i])
	InactStdErr[, i-3] 	= unlist(lapply(BinnedData[[i]], sem))
	InactMeans[,i -3] 	= unlist(lapply(BinnedData[[i]], mean))
	InactMedians[,i -3] 	= unlist(lapply(BinnedData[[i]], median))
}
colSums(PointsInBins)[CategOrdered] == table(cell_type_annot_simplest[GeneInCell])[CategOrdered]

# -------------------------------------------------------------------------------------------------------------------
llogit("## Data points per bin")
MarkDown_Table_writer_DF_RowColNames(PointsInBins)
llogit("## Bin-wise allelic means")
MarkDown_Table_writer_DF_RowColNames(ActMeans)
MarkDown_Table_writer_DF_RowColNames(InactMeans)
llogit("## Bin-wise allelic medians")
MarkDown_Table_writer_DF_RowColNames(ActMedians)
MarkDown_Table_writer_DF_RowColNames(InactMedians)
llogit("Because there are so many monoallelically detected genes, hence 0 red count from the inactive chromosome, the inact median(log2(RC-Inact))  is mostly -1")


llogit("## None but 1 of the binned data is normally distributed")
ccc = CellTypeColz_Simplest[CategOrdered]
par(mfrow = c(3, NrBins))  # 3 rows and 2 columns
for(i in 1:3) { print(names(BinnedData)[i])
	names (BinnedData[[i]]) = rowwn
	lapply(BinnedData[[i]], hist, breaks = 20, xlab = "log2(RC+1)",  xlim = c(-1,10), main =names(BinnedData)[i], col=ccc[i])
}
wplot_save_this("Active_Allele_histograms_per_bin", mdlink = T)

for(i in 4:6) { print(names(BinnedData)[i])
	names (BinnedData[[i]]) = rowwn
	lapply(BinnedData[[i]], hist, breaks = 20, xlab = "log2(RC+1)",  xlim = c(-1,10), main =names(BinnedData)[i], col=ccc[i-3])
}
wplot_save_this("Inactive_Allele_histograms_per_bin", mdlink = T)
par(mfrow = c(1, 1))  # 3 rows and 2 columns
dev.off()
llogit("- This means, unless we find true distribution, the 95%CI is not calculable.")
llogit('- ["One could still approxiamte it by bootstrapping"](http://stats.stackexchange.com/questions/112829/how-do-i-calculate-confidence-intervals-for-a-non-normal-distribution)')


# Plot Alltogether  ---------------------------------------------------------------------------------------------------------------------------------------------------------
cc= CellTypeColz_Simplest[cell_type_annot_simplest[GeneInCell]]; l(cc)
# pname ="All X genes in all cell types together"
pname = kollapse("Allelic Expression of ",GeneClass, " per Cell Type")
subb = "Standard error of the mean is shown"

plotvar = cbind("log2(read count+1) per gene from the Active Chromosome" = jitter(xact, factor = 75),
	  "log2(read count+1) per gene from the Inctive Chromosome"= jitter(xinact, factor = 75))

# plot (plotvar, col =cc, pch =20, main = pname, sub =subb, cex=0.5)
xy = c(-.5, 7)
plot (plotvar, col =cc, pch =20, main = pname, sub =subb, cex=0.5, type = "n", xlim = xy, ylim = xy,
	  xlab= expression(paste("log"[2],"(read count+1) per gene from the active chr X")),
	  ylab= expression(paste("log"[2],"(read count+1) per gene from the inactive chr X"))	)
# Add MEans
	coll = rbind(ccc,ccc,ccc)
	yerr = InactStdErr # c(sem(IDB_xi_EGC), sem(IDB_xi_LGC), sem(IDB_xi_Soma)); yerr
	xerr = ActStdErr # c(sem(IDB_xa_EGC), sem(IDB_xa_LGC), sem(IDB_xa_Soma)); xerr
	# abline(a=0,b=1, lty=3, col="grey")

# shaded areas
	for (i in c(1,3,2)) {
		midline = spline(ActMeans[,i], InactMeans[,i])
		highline = spline(ActMeans[,i], (InactMeans[,i] + yerr[,i]))
		lowline = spline(ActMeans[,i], (InactMeans[,i] - yerr[,i]))

		polygon( c(rev(highline$x), lowline$x), c(rev(highline$y), lowline$y), col = ccc[[i]], border = NA, density = 100, angle = c(45,135,45)[i])
		lines(midline, lty=1, col = 1, lwd=2 )
		lines(highline, lty=2, col = "grey33", lwd=1.5 )
		lines(lowline, lty=2, col = "grey33", lwd=1.5 )
	}
# Add mean values
	arrows(ActMeans+xerr, InactMeans, ActMeans-+xerr, InactMeans, angle=90, code=3, length=0.1, lwd = 2) 		# X dim error bar
	arrows(ActMeans, InactMeans+yerr, ActMeans,InactMeans-yerr, angle=90, code=3, length=0.1, lwd = 2) 		# Y dim error bar
	pcz =  sort(rep(c(21,22,24),3))
	PointCol = matrix_from_vector (ccc, HowManyTimes = NrBins, IsItARow = F)
	points(ActMeans, InactMeans, bg = PointCol, pch =pcz, cex=1.25, lwd = 2)
	# points(ActMeans, InactMeans, bg = coll[1:NrBins,], pch =pcz, cex=1.25, lwd = 2)
# Caclualte and Add autosomes
	source("~/analysis/DevCell_analysis/13.b3.PlotAllChrX_on1plot.AddAutosomes.R")
	lines(spline(Auto_Act_mean, Auto_Inact_mean), bg = 1, pch =23, cex=0.75, lwd = 2)
	points(Auto_Act_mean, Auto_Inact_mean, bg = 1, pch =23, cex=1)

legend("topleft", legend = c(CategOrdered, "AUTOSOMES"), col =c(ccc,1), inset = .02, pch =c( 16,15,17,18), cex =1, bty = "n")

wplot_save_this(pname, mdlink = T)
llogit("Fat rhomboids: binned medians with 95% Confidence Interval (estimated by bootstrapping)")


llogit("looking  on cell categories together, all germ cells look significantly different from control cells,
	   but not from each other, hence we could call all germ cells reactivated.")
llogit ("We later falsify this observation on the single cell level")


# EGC_Act_bin3 =(BinnedData$`EGC Act`[[3]])
# EGC_Inact_bin3 =(BinnedData$`EGC Inact`[[3]])
# whist(EGC_Act_bin3)
# -------------------------------------------------------------------------------------------------------------------
if (SecondPart) {

	# -------------------------------------------------------------------------------------------------------------------
	# Bootstrapping CI --------------------------------------------------------------------------------------------------------
	if (UseConfInterval) {

		BootstrappedCIs_I = NULL
		colln =  c ("LowerBoundary", "UpperBoundary")
		for (name in names(BinnedData)) {
			print(name)
			for (i in 1:NrBins) {
				# bootstrapping with 1000 replications
				results <- boot(data=BinnedData[[name]][[i]], statistic=Bmean, R=10000)
				aa = boot.ci(results, type=c("norm", "basic", "perc", "bca"))
				# 		BootstrappedCIs_I[[name]][[i]] = rbind( 	"norm" = (aa$norm[2:3])
				# 										 , "basic" = (aa$basic[4:5])
				# 										 , "perc" = (aa$perc[4:5])
				# 										 , "bca" = (aa$bca[4:5]))
				# colnames(BootstrappedCIs_I[[name]][[i]]) = colln; BootstrappedCIs_I[[name]][[i]]
				BootstrappedCIs_I[[name]][[i]] = aa$basic[4:5]
				names(BootstrappedCIs_I[[name]][[i]]) = colln; BootstrappedCIs_I[[name]][[i]]
				print(aa$basic[4:5])
				print('------')
			}
		}

		subb = "Bootstrapping based 95% CI around binned mean is shown"
		plot (plotvar, col =cc, pch =20, main = pname, sub =subb, cex=0.5)
		# Add booCI
			ccc = CellTypeColz_Simplest[CategOrdered]
			coll = rbind(ccc,ccc,ccc)
			xerrHi = BootstrappedCIs_Hi[, 1:3]
			xerrLo = BootstrappedCIs_Low[, 1:3]
			yerrHi = BootstrappedCIs_Hi[, 4:6]
			yerrLo = BootstrappedCIs_Low[, 4:6]

		arrows(xerrHi, InactMeans, xerrLo, InactMeans, angle=90, code=3, length=0.1) 		# X dim error bar
		arrows(ActMeans, yerrHi, ActMeans, yerrLo, angle=90, code=3, length=0.1) 		# Y dim error bar
		points(ActMeans, InactMeans, bg = coll, pch =23, cex=1.25)
		legend("topleft", legend = CategOrdered, col =ccc, inset = .02, pch =16, cex =0.8, bty = "n")

		wplot_save_this(kollapse(pname, "_CI", print = F), mdlink = F)
	} # if (UseConfInterval)

	# Plot Soma  ---------------------------------------------------------------------------------------------------------------------------------------------------------
	pname ="All X Transcripts in all somatic cells"

	plot (xact[IndexSoma],xinact[IndexSoma], col =CellTypeColz_Simple["SOMATIC"], pch =c(3,4), main = pname,
		  xlab= expression(paste("log"[2],"(read count+1) per gene from the active chr X")),
		  ylab= expression(paste("log"[2],"(read count+1) per gene from the inactive chr X")))
	wplot_save_this(pname, mdlink = T)
	llogit("You see that far the most genes come from the Active allele in somatic cells")

	# Plot Germ  ---------------------------------------------------------------------------------------------------------------------------------------------------------
	pname ="All X Transcripts in all EGCs"

	plot (xact[IndexEGC],xinact[IndexEGC], col =CellTypeColz_Simple["EGC"], pch =c(3,4), main = pname,
		  xlab= expression(paste("log"[2],"(read count+1) per gene from the active chr X")),
		  ylab= expression(paste("log"[2],"(read count+1) per gene from the inactive chr X")))
	wplot_save_this(pname, mdlink = T)
	llogit("The cloud int he middle is a bit shifted towards the acitve allele,
		   and (although saturted) there are more dos on the Active-axis")

	# Plot Late Germ  ---------------------------------------------------------------------------------------------------------------------------------------------------------
	pname ="All X Transcripts in all LGCs"

	plot (xact[IndexLGC],xinact[IndexLGC], col =CellTypeColz_Simple["LGC"], pch =c(3,4), main = pname,
		  xlab= expression(paste("log"[2],"(read count+1) per gene from the active chr X")),
		  ylab= expression(paste("log"[2],"(read count+1) per gene from the inactive chr X")))
	wplot_save_this(pname, mdlink = T)
	llogit("The proportion of biallelic genes seem the most in LGCs")

	########################################################################################################################
	# ---------------------------------------------------------------------------------------------------------------------------------------------------------

	# llogit("Since in the LE bins, are almost only 0-s and 1-s, the median will be one of these values. It is not really meaningful.")
	# llogit("That is why the median plot starts with 0-s.")
	# llogit("The Mean-based plot is more meaningful, showing more subtle differences.")

	#######################################################################################################################
	# boxplots per category  ---------------------------------------------------------------------------------------------------------------------------------------------------------
	if (UseBoxplot) {
		PActive_matrix = DP_active / (DP_active + DP_inactive)
		PActive_ALL =PP_matrix
		PActive_ALL[,MatXisActive] = 1- PActive_ALL[,MatXisActive]

		Categz = c("Autosomes", sort(names(table(cell_type_annot_simple))));Categz
		# Categz = c(sort(names(table(cell_type_annot_simple))));Categz
		GenePP_per_celltype = list (NA)
		for (i in 2:l(Categz)) {
			type = Categz[i]; print (type)
			GenePP_per_celltype[[type]] = 	as.vector(na.omit(unlist(PActive_matrix[which_names(gene_category =="X"), which_names(cell_type_annot_simple[HQ_Cellnames] == type)])))
			# names(GenePP_per_celltype)[[i]] = type
		}
		GenePP_per_celltype[[1]] = as.vector(na.omit(unlist(PActive_ALL[which_names(gene_category =="autosomal") , ])))
		names(GenePP_per_celltype)[[1]] = "Autosomes"

		str(GenePP_per_celltype)
		boxplot(GenePP_per_celltype,
				ylab = "% of reads from the active allele", main = "Active Allelic Expression", las=2)
		stripchart(GenePP_per_celltype, vertical = TRUE, method = "jitter",
				   pch = 16, add = TRUE, cex=1.5
				   , col = c("black",CellTypeColz[c("EGC", "LGC", "MGC", "AD" )])
				   # , col = c(CellTypeColz[c("EGC", "LGC", "MGC", "AD" )])
		)

		pname =("Boxplot_Active_Allelic_expression_per_gene")
		wplot_save_this(pname, mdlink = T)
		llogit("I think what we see as an active bias is technical and comes from the way we inferred the active allele")
	} # if (UseBoxplot)

	# # BOTH boxplots per category  ---------------------------------------------------------------------------------------------------------------------------------------------------------
	llogit("------")
	llogit("## Counting biallelically expressed genes show the differce between somatic and all germ cells")
	Categz = c("Autosomes", sort(names(table(cell_type_annot_simple))));Categz
	# Categz = c(sort(names(table(cell_type_annot_simple))));Categz
	PercetageOfBiallelcGenes = numeric(4)
	for (i in 2:l(Categz)) {
		type = Categz[i]; print (type)
		MAexpList = as.vector(na.omit(unlist(MonAllExp_pat[which_names(gene_category =="X"), which_names(cell_type_annot_simple[HQ_Cellnames] == type)])))
		PercetageOfBiallelcGenes[i] = 	pc_in_total_of_match(MAexpList, "Both")
	}

	PercetageOfBiallelcGenes[1] = 	pc_in_total_of_match(as.vector(na.omit(unlist(MonAllExp_pat[which_names(gene_category =="autosomal") , ]))), "Both")
	PercetageOfBiallelcGenes= 100*PercetageOfBiallelcGenes
	names(PercetageOfBiallelcGenes) = Categz
	wbarplot(PercetageOfBiallelcGenes, mdlink = T, ylab = "% of genes expressed biallelically")

} # if (SecondPart)
