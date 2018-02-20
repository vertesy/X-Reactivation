######################################################################################################
# 03b.Parental_scatter_plots
######################################################################################################
# source("~/analysis/DevCell_analysis/13b.Parental_scatter_plots.R")
# '''
# ### 13b.Parental_scatter_plots.R
# - Plot Pat vs Mat depth of SNPs in linear and log space
# - Annotate gene categories and gene names
# - Calculate sample quality statisitcs: ABR, MRR
# '''
# call monoallelic
try(dev.off(), silent =T)

# Setup ------------------------------------------------------------------------------------------------
UseNewInferredAnnotation = F
HiThr = 400

# Which plots to generate ------------------------------------------------------------------------------------------------
linear_plot = F
log_plot = T
LIN_impr_plot = F
log_impr_plot = T
SupplementaryFig2A = T 	# all cells for chr X
SupplementaryFig3A = T
AllImpr_on1_Plot = F
AllX_on1_Plot = F
GuideLines95 = T		# or draw manually a diagonal crossing c(log2(2), log2(2/20)) = 1, -3.322

# parameters ------------------------------------------------------------------------
DotSizeIndividualFigs = 1
DotSizeCombinedSupplFig = .5
PrintGeneNames = T
labelSize = .2 # param for AllImpr_on1_Plot
	DonorColor = T
	CellTypeColor = F
	GermOnly = F

# ------------------------------------------------------------------------
# NewNames = read.simple.vec("~/Google_Drive/X_react_Data/Abelz/Metadata2/NameVectorNew.vec");
# names(NewNames) =  cell_IDs

# Setup Imprinting annotation ----------------------------------------------------
if( UseNewInferredAnnotation ) {
	# gene_category_ImprUpDate = read.simple("~/Google_Drive/X_react_Data/Abelz/Metadata2/gene_category_ImprUpDate.vec")
	# names(gene_category_ImprUpDate) = names(gene_category)
	gene_category_used =  gene_category_ImprUpDate;
	OutDir_Impr = "/Scatter_PAT_Impr"
	OutDir_Impr_log2 = "/Scatter_PAT_Impr_log2"
	ImprScatterPath = kollapse( OutDirOrig, "/Scatter_PAT_Impr_log2_old/PNG/SNP_scatter_log2_", print=F)
} else { 						 gene_category_used =  gene_category
	OutDir_Impr = "/Scatter_PAT_Impr_old"
	OutDir_Impr_log2 = "/Scatter_PAT_Impr_log2_old"
	ImprScatterPath = kollapse( OutDirOrig, "/Scatter_PAT_Impr_log2_old/PNG/SNP_scatter_log2_", print=F)
}

# linear plots ------------------------------------------------------------------------------------------------------------------------
HiThr = 400
if (linear_plot) {
	OutDir = kollapse (OutDirOrig,"/Scatter_PAT"); if ( !exists (OutDir) ) {dir.create(OutDir)}
	for (j in 1:nr_cells) {
			# j =6
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		Maternal_Depth = (unlist(df_ReadCount[,mat_colz[j]]))
		Paternal_Depth = (unlist(df_ReadCount[,pat_colz[j]]))
		DP = cbind (Maternal_Depth, Paternal_Depth)
		MAE = table(MonAllExp_pat[,cell])
		ABR = percentage_formatter(pc_in_total_of_match(MAE, 'Both'))
		MRR = percentage_formatter(MAE['REF']/(MAE['REF']+MAE['ALT']))
		nr_SNPs_plotted = sum(rowSums(DP)>0)

		plot (c(0,HiThr), c(0,HiThr), xlab='Maternal Read Count', ylab='Paternal Read Count', type="n", main=Sample_Name[j]
					, sub = kollapse (  "ABR: ",ABR, ' | MRR: ',MRR, " | SNPs plotted: ", nr_SNPs_plotted, print=F) )
		points (DP, pch =20, cex =.5)
		reads = DP[gene_category == 'escapees', ]
		points (reads[,2]~reads[,1], col='green', pch = '.', cex=5)
		reads = DP[gene_category == 'escapees_het',]
		points (reads[,2]~reads[,1], col='green', pch = '.', cex=5)
		reads = DP[gene_category == 'Impr_PAT', ]
		points (reads[,2]~reads[,1], col="hotpink", pch = '.', cex=5)
		reads = DP[gene_category == 'Impr_MAT', ]
		points (reads[,2]~reads[,1], col="cyan2", pch = '.', cex=5)
		reads = DP[gene_category == 'X', ]
		points (reads[,2]~reads[,1], col=2, pch = '.', cex=5)
		xist = DP[row.names(df_ReadCount) == 'XIST',]
		if (sum(xist)) { points (xist[2]~xist[1] , col='orange', cex=2 )	}
	# 	points (x=0,y=0, col='white', pch=15, cex=1); # makeup: cover 0,0 points

		lx = c('Autosomal', 'X', 'X-esc','Xesc-pred', 'Mat. Expressed', 'Pat. Expressed','XIST')
		legend("topleft", lx,pch = 16,  col =c(1,2,'green','green',"hotpink", "cyan2", 'orange'), inset = .02, cex =.75, bty = "n")
		wplot_save_this(		kollapse("SNP_scatter_",cell, print=F )		)
	}
} # if (linear)


# log2 plots ------------------------------------------------------------------------------------------------------------------------
PlottingThrLo = 0
if (log_plot) {
	OutDir = kollapse (OutDirOrig,"/Scatter_PAT_log2"); if ( !exists (OutDir) ) {dir.create(OutDir)}
	# j =which(cell_IDs =="d2_35_gf")
	for (j in 1:nr_cells) {
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		Maternal_Depth = log2(unlist(df_ReadCount[,mat_colz[j]])); Maternal_Depth[Maternal_Depth == -Inf] = -1
		Paternal_Depth = log2(unlist(df_ReadCount[,pat_colz[j]])); Paternal_Depth[Paternal_Depth == -Inf] = -1
		DP = cbind (Maternal_Depth, Paternal_Depth)
		rownames(DP) = gene_names
		MAE = table(MonAllExp_pat[,cell])
		XBR = percentage_formatter(pc_in_total_of_match(MonAllExp_pat[ gene_category == "X" ,cell],  'Both'))
		# XBialleleicRate[j]  =XBR
		ABR = percentage_formatter(pc_in_total_of_match(MAE, 'Both'))

		nr_SNPs_plotted = 		sum(rowSums(DP)>PlottingThrLo);  # 		sort(table(rowSums(DP)))

		xlb = expression(paste("log"[2],"(Maternal Read Count)")) 		# axis label with subscript
		ylb = expression(paste("log"[2],"(Paternal Read Count)"))
		plot (c(-1, 13), c(-1,13), xlab=xlb, ylab=ylb, type="n", main=Sample_Name[j]
			  , sub = kollapse (  "ABR: ",ABR, ' | XBR: ',XBR, " | SNPs plotted: ", nr_SNPs_plotted, print=F) )
		# 95% lines ------------------------------------------------------------
		if (GuideLines95) {
			y = c(1, 20480); x = y/20
			lines (log2(y),log2(x), col="grey50", lty=2, lwd =2) 		# outside these 2 lines are the 1:20 biased
			lines (log2(x),log2(y), col="grey50", lty=2, lwd =2)
		}
		# points (jitter(DP[rowSums(DP) >-1,]), pch=".", col=1)
		points (jitter(DP[rowSums(DP) >PlottingThrLo,]), pch =20, cex =.5)
		# 		text (reads[,2]~reads[,1], labels =rownames(reads), cex = 0.5)
		reads = 	rbind(jitter(DP[gene_category == 'escapees' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col='green', pch = 16, cex =DotSizeIndividualFigs)
		reads = 	rbind(jitter(DP[gene_category == 'escapees_het' & rowSums(DP) >PlottingThrLo,], factor=0.1))
		points (reads[,2]~reads[,1], col='green', pch = 16, cex =DotSizeIndividualFigs)
		# 		reads = 	rbind(jitter(DP[gene_category == 'Impr_MAT' & rowSums(DP) >-2, ], factor=0.1))
		# 		points (reads[,2]~reads[,1], col="magenta3", pch = 3)
		# 		reads = 	rbind(jitter(DP[gene_category == 'Impr_PAT' & rowSums(DP) >-2, ], factor=0.1))
		# 		points (reads[,2]~reads[,1], col='darkblue', pch = 3)
		reads = rbind(jitter(DP[gene_category == 'X' & rowSums(DP) >PlottingThrLo, ,drop=F], factor=0.1))
		points (reads[,2]~reads[,1], col=2, pch = 16, cex =DotSizeIndividualFigs)
		if ( PrintGeneNames & dim(reads)[1]) {	print(text(reads, labels = rownames(reads), srt = 45, cex = .5, pos =4)) }

		xist = 	(DP[row.names(df_ReadCount) == 'XIST' & rowSums(DP) >PlottingThrLo,])
		if (sum(xist)) { points (xist[2]~xist[1] , col='orange', pch=16, cex =DotSizeIndividualFigs )	}
		if (sum(xist)) { points (xist[2]~xist[1] , col='orange', cex =DotSizeCombinedSupplFig*1.5 )	}

		lx = c('Autosomal', 'Chr X', 'Escapee (X)','XIST')
		legend("topleft", lx,pch = 16,  col =c(1,2,'green','orange'), inset = .02, cex =DotSizeIndividualFigs, bty = "n")

		x = c(-1, 1.1, -1); y = c(-1,-1, 1.1) 										# Draw triangle in unevaluated region
		polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))

		pname = paste("SNP_scatter_log2",cell, cell_type_annot_simple[cell], sep= "_")
		wplot_save_this( plotname = pname)
	}
} # if (log_plot)

# ----------------------------------------------------------------------------------------------------------------------------------------
# Linear Impr plots ------------------------------------------------------------------------------------------------------------------------
PlottingThrLo = 0
HiThr = 400
if (LIN_impr_plot) {
	OutDir = create_set_OutDir(OutDirOrig,OutDir_Impr)
	for (j in 1:nr_cells) {
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		Maternal_Depth = (unlist(df_ReadCount[,mat_colz[j]]))
		Paternal_Depth = (unlist(df_ReadCount[,pat_colz[j]]))
		DP = cbind (Maternal_Depth, Paternal_Depth)
		rownames(DP) = gene_names
		MAE = table(MonAllExp_pat[,cell])
		index = (gene_category_used == "Impr_MAT" |gene_category_used == "Impr_PAT")
		IBR = percentage_formatter(pc_in_total_of_match(MonAllExp_pat[ index, cell],  'Both'))
		# ImprBialleleicRate[j]  = IBR
		ABR = percentage_formatter(pc_in_total_of_match(MAE, 'Both'))

		nr_SNPs_plotted = 		sum(rowSums(DP)>PlottingThrLo);  # 		sort(table(rowSums(DP)))

		xlb = "Maternal Read Count"
		ylb = "Paternal Read Count"
		plot ( c(0, HiThr),  c(0,HiThr), xlab=xlb, ylab=ylb, type="n", main=Sample_Name[j]
			  , sub = kollapse (  "ABR: ",ABR, ' | IBR: ', IBR, " | SNPs plotted: ", nr_SNPs_plotted, print=F) )
		# points (jitter(DP[rowSums(DP) >-1,]), pch=".", col=1)
		points (jitter(DP[rowSums(DP) >PlottingThrLo,]), pch =20, cex =.5)
		# reads = rbind(jitter(DP[gene_category_used == 'X' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		# points (reads[,2]~reads[,1], col='darkred', pch = 17)
		reads = 	rbind(jitter(DP[gene_category_used == 'Impr_MAT' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col="cyan2", pch = 16)
		reads = 	rbind(jitter(DP[gene_category_used == 'Impr_PAT' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col='hotpink', pch = 16)

		# 		lx = c('Autosomal', 'Paternally Expressed','Maternally Expressed', 'Chr X')
		# 		legend("topleft", lx, pch = 16,  col =c(1, 'hotpink', 'cyan2', 'darkred'), inset = .02, cex =0.75, bty = "n")
		lx = c('Autosomal', 'Maternally Expressed','Paternally Expressed')
		legend("topleft", lx, pch = 16,  col =c(1, 'hotpink', 'cyan2'), inset = .02, cex =0.75, bty = "n")

		x = c(-1, 1.1, -1); y = c(-1,-1, 1.1) 										# Draw triangle in unevaluated region
		polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))

		pname = paste("SNP_scatter_impr",cell, cell_type_annot_simple[cell], sep= "_")
		wplot_save_this( plotname = pname)
	}

} # if (LIN_impr_plot)


# log2 Impr plots ------------------------------------------------------------------------------------------------------------------------
# OutDir = create_set_OutDir(OutDirOrig,)
PlottingThrLo = 0
if (log_impr_plot) {
	OutDir = create_set_OutDir(OutDirOrig,OutDir_Impr_log2)
	for (j in 1:nr_cells) {
		# j =which(cell_IDs =="d2_16_gf")
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		Maternal_Depth = log2(unlist(df_ReadCount[,mat_colz[j]])); Maternal_Depth[Maternal_Depth == -Inf] = -1
		Paternal_Depth = log2(unlist(df_ReadCount[,pat_colz[j]])); Paternal_Depth[Paternal_Depth == -Inf] = -1
		DP = cbind (Maternal_Depth, Paternal_Depth)
		rownames(DP) = gene_names
		MAE = table(MonAllExp_pat[,cell])
		index = (gene_category_used == "Impr_MAT" |gene_category_used == "Impr_PAT")
		IBR = percentage_formatter(pc_in_total_of_match(MonAllExp_pat[ index, cell],  'Both'))
		# ImprBialleleicRate[j]  = IBR
		ABR = percentage_formatter(pc_in_total_of_match(MAE, 'Both'))

		nr_SNPs_plotted = 		sum(rowSums(DP)>PlottingThrLo);  # 		sort(table(rowSums(DP)))

		xlb = expression(paste("log"[2],"(Maternal Read Count)")) 		# axis label with subscript
		ylb = expression(paste("log"[2],"(Paternal Read Count)"))
		plot (c(-1, 13), c(-1,13), xlab=xlb, ylab=ylb, type="n", main=Sample_Name[j]
			  , sub = kollapse (  "ABR: ",ABR, ' | IBR: ', IBR, " | SNPs plotted: ", nr_SNPs_plotted, print=F) )
		# 95% lines ------------------------------------------------------------
		if (GuideLines95) {
			y = c(1, 20480); x = y/20
			lines (log2(y),log2(x), col="grey50", lty=2, lwd =2) 		# outside these 2 lines are the 1:20 biased
			lines (log2(x),log2(y), col="grey50", lty=2, lwd =2)
		}
		# points (jitter(DP[rowSums(DP) >-1,]), pch=".", col=1)
		points (jitter(DP[rowSums(DP) >PlottingThrLo,]), pch =20, cex =.5)
		# reads = rbind(jitter(DP[gene_category_used == 'X' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		# points (reads[,2]~reads[,1], col='darkred', pch = 17)
		reads = 	rbind(jitter(DP[gene_category_used == 'Impr_MAT' & rowSums(DP) >PlottingThrLo, , drop=F], factor=0.1))
		points (reads[,2]~reads[,1], col="cyan2", pch = 16, cex =DotSizeIndividualFigs)
		if ( PrintGeneNames & dim(reads)[1]) {	print(text(reads, labels = rownames(reads), srt = 45, cex = .5, pos =4)) }

		reads = 	rbind(jitter(DP[gene_category_used == 'Impr_PAT' & rowSums(DP) >PlottingThrLo, , drop=F], factor=0.1))
		points (reads[,2]~reads[,1], col='hotpink', pch = 16, cex =DotSizeIndividualFigs)
		if ( PrintGeneNames & dim(reads)[1]) {	print(text(reads, labels = rownames(reads), srt = 45, cex = .5, pos =4)) }

		# 		lx = c('Autosomal', 'Paternally Expressed','Maternally Expressed', 'Chr X')
		# 		legend("topleft", lx, pch = 16,  col =c(1, 'hotpink', 'cyan2', 'darkred'), inset = .02, cex =0.75, bty = "n")
		lx = c('Autosomal', 'Maternally Expressed','Paternally Expressed')
		legend("topleft", lx, pch = 16,  col =c(1, 'hotpink', 'cyan2'), inset = .02, cex =DotSizeIndividualFigs, bty = "n")

		x = c(-1, 1.1, -1); y = c(-1,-1, 1.1) 										# Draw triangle in unevaluated region
		polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))

		pname = paste("SNP_scatter_impr_log2",cell, cell_type_annot_simple[cell], sep= "_")
		wplot_save_this( plotname = pname)
	}
	try(dev.off(), silent = T)
} # if (log_impr_plot)

# Supplementary Fig 2A---------------------------------------------------------------------------
OutDir = create_set_OutDir(OutDirOrig,"/Scatter_Suppl_plots")
if (SupplementaryFig2A) {
  try(dev.off(), silent = T)
	# A4: 8.27 x 11.69 inches
	pname = "SupplementaryFig2A_GeneScatter.pdf"
	fname = kollapse(OutDir,"/",pname)
	pdf(fname,width=8.27, height=11.69)
	par(mfrow = c(4, 3))  # 3 rows and 2 columns
	for(j in 1:nr_cells) {
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		Maternal_Depth = log2(unlist(df_ReadCount[,mat_colz[j]])); Maternal_Depth[Maternal_Depth == -Inf] = -1
		Paternal_Depth = log2(unlist(df_ReadCount[,pat_colz[j]])); Paternal_Depth[Paternal_Depth == -Inf] = -1
		# DP_index = log2(unlist(df_ReadCount[,pat_colz[j]] + df_ReadCount[,mat_colz[j]])); Paternal_Depth[Paternal_Depth == -Inf] = -1
		DP = cbind (Maternal_Depth, Paternal_Depth)
		rownames(DP) = gene_names
		MAE = table(MonAllExp_pat[,cell])
		XBR = percentage_formatter(pc_in_total_of_match(MonAllExp_pat[ gene_category == "X" ,cell],  'Both'))
		ABR = percentage_formatter(pc_in_total_of_match(MAE, 'Both'))

		nr_SNPs_plotted = 		sum(rowSums(DP)>PlottingThrLo);  # 		sort(table(rowSums(DP)))

		xlb = expression(paste("log"[2],"(Maternal Read Count)")) 		# axis label with subscript
		ylb = expression(paste("log"[2],"(Paternal Read Count)"))
		plot (c(-1, 13), c(-1,13), xlab=xlb, ylab=ylb, type="n", main=Sample_Name[j]
			  , sub = kollapse (  "ABR: ",ABR, ' | XBR: ',XBR, " | SNPs plotted: ", nr_SNPs_plotted, print=F) )
		# points (jitter(DP[rowSums(DP) >-1,]), pch=".", col=1)
		points (jitter(DP[rowSums(DP) >PlottingThrLo,]), pch =20, cex =.5)
		if (GuideLines95) {
			y = c(1, 20480); x = y/20
			lines (log2(y),log2(x), col="grey50", lty=2, lwd =2) 		# outside these 2 lines are the 1:20 biased
			lines (log2(x),log2(y), col="grey50", lty=2, lwd =2)
		}
		# 		text (reads[,2]~reads[,1], labels =rownames(reads), cex = 0.5)
		reads = 	rbind(jitter(DP[gene_category == 'escapees' & rowSums(DP) > PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col='green', pch = 16, cex =1)
		reads = 	rbind(jitter(DP[gene_category == 'escapees_het' & rowSums(DP) > PlottingThrLo,], factor=0.1))
		points (reads[,2]~reads[,1], col='green', pch = 16, cex =1)
		reads = rbind(jitter(DP[gene_category == 'X' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col=2, pch = 16, cex =1)
		xist = 	(DP[row.names(df_ReadCount) == 'XIST' & rowSums(DP) >PlottingThrLo,])
		if (sum(xist)) { points (xist[2]~xist[1] , col='orange', pch=16, cex =DotSizeCombinedSupplFig )	}
		if (sum(xist)) { points (xist[2]~xist[1] , col='orange', cex =1.5 )	}

		lx = c('Autosomal', 'Chr X', 'Escapee (X)','XIST')
		legend("topleft", lx,pch = 16,  col =c(1,2,'green','orange'), inset = .02, cex =DotSizeCombinedSupplFig, bty = "n")

		x = c(-1, 1.1, -1); y = c(-1,-1, 1.1) 										# Draw triangle in unevaluated region
		polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))
	} # for
	try(dev.off(), silent = T)
} #if

# Supplementary Fig 3A---------------------------------------------------------------------------
if (SupplementaryFig3A) {
  try(dev.off(), silent = T)
	# A4: 8.27 x 11.69 inches
	pname = "SupplementaryFig3A_GeneScatter_Impr.pdf"
	fname = kollapse(OutDir,"/",pname)
	pdf(fname,width=8.27, height=11.69)
	par(mfrow = c(4, 3))  # 3 rows and 2 columns
	for(j in 1:nr_cells) {
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		Maternal_Depth = log2(unlist(df_ReadCount[,mat_colz[j]])); Maternal_Depth[Maternal_Depth == -Inf] = -1
		Paternal_Depth = log2(unlist(df_ReadCount[,pat_colz[j]])); Paternal_Depth[Paternal_Depth == -Inf] = -1
		DP = cbind (Maternal_Depth, Paternal_Depth)
		rownames(DP) = gene_names
		MAE = table(MonAllExp_pat[,cell])
		index = (gene_category_used == "Impr_MAT" | gene_category_used == "Impr_PAT")
		IBR = percentage_formatter(pc_in_total_of_match(MonAllExp_pat[ index, cell],  'Both'))
		ABR = percentage_formatter(pc_in_total_of_match(MAE, 'Both'))

		nr_SNPs_plotted = 		sum(rowSums(DP)>PlottingThrLo);  # 		sort(table(rowSums(DP)))

		xlb = expression(paste("log"[2],"(Maternal Read Count)")) 		# axis label with subscript
		ylb = expression(paste("log"[2],"(Paternal Read Count)"))
		plot (c(-1, 13), c(-1,13), xlab=xlb, ylab=ylb, type="n", main=Sample_Name[j]
			  , sub = kollapse (  "ABR: ",ABR, ' | IBR: ', IBR, " | SNPs plotted: ", nr_SNPs_plotted, print=F) )
		if (GuideLines95) {
			y = c(1, 20480); x = y/20
			lines (log2(y),log2(x), col="grey50", lty=2, lwd =2) 		# outside these 2 lines are the 1:20 biased
			lines (log2(x),log2(y), col="grey50", lty=2, lwd =2)
		}
		# points (jitter(DP[rowSums(DP) >-1,]), pch=".", col=1)
		points (jitter(DP[rowSums(DP) >PlottingThrLo,]), pch =20, cex =.5)
		# reads = rbind(jitter(DP[gene_category_used == 'X' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		# points (reads[,2]~reads[,1], col='darkred', pch = 17)
		reads = 	rbind(jitter(DP[gene_category_used == 'Impr_MAT' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col="cyan2", pch = 16, cex =1)
		reads = 	rbind(jitter(DP[gene_category_used == 'Impr_PAT' & rowSums(DP) >PlottingThrLo, ], factor=0.1))
		points (reads[,2]~reads[,1], col='hotpink', pch = 16, cex =1)

		# lx = c('Autosomal', 'Paternally Expressed','Maternally Expressed', 'Chr X')
		# legend("topleft", lx, pch = 16,  col =c(1, 'hotpink', 'cyan2', 'darkred'), inset = .02, cex =0.75, bty = "n")
		lx = c('Autosomal', 'Paternally Expressed','Maternally Expressed')
		legend("topleft", lx, pch = 16,  col =c(1, 'cyan2', 'hotpink'), inset = .02, cex =DotSizeCombinedSupplFig, bty = "n")

		x = c(-1, 1.1, -1); y = c(-1,-1, 1.1) 										# Draw triangle in unevaluated region
		polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))
		} # for
	try(dev.off(), silent = T)
} #if

# "AllImpr_on1_Plot_Names" ---------------------------------------------------------------------------
if (AllImpr_on1_Plot) {
  try(dev.off(), silent = T)
	jitfactor=1
	pname = "AllImpr_on1_Plot"
	# fname = kollapse(OutDir,"/",pname, ".pdf")
	xlb = expression(paste("log"[2],"(Maternal Read Count)")) 		# axis label with subscript
	ylb = expression(paste("log"[2],"(Paternal Read Count)"))
	plot (c(-1, 13), c(-1,13), xlab=xlb, ylab=ylb, type="n", main="AllImpr_on1_Plot"
		  , sub = "All Impr SNPs on one plot with 10% and 2% allelic bias boundaries. \n D1 black, D2 red, D3 green, D4 blue, D5 baby-blue" )
	for(j in 1:nr_cells) {
		cell = substr(colnames(df_ReadCount)[mat_colz[j]],1,8); 	cell
		if (GermOnly){ if (cell %in% which_names(Control)) {next} }
		Maternal_Depth = log2(unlist(df_ReadCount[,mat_colz[j]])); Maternal_Depth[Maternal_Depth == -Inf] = -1
		Paternal_Depth = log2(unlist(df_ReadCount[,pat_colz[j]])); Paternal_Depth[Paternal_Depth == -Inf] = -1
		DP = cbind (Maternal_Depth, Paternal_Depth)
		rownames(DP) = gene_names

		reads = 	rbind((DP[gene_category_used == 'Impr_MAT' & rowSums(DP) >PlottingThrLo, , drop=F]))
		if(DonorColor){		points (reads[,2]~reads[,1], col=as.numeric(as.factor(donor)[cell]), pch = 21, cex =.66)}
		if(CellTypeColor){	points (reads[,2]~reads[,1], col=cell_type_col_simple[cell], pch = 21, cex =.66)}
		points (reads[,2]~reads[,1], col="cyan2", pch = c(10,13), cex =.3)
		if ( PrintGeneNames & dim(reads)[1] > 0) {	print(text(reads, labels = rownames(reads), srt = 45, cex = labelSize, pos =4)) }

		reads = 	rbind((DP[gene_category_used == 'Impr_PAT' & rowSums(DP) >PlottingThrLo, , drop=F]))
		if(CellTypeColor){	points (reads[,2]~reads[,1], col=cell_type_col_simple[cell], pch = 21, cex =.66)}
		if(DonorColor){		points (reads[,2]~reads[,1], col=as.numeric(as.factor(donor)[cell]), pch = 21, cex =.66)}
		points (reads[,2]~reads[,1], col='hotpink', pch = c(10,13), cex =.3)
		if ( PrintGeneNames & dim(reads)[1] > 0) {	print(text(reads, labels = rownames(reads), srt = 45, cex = labelSize, pos =4)) }

	} # for

	lx = c('Autosomal', 'Paternally Expressed','Maternally Expressed')
	legend("topleft", lx, pch = 16,  col =c(1, 'cyan2', 'hotpink'), inset = .02, cex =0.75, bty = "n")
	x = c(-1, 1.1, -1); y = c(-1,-1, 1.1) 										# Draw triangle in unevaluated region
	polygon(x, y, density = NULL, angle = 45, border = NA, col = rgb(0, 0, 0, 0.5))
	# lines ------
	y = c(1,4,16,64,256,2048, 20480)
	x = y/10
	lines (log2(y),log2(x), col=3, lty=3) 		# outside these 2 lines are the 1:10 biased
	lines (log2(x),log2(y), col=3, lty=3) 		# at -1 it turns parallel with the axis
	x = y/50
	lines (log2(y),log2(x), col=4, lty=3) 		# outside these 2 lines are the 1:50 biased, Monoallelically called genes
	lines (log2(x),log2(y), col=4, lty=3)

	# try(dev.off(), silent = T)
fname = "AllImpr_on1_Plot_Names"

if(DonorColor){		fname = kollapse(fname, "_DonorColor")}
if(CellTypeColor){	fname = kollapse(fname, "_CellTypeColor")}
if(GermOnly){	fname = kollapse(fname, "_GermOnly")}
wplot_save_this(fname)

} #if

# "AllX_on1_Plot" ---------------------------------------------------------------------------
# Basic discovery of which genes tend to lag behind in reactivation
if (AllX_on1_Plot) {
	source ("~/analysis/DevCell_analysis/13.zz.PlotAllChrX_on1plot_NAMES.R")
} #if

# Finish
OutDir = OutDirOrig


# ShowLineOfMonoAllelicity =F
# if (ShowLineOfMonoAllelicity) {
# 	y = c(1,4,16,64,256,2048, 20480)
# 	x = y/10
# 	lines (log2(y),log2(x), col=3, lty=3) 		# outside these 2 lines are the 1:10 biased
# 	lines (log2(x),log2(y), col=3, lty=3) 		# at -1 it turns parallel with the axis
# 	x = y/50
# 	lines (log2(y),log2(x), col=4, lty=3) 		# outside these 2 lines are the 1:50 biased, Monoallelically called genes
# 	lines (log2(x),log2(y), col=4, lty=3)
# }
