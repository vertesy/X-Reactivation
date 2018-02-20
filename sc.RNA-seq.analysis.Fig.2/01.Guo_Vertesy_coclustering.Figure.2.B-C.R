######################################################################################################
# 01.Guo_Vertesy_coclustering.Figure.2.B-C.R
######################################################################################################
# source ("~/analysis/Scripts_Transcriptome_Analysis/NewTranscriptomeAnalysis/01.CoClustering/01.Guo_Vertesy_coclustering.Figure.2.B-C.R")
try(dev.off(),silent = T)


# Functions ------------------------
wlegend.old <- function(x=c("topleft", "topright", "bottomright", "bottomleft")[4],
                        legend, fill = NULL, ..., bty = "n", OverwritePrevPDF =T) { # Add a legend, and save the plot immediately
  legend(x=x,legend=legend,fill=fill, ..., bty=bty)
  if (OverwritePrevPDF) {   wplot_save_this(plotname = plotnameLastPlot)  }
}

# Setup ----------------------------
setup_MarkdownReports(OutDir = OutDir, scriptname = "01.Guo_Vertesy_coclustering.Figure.2.B-C.R")


# Parameters ----------------------------
distance_metric4MDS = "euclidean"

# Metadata  ------------------------------------------------------------------------------------------------------------



# MDS ------------------------------------------------------------------------------------------------------------
llprint("## Multidimensional Scaling")

TPM_subset = Expression_All_Filtered

FPKM_log2 = cor(log10(TPM_subset+1))

FPKM_log.dist <- dist(t(FPKM_log2), method = distance_metric4MDS) # Eucledian Distance Matrix Calculation
FPKM_log2.fit <- cmdscale(FPKM_log.dist, eig=TRUE, k = 2) # Classical (Metric) Multidimensional Scaling

VariationExplainedByPCs = 100 *FPKM_log2.fit$eig/sum(FPKM_log2.fit$eig)
wbarplot(VariationExplainedByPCs, ylab="Variation Explained (%)")
ve = iround(VariationExplainedByPCs[1:2])
txt = c("PC1", "PC2")
namez =paste0(txt, " (",ve,  "%)")

x.goi.log2 <- FPKM_log2.fit$points[,1]
y.goi.log2 <- FPKM_log2.fit$points[,2]

FemaleEmbryonicSingleCells = data.frame(
	"PC1 (x %)" = FPKM_log2.fit$points[,1],
	"PC2 (y %)" = -FPKM_log2.fit$points[,2]
)
colnames(FemaleEmbryonicSingleCells)= namez

FindBack = Sample_Name[rownames(FemaleEmbryonicSingleCells)]
OrigNames = rownames(FemaleEmbryonicSingleCells)


xlim = range(FemaleEmbryonicSingleCells[,1])
ylim = range(FemaleEmbryonicSingleCells[,2])

LQ_cells = sort(LQ_cells)
ccc= clusterCol

rownames(FemaleEmbryonicSingleCells)
wplot(FemaleEmbryonicSingleCells, plotname = "Fig1.C.PCA_CellState", pch =16, h=5, col = ccc, axes=F,frame.plot=T, mdlink = T)
wlegend.old("bottomleft", legend =names(ColorList$CellType),  fill = ColorList$CellType, bty = "n")


wplot(FemaleEmbryonicSingleCells, plotname = "PCA_Names", type ="n")
text(FemaleEmbryonicSingleCells, rownames(FemaleEmbryonicSingleCells), cex = .3, srt =45)
wplot_save_this(plotname = plotnameLastPlot)

# oo(xx= rownames(FemaleEmbryonicSingleCells)
ccc= DonorCol
wplot(FemaleEmbryonicSingleCells, plotname = "Fig1.B.PCA_EmbryonicAge", pch =16, h=5, col = ccc, axes=F,frame.plot=T)
wlegend.old("bottomleft", legend = names(ColorList$Donor),  fill = ColorList$Donor, bty = "n")


#  ------------------------------------------------------------------------------------------------------------

