######################################################################################################
# Genes in non-reactivated
######################################################################################################
# source("./analysis/02_AlleleSpecificExpression/13.z.Reviewer1_q26.R")

# Setup ----------------------------
require("DESeq2"); # source("https://bioconductor.org/biocLite.R"); biocLite("DESeq2")

#plot differentially expressed genes

# if(!require(package)) install.packages(package) } # install package if cannot be loaded
setup_MarkdownReports(OutDir = "~/Google_Drive/X_react_Data/Abelz/X_react_Paper2/13.h.X-Reactivation/DiffExp_RaceID",
                      scriptname = "13.z.Reviewer1_q26.raceID.R")
# create_set_SubDir("DeSeq")

if (!exists("TPM")) { TPM = read.simple.tsv("~/Google_Drive/X_react_Data/cout_tables/AllCells/rsem.tpm_table_all.tsv") }

# Functions -----------------------------------------------------------------------

plotdiffgenesnb <- function(diffexpnb_object, pthr=.05, lthr=1, mthr=0, xname="A", yname="B", bgcol="grey", ppch = 20,
                            lcol = "grey33", show_names=TRUE, padj=TRUE, draw_mthr =T, draw_lthr =T, dontplotbelow=1, ...){
  y <- diffexpnb_object$res

  xlname =  paste0("log2 ((mRNA[",xname,"] + mRNA[",yname,"])/2)")
  ylname =  paste0("log2 (mRNA[",yname,"] / mRNA[",xname,"])")
  meanExpr = log2( (y$baseMeanA + y$baseMeanB)/2 )
  if (!is.na(dontplotbelow)) { meanExpr[meanExpr<dontplotbelow] = NA} # exlcude lowly expressed

  plot(meanExpr, y$log2FoldChange, pch=ppch, xlab=xlname , ylab=ylname, col=bgcol, ...)
  if (draw_mthr) { abline(v=mthr, lty=3, col=lcol) }
  if (draw_lthr) { abline(h=c(lthr, -lthr), lty=3, col=lcol)  }

  if ( ! is.null(pthr) ){
    if ( padj ) f <- y$padj < pthr else f <- y$pval < pthr
    points(meanExpr[f], y$log2FoldChange[f], col="red",pch=20)
  }
  if ( !is.null(lthr) ) f <- f & abs( y$log2FoldChange ) > lthr
  if ( !is.null(mthr) ) f <- f &  (log2( (y$baseMeanA + y$baseMeanB)/2 ) ) > mthr
  if ( show_names )  text(meanExpr[f],y$log2FoldChange[f],id2name(rownames(y))[f],cex=.5)

}

# Parameters ----------------------------
EGC_only = T

GermNamez_HQ = c( 'd1_01_gf', 'd1_02_gf', 'd1_03_gf', 'd1_04_gf', 'd1_05_gf', 'd1_06_gf', 'd1_27_gf', 'd1_28_gf', 'd1_29_gf', 'd2_12_gf', 'd2_13_gf', 'd2_14_gf', 'd2_15_gf', 'd2_16_gf', 'd2_32_gf', 'd2_33_gf', 'd2_52_gf', 'd2_53_gf', 'd2_54_gf', 'd2_56_gf', 'd2_58_gf', 'd4_01_gf', 'd4_03_gf', 'd4_04_gf', 'd4_05_gf', 'd4_06_gf', 'd4_07_gf', 'd4_08_gf', 'd4_09_gf', 'd4_10_gf', 'd4_31_gf', 'd4_32_gf', 'd4_33_gf', 'd4_34_gf', 'd4_35_gf', 'd4_36_gf', 'd4_37_gf', 'd4_38_gf', 'd5_14_gf', 'd5_15_gf', 'd5_16_gf', 'd5_17_gf', 'd5_18_gf', 'd5_19_gf', 'd5_20_gf', 'd5_21_gf', 'd5_22_gf', 'd5_23_gf', 'd5_66_gf', 'd5_67_gf', 'd5_69_gf', 'd5_70_gf', 'd5_72_gf')
GermNamez_HQ_PGCs = c( 'd1_01_gf', 'd1_02_gf', 'd1_03_gf', 'd1_04_gf', 'd1_05_gf', 'd1_06_gf', 'd1_27_gf', 'd1_28_gf', 'd1_29_gf', 'd2_12_gf', 'd2_13_gf', 'd2_14_gf', 'd2_15_gf', 'd2_16_gf', 'd2_32_gf', 'd2_33_gf', 'd2_52_gf', 'd2_53_gf', 'd2_54_gf', 'd2_56_gf', 'd2_58_gf', 'd4_01_gf', 'd4_05_gf', 'd4_06_gf', 'd4_07_gf', 'd4_08_gf', 'd4_09_gf', 'd4_31_gf', 'd4_32_gf', 'd4_33_gf', 'd4_36_gf', 'd4_37_gf', 'd4_38_gf', 'd5_15_gf', 'd5_16_gf', 'd5_18_gf', 'd5_19_gf', 'd5_23_gf', 'd5_70_gf')
GermNamez_HQ_AGCs = c( 'd4_03_gf', 'd4_04_gf', 'd4_10_gf', 'd4_34_gf', 'd4_35_gf', 'd5_14_gf', 'd5_17_gf', 'd5_20_gf', 'd5_21_gf', 'd5_22_gf', 'd5_66_gf', 'd5_67_gf', 'd5_69_gf', 'd5_72_gf')
chrX_bias = c( 'd1_02_gf', 'd1_03_gf', 'd1_04_gf', 'd1_06_gf', 'd2_12_gf', 'd2_15_gf', 'd2_16_gf', 'd2_58_gf', 'd4_05_gf', 'd4_07_gf', 'd4_31_gf')

Cells = GermNamez_HQ
Cells = if(EGC_only) Cells[cell_type_annot[Cells] == "PGC"]
no_chrX_bias = setdiff(Cells, chrX_bias)

sc2 =  sc
dim(sc@fdata)
setdiff(no_chrX_bias,colnames(sc@fdata))
deo <-diffexpnb(sc@fdata[, GermNamez_HQ], cells_of_interest = chrX_bias, cells_background = no_chrX_bias, norm=FALSE, vfit=sc@background$vfit, method="pooled") # , logreg=FALSE,

varname = paste("Diff. Expressed genes between XiXa and XaXa")
p$pval = .05
plotdiffgenesnb(deo, xname="XaXa",yname="XiXa",pthr = p$pval , lthr=0, mthr=p$mthr, dontplotbelow = 0, show_names=T, padj=T, main=varname) # p$lthr LOG fold change theshold;  p$mthr Min LOG mean expression thr
wplot_save_this(plotname = varname, mdlink = T, w=14,h=14)

DE <- iround(deo$res[which(deo$res$"padj" < p$pval ),]); dim(DE) #select significant genes
View(DE)

thr_fc_GeoM = 2
thr_fc_Mean = 2
deo$res$foldChange_GeoMean >thr_fc_GeoM

write.simple.tsv(DE)

idx= grep("SNAR", rownames(DE), invert = T, value = T)


pdfA4plot_on(pname = "DiffExpGenes")
for (i in 1:l(idx)) {
  g = idx[i]
  ixx = as.numeric(GermNamez_HQ_PGCs %in% chrX_bias);l(ixx)
  ixx = translate(ixx, oldvalues = 1:0, newvalues = c("XiXa", "XaXa"))
  ls_Ex = split(as.numeric(sc@fdata[g, GermNamez_HQ_PGCs]),f = ixx )
  wstripchart(ls_Ex, plotname = g, jitter = .4, colorbyColumn = T, col = c(2:3), tilted_text = T, savefile = F)
}
pdfA4plot_off()

