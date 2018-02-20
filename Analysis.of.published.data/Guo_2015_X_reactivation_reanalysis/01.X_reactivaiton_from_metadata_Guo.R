######################################################################################################
# 01.X_reactivaiton_from_metadata_Guo.R
######################################################################################################

"This script writes in two directories"

# Functions  ------------------------------------------------------------------------------------
require(MarkdownReports) # https://vertesy.github.io/MarkdownReports/
try (source ('~/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F)


# Setup  ------------------------------------------------------------------------------------

InputDir = "~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Guo_X_reactivation/data/"
setup_MarkdownReports(OutDir = "~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Guo_X_reactivation/01.a.Transcript_Barplots",
                      scriptname = "01.X_reactivaiton_from_metadata_Guo.R",
                      title = "X reactivaiton from metadata of Guo et al")
OutDirOrig =OutDir
ls_F = list.files(InputDir, pattern = "^X_ratios_w"); ls_F
nrf = l(ls_F); nrf

# Plot Each transcript  ------------------------------------------------------------------------------------

llprint("### Per gene barplots (as in paper) recreated for each transcript.")
Ls_emb= NULL
for (f in 1:nrf) {
  print(f)
  Ls_emb[[f]] = read.simple.tsv(InputDir,ls_F[f])
  X_ratio = Ls_emb[[f]][, -(2:3)]

  cellzz=colnames(X_ratio)[-1]
  cellzz[grep("Soma",cellzz)] =  "SOMA"
  cellzz[grep("PGC",cellzz)] =  "PGC"
  # colnames(X_ratio) = c("Gene", cellzz)
  r=1
  for (r in 1:nrow(X_ratio) ) {
    Expr = 100*as.numeric(X_ratio[r, -1])
    names(Expr) = cellzz
    Expr = na.omit.strip(Expr)
    pname= kollapse( "Transcript.", X_ratio$Gene[r], " in Embryo ",substr(ls_F[f], 10,11))
    wbarplot(rbind(Expr,100-Expr), col = 3:2, plotname = pname, main =pname, mdlink = T)
  }
}


# No_Difference_Between_Escapees_and_X_genes ------------------------------------------------------------------------------------
setup_MarkdownReports(OutDir = "~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Guo_X_reactivation/01.b.Transcript_Barplots_QC",
                      scriptname = "01.X_reactivaiton_from_metadata_Guo.R",
                      title = "X reactivaiton from metadata of Guo et al")

llprint("## No_Difference_Between_Escapees_and_X_genes ")


# inline_vec.num(na.omit.strip(fromClipboard()))

Esc = list(
  "w4" = c(c( 1, 0.020833, NaN, NaN, 0, 0 ), c( 0, 0, 0.271845, 0.344828, 0, 0, 1, 1, 1, 1, 1, 0.027778 )),
  "w8" = c(c( 1, 0, 1 ), c( 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1 )),
  "w10" = c( 0, 0, 0, 0 ),
  "w11" = c( 0, 1, 0.005988, 1, 0.006944, 0, 1, 0, 0.003636 ),
  "w17" = c( 0, 1, 0, 0.998066, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 )
)

Proper_X = list(
  "w4" = NaN,
  "w8" = c( 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0.994792, 1, 1, 1, 1, 1, 1, 1, 0.011173, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 ),
  "w10" = NaN ,
  "w11" = c( 1, 0, 1 ),
  "w17" = c( 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0.5, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 )
)


No_Difference_Between_Escapees_and_X_genes = list( "Escapees" = 100*unlist(Esc),
                                                   "Proper X genes" = 100*unlist(Proper_X))

llprint("## Escapee genes do not escape X inactivation in somatic cells in Guo et al")
llprint("There is one cell with 2 SNPs in one gene that show biallelic expression")
wstripchart(No_Difference_Between_Escapees_and_X_genes, incrBottMarginBy = 2, col= 3:4, colorbyColumn = T, jitter = .4,
            ylab="% reads from the alternative allele", mdlink = T)


# Xinactivaiton_is_not_Random ------------------------------------------------------------------------------------------------
llprint("## There are 6 (20%) cells that express at least one allele inconsistent with haplotype expression of other cells")
llprint("I counted all cells with at least 2 genes. At 2 genes we probably have little chance discover inconsistency")

Xinactivaiton_is_not_Random = matrix(data = NaN, nrow = 5, ncol = 4, dimnames = list(c("w4", "w8", "w10", "w11", "w17"), c("Haplotype1","Haplotype2","InconsistentHapl", "TotalSomaCellsWithMoreThan1MonoAllGene")))
# Inconsistemnt if there is one gene that is not expressed from the same haplotype

# SomaCellsWithMoreThan1MonoAllGene

Xinactivaiton_is_not_Random[1,] = c( 2,2,1,5 )
Xinactivaiton_is_not_Random[2,] = c( 6,0,4,10 )
Xinactivaiton_is_not_Random[3,] = c( 0,2,0,2 )
Xinactivaiton_is_not_Random[4,] = c( 2,1,0,3 )
Xinactivaiton_is_not_Random[5,] = c( 6,2,1,9)

Xinactivaiton_is_not_Random_ = t(Xinactivaiton_is_not_Random[,1:3 ] / Xinactivaiton_is_not_Random[,4 ])
colnames(Xinactivaiton_is_not_Random_) = paste(colnames(Xinactivaiton_is_not_Random_), "(",Xinactivaiton_is_not_Random[,4 ], ")")
wbarplot(Xinactivaiton_is_not_Random_, col = 2:4, ylimits = c(0,1.5) , sub = "Number of cells in bracket", mdlink = T)
legend("topleft", legend = c("Haplotype1","Haplotype2","InconsistentHapl"), fill = 2:4, bty="n")
wplot_save_this(plotname = plotnameLastPlot)


llprint("## There is some synchrony of allelic usage among the cells, but its only significant at w8, without correction for multiple testing")
llprint("I would consider it non-siginificant in total")

Xinactivaiton_is_not_Random_2 = t(Xinactivaiton_is_not_Random[,1:2 ] / rowSums(Xinactivaiton_is_not_Random[,1:2 ]))
colnames(Xinactivaiton_is_not_Random_2) = paste(colnames(Xinactivaiton_is_not_Random_2), "(",rowSums(Xinactivaiton_is_not_Random[,1:2 ]), ")")
wbarplot(Xinactivaiton_is_not_Random_2, col = 2:4, ylimits = c(0,1.5), hline = 1, filtercol = F, sub = "Number of cells in bracket", mdlink = T)
legend("topleft", legend = c("Haplotype1","Haplotype2"), fill = 2:3, bty="n")
wplot_save_this(plotname = plotnameLastPlot)

# Pie Charts per embryo ------------------------------------------------------------------------------------

llprint("## Escapee and proper X genes have different biases")

Expressed_Genotypes_17_in_X_Genes = table(round(unlist(Proper_X)))
Expressed_Genotypes_11_in_Esc_Genes = table(round(unlist(Esc)))
Expressed_Genotypes_All_28_Genes_on_X = Expressed_Genotypes_17_in_X_Genes+ Expressed_Genotypes_11_in_Esc_Genes
names(Expressed_Genotypes_17_in_X_Genes) = names(Expressed_Genotypes_11_in_Esc_Genes) = names(Expressed_Genotypes_All_28_Genes_on_X) =  c("Ref", "Alt")


wpie(Expressed_Genotypes_17_in_X_Genes, mdlink = T)
wpie(Expressed_Genotypes_11_in_Esc_Genes, mdlink = T)
wpie(Expressed_Genotypes_All_28_Genes_on_X)

#   ------------------------------------------------------------------------------------






