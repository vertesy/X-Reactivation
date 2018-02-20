##########################################################################################
### 01.TissueSpecificImprintedGenes.R
##########################################################################################
# source("~/analysis/Baran_2015_TissueSpecificImprinting/01.TissueSpecificImprintedGenes.R")

# Functions  ------------------------------------------------------------------------------------


# Setup  ------------------------------------------------------------------------------------

setup_MarkdownReports(OutDir = "~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Imprinting/Baran_2015_TissueSpecificImprinting",
                      scriptname = "01.TissueSpecificImprintedGenes.R",
                      title = "Tissue Specific Imprinting is more common among Orphan imprinted genes, while IC member genes are typically body-wide imprinted" )


llprint("#### Reanalyzing  results from:")
llprint( "> Baran, Y., Subramaniam, M., Biton, A., Tukiainen, T., Tsang, E. K., Rivas, M. A., … Lappalainen, T. (2015). **The landscape of genomic imprinting across diverse adult human tissues.** Genome Research, 25(7). http://doi.org/10.1101/gr.192278.115")

MetadataDir = "./metadata/Imprinting/"


## Metadata Run once ------------------------
# IC_Membership_All = read.simple.tsv.named.vector("./metadata/Imprinting/IC_Membership_All.tsv")
# IC_Membership_All = names(IC_Membership_All)
# Imprinted_confirmed_union_HUGO_MatExpr = read.simple.vec("./metadata/Imprinting/Imprinted_confirmed_union_HUGO_MatExpr.vec")
# Imprinted_confirmed_union_HUGO_PatExpr = read.simple.vec("./metadata/Imprinting/Imprinted_confirmed_union_HUGO_PatExpr.vec")
#
# OrphanGenes = symdiff(union(Imprinted_confirmed_union_HUGO_MatExpr, Imprinted_confirmed_union_HUGO_PatExpr), IC_Membership_All)
# OrphanGenes = sort(OrphanGenes[[1]])

# inline_vec.char(OrphanGenes)
# Manually excluded: , 'KCNQ1DN'
# inline_vec.char(IC_Membership_All)

# Metadata ------------------------
OrphanGenes = c( 'AIM1', 'ANO1', 'ATP10A', 'CPA4', 'DIRAS3', 'DLGAP2', 'DLX5', 'DNMT1', 'FAM50B', 'GDAP1L1', 'GLIS3', 'GNAS', 'GNAS-AS1', 'GPR1', 'HYMAI', 'INPP5F', 'KCNK9', 'KLF14', 'LIN28B', 'LRRTM1', 'MAGI2', 'MCTS2P', 'MEST', 'MESTIT1', 'MIMT1', 'MIR296', 'MIR298', 'MIR371A', 'NAA60', 'NAP1L5', 'NLRP2', 'NNAT', 'NTM', 'PEG10', 'PLAGL1', 'PPP1R9A', 'RB1', 'RBP5', 'SGCE', 'SGK2', 'SLC22A2', 'SLC22A3', 'SNORD107', 'SNORD108', 'SNORD109A', 'SNORD109B', 'SNORD115-45', 'SNORD64', 'TCEB3C', 'TFPI2', 'TP73', 'UBE3A', 'WT1', 'WT1-AS', 'ZC3H12C', 'ZDBF2', 'ZFAT', 'ZNF597')
IC_Membership_All = c( 'H19', 'IGF2', 'IGF2-AS', 'INS', 'KCNQ1', 'KCNQ1OT1', 'CDKN1C', 'PHLDA2', 'DLK1', 'MEG3', 'RTL1', 'MEG8', 'SNORD111', 'SNORD112', 'SNORD113', 'SNORD114', 'MIRG', 'DIO3', 'MKRN3', 'MAGEL2', 'NDN', 'PWRN1', 'C15orf2', 'SNRPN', 'SNURF', 'IPW', 'USP29', 'ZNF264', 'PEG3', 'PEG3-AS1', 'ZIM1', 'ZIM2', 'ZIM3')


Baran_Table2 = read.simple.tsv(MetadataDir,"Baran_2015.TissueSpecificImprinting.tsv")
attach_w_rownames(Baran_Table2)
GeneNamesGTEx = names(Percent_Imprinted)

OrphanGenes_found = intersect(GeneNamesGTEx, OrphanGenes)
IC_Members_found = intersect(GeneNamesGTEx, IC_Membership_All)

PercentageOfTissuesImprinted_OrphanGenes 	= 100 * Percent_Imprinted[OrphanGenes_found]
PercentageOfTissuesImprinted_IC_Members 	= 100 *Percent_Imprinted[IC_Members_found]



# Stripchart comparison of the categories -----------------------------------
llprint("## IC member genes tend to be more body-wide imprinted" )

PercentageOfTissuesImprinted = list(
  "OrphanGenes" = PercentageOfTissuesImprinted_OrphanGenes,
  "IC_Members" = PercentageOfTissuesImprinted_IC_Members
)

# Define the color as the number of total tissues where the observation is made
TotalUsefulObservations = biallelic+imprinted
EvidenceLevel = val2col(TotalUsefulObservations); names(EvidenceLevel) = names(TotalUsefulObservations)

EvidenceLevel_ls = list(  EvidenceLevel[OrphanGenes_found],
                          EvidenceLevel[IC_Members_found])

MWW = wilcox.test( PercentageOfTissuesImprinted[[1]], PercentageOfTissuesImprinted[[2]] , alternative = "less")
llprint ( 'TSI in IC VS non-IC genes. p-value @ MWW: ', iround(as.numeric(MWW[3]),2) )


rng = range(TotalUsefulObservations[union(OrphanGenes_found, IC_Members_found)])
subb =  paste0('p-value @ MWW: ', iround(as.numeric(MWW[3]),2), " | Red: more evidence/nr of tissues. Range:",rng[1]," to ",rng[2] )
wstripchart_list(PercentageOfTissuesImprinted, jitter = .35, bg = EvidenceLevel_ls, tilted_text = T, sub=subb, incrBottMarginBy = 1, mdlink = T)

TSIs = iround(100- unlapply(PercentageOfTissuesImprinted, median))
MarkDown_Table_writer_NamedVector(TSIs, title_of_table = "% of tissues leaking")


# Evidence level does not explain global/tissue spec. imprinting -----------------------------------
llprint("## Evidence level does not explain global/tissue spec. imprinting" )

PercentageImprinted = c(PercentageOfTissuesImprinted_OrphanGenes, PercentageOfTissuesImprinted_IC_Members)
DispOrder = names(PercentageImprinted)

TissueSpecificity = cbind ("Number of Organs Evaluated" = TotalUsefulObservations[DispOrder],
                      "Percentageof Organs Imprinted " = PercentageImprinted)

colz = c(rep(2,l(OrphanGenes_found)),rep(4,l(IC_Members_found)))
wplot(TissueSpecificity, col = colz, mdlink = T)
legend("bottomleft", fill = c(2,4), legend = c("OrphanGenes","IC_Members"), bty = "n")
wplot_save_this(plotname = "TissueSpecificity_in_IC_and_Orphan_genes", mdlink = T)

text(jitter(TissueSpecificity, amount = 1), labels = rownames(TissueSpecificity), srt=45, cex=.5)
wplot_save_this(plotname = "TissueSpecificity_in_IC_and_Orphan_genes.wNames", w=14, h=14)


# Evidence level does not explain global/tissue spec. imprinting -----------------------------------
# Reviewer1q19 =T
# if (Reviewer1q19) {
#   SomaPerWeek = splititsnames_byValues(weeks[ControlNamez_HQ])
#   FigS.Y.NoIncreaseInSomaticImprinting = SomaPerWeek
#   for (j in 1:l(SomaPerWeek)) {
#     FigS.Y.NoIncreaseInSomaticImprinting[[j]] = 100*na.omit.strip(iPE_matrix[ OrphanGenes, SomaPerWeek[[j]]])
#   }
#   wstripchart(FigS.Y.NoIncreaseInSomaticImprinting, col = 1:5, colorbyColumn = T, border = "grey50", jitter = .4,
#               ylb = "% reads from the expected allele", tilted_text = T)
# } #


Reviewer1q20 =T
if (Reviewer1q20) {
  iPE_matrix

  lss = list(
	"Germ Cells - Orphan Genes" = 		getRows(iPE_matrix[ , GermNamez_HQ], OrphanGenes, removeNAonly = T) ,
	"Germ Cells - IC Genes" = 			getRows(iPE_matrix[ , GermNamez_HQ], IC_members, removeNAonly = T),
	"Somatic Cells - Orphan Genes" =	getRows(iPE_matrix[ , ControlNamez_HQ], OrphanGenes, removeNAonly = T),
	"Somatic Cells - IC Genes" = 		getRows(iPE_matrix[ , ControlNamez_HQ], IC_members, removeNAonly = T)
	)

  try.dev.off()
  pdfA4plot_on("FigS.X.Bias_Expression", rows = 3, cols = 2)
  for (i in 1:4) {

    x = lss[[i]]

    igg = rownames(x)
    pgg=colsplit(t(x), igg); names(pgg)=igg
    pgg =lapply(pgg, na.omit.strip)
    GermGenes = cbind(
      "Number of cells with measurable alllelic expression" = jitter(unlapply(pgg, l),amount = .5),
      "% reads from paternal allele"= jitter(100*unlapply(pgg, mean),amount = 2.5)
    )

    val=iround(log10(rowMeans(DP_SUM_parental_HQ[igg, colnames(x)], na.rm = T)+1))
    ccc = val2col( val) # , col = terrain.colors(10)
    wplot(GermGenes, plotname = names(lss)[i], savefile = F, bg = ccc, pch = 23)
    abline(v = 5, lty=2, col="grey")
    # rect(-1,-10,5,110, col = "grey")
    text(GermGenes, labels = igg, srt =45, cex=.5, pos = 2)
    llb = round(range(rowMeans(DP_SUM_parental_HQ[igg, colnames(x)])), digits = 1)
    legend("bottomright",legend = llb,  fill = val2col(llb), bty = "n", title = "Mean Expression")
  }

  par(mfrow = c(3, 1)) # full width plot
  x = 100*iPE_matrix[ , GermNamez_HQ]
  igg = rownames(x)
  pgg=colsplit(t(x), igg); names(pgg)=igg
  pgg =lapply(pgg, na.omit.strip)
  pgg =lapply(pgg, jitter, amount = 2.5)

  ccc = (igg %in% IC_members)+2
  ord = names(sort(unlapply(pgg, mean)) )
  wstripchart(pgg[ord], savefile = F, pchcex = .5, colorbyColumn = F, col = ccc, tilted_text = T,pch = c(24,25), bg = NULL
              , ylb = "% reads from paternal allele", main = "Raw allelic bias for each imprinted gene", jitter = .4 )

  lbl = unlapply(pgg, l)[ord]
  barplot_label(barplotted_variable = unlapply(pgg[ord], mean), labels = lbl, bottom = T)

  pdfA4plot_off()
  }



