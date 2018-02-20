Fig.3c =T
if (Fig.3c) {
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
