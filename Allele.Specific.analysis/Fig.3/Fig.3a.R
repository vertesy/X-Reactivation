Fig.3a.reviewer3_q1 =T
if (Fig.3a.reviewer3_q1) {

  # x_mat = iPE_matrix
  ls_gene_category =  splititsnames_byValues(gene_category)
  iGenes_ = c(ls_gene_category$Impr_PatExpr, ls_gene_category$Impr_MatExpr)

  x_mat = DP_SUM_parental_HQ[ iGenes_, GermNamez_HQ]
  x_mat[x_mat==0] = NA


  x = colsplit(t(x_mat+1), f=rownames(x_mat))
  x = lapply(x, na.omit.strip)
  x = lapply(x, log2)
  l_X = sort(unlapply(x, length), decreasing = T)
  l_X = l_X[l_X>2]
  FigS.All_Imprinted_Genes.PGC = x[names(l_X)]
  wstripchart(FigS.All_Imprinted_Genes.PGC, pchcex = 1, jitter = .4, border = "grey", ylb = "log2(Allleic Read Count +1)", tilted_text = T, w=14)
  barplot_label(barplotted_variable = unlapply(FigS.All_Imprinted_Genes.PGC, mean), labels = unlapply(FigS.All_Imprinted_Genes.PGC, length), bottom = T, w=14)

  x_mat_m = DP_mat[ iGenes_, GermNamez_HQ]
  x_mat_p = DP_pat[ iGenes_, GermNamez_HQ]


  x = colsplit(t(x_mat+1), f=rownames(x_mat))
  x = lapply(x, na.omit.strip)
  x = lapply(x, log2)

  xm_mat_m = DP_mat[ iGenes_, GermNamez_HQ]
  xm_mat_m[is.na(x_mat)] = NA

  xm = colsplit(t(xm_mat_m+1), f=rownames(xm_mat_m))
  xm = lapply(xm, na.omit.strip)
  xm = lapply(xm, log2)
  # l_Xm = sort(unlapply(xm, length), decreasing = T)

  xp_mat_p = DP_pat[ iGenes_, GermNamez_HQ]
  xp_mat_p[is.na(x_mat)] = NA
  xp = colsplit(t(xp_mat_p+1), f=rownames(xp_mat_p))
  xp = lapply(xp, na.omit.strip)
  xp = lapply(xp, log2)
  # l_xp = sort(unlapply(xp, length), decreasing = T)

  FigS.All_Imprinted_Genes.PGC_paired = intermingle2lists(xm[names(l_X)], xp[names(l_X)])
  wstripchart(FigS.All_Imprinted_Genes.PGC_paired, jitter = .4, border = "grey", ylb = "log2(Allleic Read Count +1)", tilted_texmt = T, w=14, pch=18, pchcex = 1, col=c(2,4),colorbyColumn = T)
  barplot_label(barplotted_variable = unlapply(FigS.All_Imprinted_Genes.PGC_paired, mean), labels = unlapply(FigS.All_Imprinted_Genes.PGC_paired, length), bottom = T, w=14)

}


Fig.3a.reviewer3_q1b =T
if (Fig.3a.reviewer3_q1b) {

  # x_mat = iPE_matrix
  ls_gene_category =  splititsnames_byValues(gene_category)
  iGenes_ = c(ls_gene_category$Impr_PatExpr, ls_gene_category$Impr_MatExpr)

  x_mat = DP_SUM_parental_HQ[ iGenes_, ControlNamez_HQ]
  x_mat[x_mat==0] = NA


  x = colsplit(t(x_mat+1), f=rownames(x_mat))
  x = lapply(x, na.omit.strip)
  x = lapply(x, log2)
  l_X = sort(unlapply(x, length), decreasing = T)
  l_X = l_X[l_X>2]
  FigS.All_Imprinted_Genes.Soma = x[names(l_X)]
  wstripchart(FigS.All_Imprinted_Genes.Soma, pchcex = 1, jitter = .4, border = "grey", ylb = "log2(Allleic Read Count +1)", tilted_text = T, w=14)
  barplot_label(barplotted_variable = unlapply(FigS.All_Imprinted_Genes.Soma, mean), labels = unlapply(FigS.All_Imprinted_Genes.Soma, length), bottom = T, w=14)

  x_mat_m = DP_mat[ iGenes_, ControlNamez_HQ]
  x_mat_p = DP_pat[ iGenes_, ControlNamez_HQ]


  x = colsplit(t(x_mat+1), f=rownames(x_mat))
  x = lapply(x, na.omit.strip)
  x = lapply(x, log2)

  xm_mat_m = DP_mat[ iGenes_, ControlNamez_HQ]
  xm_mat_m[is.na(x_mat)] = NA

  xm = colsplit(t(xm_mat_m+1), f=rownames(xm_mat_m))
  xm = lapply(xm, na.omit.strip)
  xm = lapply(xm, log2)
  # l_Xm = sort(unlapply(xm, length), decreasing = T)

  xp_mat_p = DP_pat[ iGenes_, ControlNamez_HQ]
  xp_mat_p[is.na(x_mat)] = NA
  xp = colsplit(t(xp_mat_p+1), f=rownames(xp_mat_p))
  xp = lapply(xp, na.omit.strip)
  xp = lapply(xp, log2)
  # l_xp = sort(unlapply(xp, length), decreasing = T)

  FigS.All_Imprinted_Genes_paired.Soma = intermingle2lists(xm[names(l_X)], xp[names(l_X)])
  wstripchart(FigS.All_Imprinted_Genes_paired.Soma, jitter = .4, border = "grey", ylb = "log2(Allleic Read Count +1)", tilted_texmt = T, w=14, pch=18, pchcex = 1, col=c(2,4),colorbyColumn = T)
  barplot_label(barplotted_variable = unlapply(FigS.All_Imprinted_Genes_paired.Soma, mean), labels = unlapply(FigS.All_Imprinted_Genes_paired.Soma, length), bottom = T, w=14)
  oo()
}
