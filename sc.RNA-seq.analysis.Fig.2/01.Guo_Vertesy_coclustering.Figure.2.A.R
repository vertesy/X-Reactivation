######################################################################################################
# 01.Guo_Vertesy_coclustering.Figure2.A.R
######################################################################################################
# source ("~/analysis/Scripts_Transcriptome_Analysis/NewTranscriptomeAnalysis/01.CoClustering/01.Guo_Vertesy_coclustering.Figure2.A.R")
rm(list=ls(all.names  = TRUE))


# Functions ------------------------
try (source ('~/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F); try.dev.off()
require(pheatmap)
require(gplots)


# Setup ----------------------------
OutDir = "~/Google_Drive/X_react_Data/RNA/Expression_analysis/Transcriptome_Analysis/01.Guo_Vertesy_coclustering"
OutDirOrig =OutDir
setup_MarkdownReports(OutDir = OutDir, scriptname = "01.Guo_Vertesy_coclustering.Figure2.A.R")

# Metadata  ------------------------------------------------------------------------------------------------------------
input = "FPKM"
fname =  "rsem.FPKM_table_all.tsv"
fname_Guo = "FPKM_Guo_f.rounded.tsv"

metadata_cells =  read.simple.tsv("metadata_cells.tsv"); attach_w_rownames(metadata_cells)
cell_IDs = rownames(metadata_cells)

Weeks_11donors= c("D1" = 9.1, "D2" = 8, "D4" = 10, "D5" = 14.4, "D6" = 4, "D7" = 4.1, "D8" = 8, "D9" = 10, "D10" = 11, "D11" = 17)
NrDonors = l(Weeks_11donors);NrDonors

# Guo
GuoMetadata = read.simple.tsv("FPKM_Guo_f.rounded.metadata.tsv")

# Genes
pgcGenes = read.simple.tsv("20160309_GermCellDevelopmentGenes.tsv")
pgcGeneNames = rownames(pgcGenes)


# Parameters  ------------------------------------------------------------------------------------------------------------
usePearson = F
MinGeneCount = 5000
UsePGCGenes = T
MinEx = 50
MinCount = 5
GermOnly =F
ReadFiles =T
RemoveBatchEffect =F
k_spearman = 4
k_pearson = k_spearman
ManualAnnot = T
SecondPart =F

CorrPerEmbryo = F
NewCorr =T
CorrPerCellType = T

log_settings_MarkDown(usePearson, input, MinGeneCount, UsePGCGenes, MinEx, MinCount, GermOnly, ReadFiles, RemoveBatchEffect, k_pearson, k_spearman, ManualAnnot, SecondPart)


# Read the data  ------------------------------------------------------------------------------------------------------------
if(ReadFiles){
  Expression = read.simple.tsv(fname)
  Expression_Guo = read.simple.tsv(fname_Guo)
  Cell_Guo_Orig = colnames(Expression_Guo)
  sstrsplit(Cell_Guo_Orig, n = 3)
  # write.simple.tsv(round(Expression_Guo, digits = 0), ManualName = "~/Google_Drive/X_react_Data/cout_tables/Guo_Expression/FPKM_Guo_f.rounded.tsv")
}

# inferrednames = orignames = colnames(Expression)

gn_Vertesy = rownames(Expression); l(gn_Vertesy)
gn_Guo = rownames(Expression_Guo); l(gn_Guo)
cell_Guo = colnames(Expression_Guo) = GuoMetadata$SampleName_Guo

CellTypesGuo = table(GuoMetadata$CellType_Guo)
  wpie(CellTypesGuo, percentage = F)
WeeksGuo = table(GuoMetadata$Week_Guo)
names(WeeksGuo) = paste("wk.", names(WeeksGuo) )
  wpie(WeeksGuo, percentage = F)

cell_Vertesy = colnames(Expression) = Sample_Name[colnames(Expression)]

l(union(gn_Guo,gn_Vertesy))
gn_shared = intersect(gn_Guo,gn_Vertesy); l(gn_shared)
# Merge the data  ------------------------------------------------------------------------------------------------------------
FemaleCells = grep (colnames(Expression), pattern = "^D3 |^d3_", value = T, invert = T) # Exclude Male
xHQCells = grep (colnames(Expression), pattern = "LQ$", value = T, invert = T) # Exclude LQ

# grep("D6",colnames(Expression_Guo))
Expression_All =merge_numeric_df_by_rn(x = Expression[gn_shared,intersect(FemaleCells,xHQCells)], y = Expression_Guo[gn_shared,])
llprint("We analyze", ncol(Expression_All), "female gondal cells")
# toClipboard(colnames(Expression_All))

# Cell Filtering  ------------------------------------------------------------------------------------------------------------
MinGeneCount=5000

NrExpressedGenes = colSums(Expression_All>0)
wbarplot(NrExpressedGenes, hline = MinGeneCount, w=14)

HQ_cells = which_names(NrExpressedGenes >= MinGeneCount); Nr_HQ_fem = l(HQ_cells); Nr_HQ_fem
Expression_All = Expression_All[ , HQ_cells]; dim(Expression_All)

LQ_cells = which_names(NrExpressedGenes < MinGeneCount)

Batch = rep(0, Nr_HQ_fem)
  LastVertesyCell = max(which(substr(HQ_cells,start = 1,2)=="D5"))
  Batch[(LastVertesyCell+1):Nr_HQ_fem]=1

HQ_soma = grep("Soma|AD$|GO$|Control$", HQ_cells, value = T); l(HQ_soma)
HQ_germ_cells = grep("Soma|AD$|GO$|Control$", HQ_cells, invert =T, value = T); l(HQ_germ_cells)
  Nr_HQ_germ = l(HQ_germ_cells); Nr_HQ_germ

Control = rep(0, Nr_HQ_fem); names(Control) = HQ_cells; Control[HQ_soma] = 1
  SomaAndGermCells = table(Control); names(SomaAndGermCells) = c( "Germ", "Soma")
  wpie(SomaAndGermCells,both_pc_and_value = T, col=2:3)

if (GermOnly) {   Expression_All = Expression_All[,HQ_germ_cells]; dim(Expression_All) }


if (RemoveBatchEffect) {
  BatchEffectGenes = read.simple.vec("BatchEffectGenes.vec")
  idx_common_fd = setdiff(rownames(Expression_All_Filtered), BatchEffectGenes); l(idx_common_fd)
  Expression_All_Filtered = Expression_All_Filtered[idx_common_fd, ]; dim(Expression_All_Filtered)
}

# Normalize  ------------------------------------------------------------------------------------------------------------
write.simple.tsv(Expression_All)
Expression_All_norm = median_normalize(Expression_All)
write.simple.tsv(Expression_All_norm)

# Gene Filtering  ------------------------------------------------------------------------------------------------------------
gn_HE_idx = rowSums(Expression_All >= MinEx)>=MinCount
gn_HE = which_names(gn_HE_idx)
Expression_All_HE = Expression_All_norm[gn_HE, ]; dim(Expression_All_HE)

FOUNDPGCGENES = intersect(rownames(Expression_All_norm), pgcGeneNames); l(FOUNDPGCGENES)
Expression_All_PGCgenes = Expression_All_norm[FOUNDPGCGENES, ]

if (UsePGCGenes) {
  Expression_All_Filtered = Expression_All_PGCgenes
  llprint("### We filter on", l(FOUNDPGCGENES), "PGC genes.")
} else {
  Expression_All_Filtered = Expression_All_HE
  llprint("### We keep the most highly expressed ",sum(gn_HE_idx),"or", pc_TRUE(gn_HE_idx), " of the genes")
}

barplot(colSums(Expression_All_Filtered))
dim(Expression_All_Filtered)

# Pearson ------------------------------------------------------------------------------------------------------------


if (usePearson) {
  cormethod ="pearson"
  ExprCorr = cor(Expression_All_Filtered, method = cormethod)
  dev.off(); x= pheatmap(ExprCorr)
  pname = kollapse("CoCluser_All_HE_Genes_", cormethod)
  wplot_save_this(plotname = pname, w=20, h=20)

    DisplayOrder = x$tree_row$labels[x$tree_row$order]
    clusterID = cutree(x$tree_row, k = k_pearson)
    gapz = which(!duplicated(clusterID[DisplayOrder])[-1])

    Expression_All_PGC_genes = Expression_All_norm[intersect(rownames(Expression_All_norm), pgcGeneNames), DisplayOrder]
    Expression_All_PGC_genes_log10 = log10(Expression_All_PGC_genes+1)
    pheatmap(Expression_All_PGC_genes_log10, cluster_cols = F, gaps_col = gapz)
    wplot_save_this(plotname = paste0("AllCells_PGC_genes_",cormethod), w = 20, h=15)

}

# Speraman------------------------------------------------------------------------------------------------------------

if (!usePearson) {
  cormethod ="spearman"
  llprint("## We cluster on",cormethod,"correlation")

  ExprCorr = cor(Expression_All_Filtered, method = cormethod)
  dim(Expression_All_Filtered)

  try.dev.off(); x=pheatmap(ExprCorr, clustering_method = "ward.D2",cutree_rows = 4, cutree_cols = 4, treeheight_row = 0)
  pname = kollapse("CoCluser_All_HE_Genes_", cormethod)
  wplot_save_this(plotname = pname, w=20, h=20)

  DisplayOrder = x$tree_row$labels[x$tree_row$order]
  clusterID = cutree(x$tree_row, k = k_spearman)
  gapz = which(!duplicated(clusterID[DisplayOrder])[-1])

  Expression_All_PGC_genes = Expression_All_norm[intersect(rownames(Expression_All_norm), pgcGeneNames), DisplayOrder]
  Expression_All_PGC_genes_log10 = log10(Expression_All_PGC_genes+1)
    try.dev.off()
    pheatmap(Expression_All_PGC_genes_log10, cluster_cols = F, gaps_col = gapz)
    wplot_save_this(plotname = paste0("AllCells_PGC_genes_",cormethod), w = 20, h=15)

}


# UpdateMETADATA ------------------------------------------------------------------------------------------------------------
UpdateAnnot =T

cellNamez = names(clusterID)
clusterID = cutree(x$tree_row, k = 4); clusterID

if (UpdateAnnot) {

  ClusteredCells = splititsnames_byValues(clusterID); ClusteredCells
  clusterID = translate(vec = clusterID, oldvalues = 1:4, newvalues = c(2,1,3,4)) # So that SOMA is one
    ClusterIDnew = rep(NA, Nr_HQ_fem); ClusterIDnew[!Control] = clusterID[HQ_germ_cells];

  clusterNames = c("SOMATIC", "PGC", "LGC", "MGC")
  clusterNames_vec = clusterNames[clusterID]; names(clusterNames_vec) = names(clusterID)
  clusterNames_vec_simplest =  translate(clusterNames_vec, oldvalues = c('LGC','MGC'), newvalues = "AGC")
    CellTypeDistribution = table(clusterNames_vec)
    wpie(CellTypeDistribution, percentage = F)

  AbelsCells_ID = na.omit.strip(value2name_flip(Sample_Name)[names(clusterID)])
  AbelsCells= value2name_flip(AbelsCells_ID)

  # Update ------------------------------------------------
  metadata_cells_CoClustering = metadata_cells

  idx = intersect(rownames(metadata_cells_CoClustering), AbelsCells_ID)

  germz_ = grep(pattern = "GC$", AbelsCells, value = T)
  idx_germ = names(germz_)
  metadata_cells_CoClustering[idx_germ, "cell_type_annot"] = clusterNames_vec[germz_]
  metadata_cells_CoClustering[idx, "cell_type_annot_simple"] = clusterNames_vec[AbelsCells]
  metadata_cells_CoClustering[idx, "cell_type_annot_simplest"] = clusterNames_vec_simplest[AbelsCells]
  metadata_cells_CoClustering$Sample_Name = paste0(substr(metadata_cells_CoClustering$Sample_Name, 1,6), metadata_cells_CoClustering$cell_type_annot)

  # Colors ---
  MetaDataDir = "./metadata2/"

  CellTypeColz = read.simple.tsv.named.vector(MetaDataDir, "ColorTheme.final.tsv")
  cell_type_col		 = CellTypeColz[metadata_cells_CoClustering$"cell_type_annot"]; names (cell_type_col) = cell_IDs

  CellTypeColz_Simple = read.simple.tsv.named.vector(MetaDataDir, "ColorTheme.final.simplest.tsv")
  cell_type_col_simple = CellTypeColz_Simple[metadata_cells_CoClustering$"cell_type_annot_simplest"]; names (cell_type_col_simple) = cell_IDs

  metadata_cells_CoClustering$"cell_type_annot_simplest" = cell_type_col_simple
  metadata_cells_CoClustering$"cell_type_annot"          = cell_type_col

  # Color_Check(value2name_flip(table(cell_type_col)), savefile = T)
  # Color_Check(value2name_flip(table(cell_type_col_simple)), savefile = T)

  write.simple.tsv(metadata_cells_CoClustering)

}

# ------------------------------------------------------------------------------------------------------------------------
stopif(sum(HQ_cells != colnames(ExprCorr)), message = "Unmatched cell names")

replotW_annot = T
if (replotW_annot) {

  InferCelltypeAnnot =T
  if (InferCelltypeAnnot) {
    CellType = clusterNames_vec[HQ_cells]
    l(HQ_cells)
    prefix = stringr::str_split_fixed(HQ_cells, pattern = " ", n=3)
    prefix[,3] =""

    prefix_vec = apply(prefix, 1, paste, collapse=" ")
    DisplayNames  = paste0(prefix_vec, CellType);
  }

  colnames(ExprCorr) = DisplayNames
  names(DisplayNames) = HQ_cells; names(HQ_cells) = DisplayNames

  try.dev.off(); x=pheatmap(ExprCorr, clustering_method = "ward.D2",cutree_rows = 4, cutree_cols = 4, treeheight_row = 0)
  pname = kollapse("CoCluser_All_HE_Genes_", cormethod)
  wplot_save_this(plotname = pname, w=20, h=20)


  # Gene Expression ----------------------------
  DisplayOrder = x$tree_row$labels[x$tree_row$order]
  Expression_All_PGC_genes = Expression_All_norm[intersect(rownames(Expression_All_norm), pgcGeneNames), DisplayOrder]
  Expression_All_PGC_genes_log10 = log10(Expression_All_PGC_genes+1)
  xx= DisplayNames[colnames(Expression_All_PGC_genes)]
  colnames(Expression_All_PGC_genes_log10) = xx

  # Colors ----------------------------

  names(CellType) =     DisplayNames
  CellAnnotation =            as.data.frame(CellType)
  CellAnnotation$Donor   =    Donor = prefix[,1]
  CellAnnotation$Weeks_All =  Weeks_All = Weeks_11donors[CellAnnotation$Donor]
  CellAnnotation$OldID =      HQ_cells
  CellAnnotation$clusterID =  ClusterIDnew
  CellAnnotation$Control =    Control
  CellAnnotation$Batch =      Batch
  CellAnnotation$DonorDisplay =    paste0(Donor[DisplayNames], " (",Weeks_All[DisplayNames], "wk)")
    attach_w_rownames(CellAnnotation)
  CellAnnotation$DonorDisplay =    paste0(Donor[DisplayNames], " (",Weeks_All[DisplayNames], "wk)")


  ShowDis =  c("CellType", "Donor", "Weeks_All")
  GeneClass = pgcGenes; colnames(GeneClass) ="GeneClass"
  GeneClass_col = terrain.colors(l(table(GeneClass))) ;  names(GeneClass_col) =unlist(unique(GeneClass))

  ColorList = list("CellType" = c(CellTypeColz[c( 'MGC', 'PGC', 'LGC')], CellTypeColz_Simple[c( 'SOMATIC')]),
             "Donor" = icolor_categories(prefix[,1]),
             "GeneClass"=GeneClass_col)

  try.dev.off()
  pheatmap(Expression_All_PGC_genes_log10, cluster_cols = F, gaps_col = gapz, annotation_col = CellAnnotation[,ShowDis], annotation_colors = ColorList, annotation_row = GeneClass)
  wplot_save_this(plotname = paste0("Fig1.A.AllCells_PGC_genes_",cormethod), w = 20, h=15)

  # Colors ----------------------------
  colnames(ExprCorr) = DisplayNames
  try.dev.off(); x=pheatmap(ExprCorr, clustering_method = "ward.D2",cutree_rows = 4, cutree_cols = 4, treeheight_row = 0, annotation_col = CellAnnotation[,ShowDis], annotation_colors = ColorList)
  pname = kollapse("FigS1.S.CoCluser_PGC_Genes_", cormethod)
  wplot_save_this(plotname = pname, w=20, h=20)


}


reviewer2_q2 = T
if (reviewer2_q2) {
  DispOrd = c(1, 4:9,2:3)
  x = table(Donor[DisplayNames])[DispOrd]
  nx = names(x)


  Classes = c("PGC", "LGC", "MGC" )
  GC_distr_mat = matrix.fromNames(nx,Classes)
  for (cl in Classes) {
    COI = which_names(CellType[DisplayNames[HQ_germ_cells]] == cl)
    GC_distr_mat[,cl] = table_fixed_categories(  Donor[COI], categories_vec = nx)
  }


  Fig2.F.Table = cbind(
    "Week" = Weeks_11donors[nx],
    "Cells" = table_fixed_categories(sstrsplit(c(LQ_cells, HQ_cells), pattern = " ")[,1], nx),
    "HQ Cells" = x,
    "Somatic" = table_fixed_categories(Donor[DisplayNames[HQ_soma]], categories_vec = nx),
    "GC" = table_fixed_categories(Donor[DisplayNames[HQ_germ_cells]], categories_vec = nx),
    GC_distr_mat
  )
  any_print("D6 has only 2 germ cells, they dont pass the threshold")
  write.simple.tsv(Fig2.F.Table)
}


# --------------------------------------------
# clusterCol_unique = gplots::rich.colors(max(clusterID))
clusterCol =ColorList$CellType[CellType]
names(clusterCol) =names(CellType)
# clusterCol = translate(clusterID, oldvalues = 1:4, newvalues = clusterCol_unique)

DonorCol = ColorList$Donor[prefix[,1]]
names(DonorCol) = names(CellType)

ClusteredCells = splititsnames_byValues(clusterID); ClusteredCells
llprint("## The cells are seaprated as:")
llwrite_list(ClusteredCells)

attach_w_rownames(GuoMetadata)
names(Week_Guo)

UseManualAnnot =F
if (UseManualAnnot) {
  # week_cocl = read.simple.tsv.named.vector("~/analysis/Scripts_Transcriptome_Analysis/NewTranscriptomeAnalysis/01.CoClustering/01.Guo_Vertesy_coclustering.R.Figure.2-3.panel.Metadata.weeks.tsv")
} else {
  hqGuo = SampleName_Guo[SampleName_Guo %in% HQ_germ_cells]
  hqVertesy = Sample_Name[Sample_Name %in% HQ_germ_cells]

  wgg = Week_Guo[names(hqGuo)]
  wgx = paste0(wgg,"wk"); names(wgx) =names(wgg)
  week_cocl = c(wgx, weeks[names(hqVertesy)])
  # names(week_cocl) = c(hqGuo, hqVertesy)

  # DonorsVertesy = c( 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5')
  # wkk = weeks[flip_value2name(donor)[DonorsVertesy]]
  # toClipboard(wkk)
}


# SecondPart ------------------------------------------------------------------------------------------------------------
if (SecondPart){
  # ## Clustering ------------------------------------------------------------------------------------------------------------

  # Clustering methods comparison ------------------------------------------------------------------------------------------------------------
  Agglomeration_Methods = c("ward.D", "ward.D2", "single", "complete", "UPGMA"= "average" , "WPGMA"= "mcquitty", "WPGMC" = "median" , "UPGMC" = "centroid" )

  for (m in 1:l(Agglomeration_Methods)) {
    Linkage_criterion = Agglomeration_Methods[m]
    # Linkage_criterion = "average"
    hc = hclust(d = dist((Expression_All_Filtered),method = "manhattan"), method = Linkage_criterion )
    pname = kollapse("Heatmap_FPKM_",cormethod ,"_", Linkage_criterion)

    dev.off()
    pheatmap(ExprCorr[hc$order, hc$order], cluster_cols = F, main = pname, clustering_method = Linkage_criterion, treeheight_row = 15 )
    wplot_save_this(plotname = pname, w=20, h=15)
    Linkage_criterion = Agglomeration_Methods[m]; print(Linkage_criterion)
  }


  pgcGeneNames = rownames(GermCellDevelopmentGenes)
  PGC_gn_found = intersect(pgcGeneNames, rownames(Expression_All_Filtered)); l(PGC_gn_found)

  for (m in 1:l(Agglomeration_Methods)) {
    Linkage_criterion = Agglomeration_Methods[m]
    # dist_metric= "euclidean"
    dist_metric= "manhattan"
    Expr = t(log10(Expression_All_Filtered+1))
    hc = hclust(d = dist(Expr, method = dist_metric), method = Linkage_criterion )
    pname = kollapse("Heatmap_FPKM_",dist_metric ,"_", Linkage_criterion)

    dev.off()
    pheatmap(t(Expr)[PGC_gn_found, hc$order], cluster_cols = F, main = pname, clustering_method = Linkage_criterion, treeheight_row = 15 )
    wplot_save_this(plotname = pname, w=20, h=15)
    Linkage_criterion = Agglomeration_Methods[m]; print(Linkage_criterion)
  }


  GermCellDevelopmentGenes = read.simple.tsv("~/Google_Drive/X_react_Data/Abelz/Metadata2/Manual_Metadata/20160309_GermCellDevelopmentGenes.tsv")
  Expression_Guo_norm = median_normalize(Expression_Guo)

  gn_PGC= intersect(rownames(Expression_Guo_norm), rownames(GermCellDevelopmentGenes))
  l(gn_PGC)

  Expression_Guo_Genes = Expression_Guo_norm[gn_PGC,]
  Expression_Guo_Genes[is.na(Expression_Guo_Genes)] <- 0
  Expression_Guo_Genes[Expression_Guo_Genes ==0] <- 1
  Expression_Guo_Genes_log10 = log10(Expression_Guo_Genes)

  pheatmap(Expression_Guo_Genes_log10)
  pname="Expression_Guo_Genes_log10"
  wplot_save_this(plotname = pname, w=15, h=15)

} # SecondPart


# Correlation within embryos ---------------------------------------------------------------------------------
HQ_Donors = c( '4.1', '8', '8', '9.1', '10', '10', '14.4', '17')
names(HQ_Donors) = c( 'D7', 'D2', 'D8', 'D1', 'D4', 'D9', 'D5', 'D11')

"You have to run both with T and F"
suffix = if (NewCorr) ".TrWide" else ""
rownames(ExprCorr) = colnames(ExprCorr) = DisplayNames[HQ_cells]

if (NewCorr) {
  DisplayNames_germ = DisplayNames[HQ_germ_cells]
  # xx = median_normalize(Expression_All_HE[, DisplayNames_germ]); dim(xx)
  xx = median_normalize(Expression_All_HE[, HQ_germ_cells]); dim(xx)
  TotalExprCorr = cor(xx, method = cormethod); dim(TotalExprCorr)
  colnames(TotalExprCorr) =rownames(TotalExprCorr) = DisplayNames_germ
}

if (CorrPerEmbryo) {
  sum(Weeks_All != Weeks_11donors[Donor[DisplayNames]])
  Trim1 = which_names(Weeks_All < 10 &  !Control); Trim1
  Trim2 = which_names(Weeks_All >= 10 &  !Control); Trim2

  Donor_ls = splititsnames_byValues(Donor[DisplayNames_germ])
  Corrz = NULL
  for (d in 1:l(HQ_Donors)) {
    cz = Donor_ls[[ names(HQ_Donors)[d] ]]

    cormat = if (NewCorr) TotalExprCorr[cz,cz] else ExprCorr[cz,cz]
    diag(cormat) = NaN
    Corrz[[d]] = na.omit.strip(as.numeric(cormat))
  }
  names(Corrz) = paste0(names(HQ_Donors), " (",HQ_Donors, "wk)")

  # Corrz = reorder.list(Corrz, namesOrdered = names(sort(Weeks_11donors)))

  yll = "Spearman Correlation"
  pname = kollapse(yll,suffix)
  wstripchart(Corrz, cex = .1, colorbyColumn = T, col = ColorList$Donor, tilted_text = T, pchlwd = 5, pchcex = .1, plotname = pname)

  pname = kollapse("Spearman Correlation within embryos",suffix)
  wvioplot_list(Corrz, coll = ColorList$Donor, tilted_text = T, yoffset = -.1, ylb = yll, plotname = pname)

  BeforeAndAfterWeek10 = list(
    UpTo_Wk10 = unlist(Corrz[1:4]),
    From_Wk10 = unlist(Corrz[5:8])
  )

  MWW = wilcox.test(BeforeAndAfterWeek10$"UpTo_Wk10", BeforeAndAfterWeek10$"From_Wk10")
  subb =paste('p-value @ MWW: ', iround(as.numeric(MWW[3]),2) ); #llprint ( subb)

  pname = kollapse("Spearman Correlation Before and After wk10",suffix)
  wstripchart(BeforeAndAfterWeek10, jitter = .4, colorbyColumn = T, col=2:3, pchlwd = 5, pchcex = .1
              , tilted_text = T, ylb = yll, sub = subb, plotname = pname)
  wvioplot_list(BeforeAndAfterWeek10, tilted_text = T, ylb = yll, sub = subb, plotname = pname)

}
# Correlation within classes ---------------------------------------------------------------------------------

if (CorrPerCellType) {
  names(Batch) = names(ClusterIDnew) = DisplayNames

  VerC = which_names(Batch ==0)
  GuoC = which_names(Batch ==1)

  Corrz = NULL
  for (d in 1:3) {

    # cz =DisplayNames[ClusteredCells[[d+1]]]
    cz =ClusteredCells[[d+1]]
    separateBatches=F
    if (separateBatches) {
      cVer = intersect(cz , VerC) # separete the cells : otherwise you measure batch effects
      cGuo = intersect(cz , GuoC)
      cormat_V = if (NewCorr) TotalExprCorr[cVer, cVer] else ExprCorr[cVer, cVer]
      cormat_G = if (NewCorr) TotalExprCorr[cGuo, cGuo] else ExprCorr[cGuo, cGuo]
      diag(cormat_V) = NaN; diag(cormat_G) = NaN;
      Corrz[[d]] = na.omit.strip(c(as.numeric(cormat_V), as.numeric(cormat_G)))
    } else {
      cormat = if (NewCorr) TotalExprCorr[cz,cz] else ExprCorr[cz,cz]
      diag(cormat) = NaN
      Corrz[[d]] = na.omit.strip(as.numeric(cormat))
    }


  }
  names(Corrz)  = c("PGC", "LGC", "MGC")

  yll = "Spearman Correlation"
  pname = kollapse(yll,suffix)
  wstripchart(Corrz, cex = .1, colorbyColumn = T, col = ColorList$CellType, tilted_text = T, pchlwd = 5, pchcex = .1, plotname = pname)

  pname = kollapse("Spearman Correlation per cell type",suffix)
  wvioplot_list(Corrz, coll = ColorList$CellType, tilted_text = T, yoffset = -.1, ylb = yll, plotname = pname)

}

reviewer1_q12 = T
if (reviewer1_q12) {
  GC_types = c("PGC","LGC","MGC" )
  dim(CellAnnotation)
  xx =split(CellAnnotation$CellType, f = CellAnnotation$DonorDisplay)
  yy = lapply(xx, table_fixed_categories, categories_vec= GC_types)

  DispOrder = c('D7 (4.1wk)', 'D2 (8wk)', 'D8 (8wk)', 'D1 (9.1wk)', 'D4 (10wk)', 'D9 (10wk)', 'D10 (11wk)', 'D5 (14.4wk)', 'D11 (17wk)')
  yy = reorder.list(yy, namesOrdered = DispOrder)

  PlotName = "Fig.S1.X.CellType_per_donor"
  pdfA4plot_on(pname = PlotName, cols=5)
  for (i in 1:l(yy)) {
    wpie(unlist(yy[[i]]), col = CellTypeColz[GC_types] ,savefile = F, percentage = F, plotname = names(yy)[i], radius = .75*log10(unlapply(yy, sum)[i]) )#, border =ColorList$Donor[i]
  }
  pdfA4plot_off()
}


write.simple.tsv(CellAnnotation)
