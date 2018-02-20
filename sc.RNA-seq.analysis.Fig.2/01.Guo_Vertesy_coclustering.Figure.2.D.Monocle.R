######################################################################################################
# 01.Guo_Vertesy_coclustering.Figure.2.D.Monocle.R
######################################################################################################
# source ("~/analysis/Scripts_Transcriptome_Analysis/NewTranscriptomeAnalysis/01.CoClustering/01.Guo_Vertesy_coclustering.Figure.2.D.Monocle.R")
try.dev.off()

# source ("~/analysis/Scripts_Transcriptome_Analysis/NewTranscriptomeAnalysis/01.CoClustering/01.Guo_Vertesy_coclustering.Figure2.A.R")


# Functions ------------------------
require("ggplot2")
require("reshape2")
require ("monocle"); # source("https://bioconductor.org/biocLite.R"); biocLite("monocle")


# Setup ----------------------------
OutDir = kollapse(OutDirOrig, "/TrRankings_Monocle_RaceID")
OutDir = "~/Google_Drive/X_react_Data/RNA/Expression_analysis/Transcriptome_Analysis/TrRankings_Monocle.CoClustering"
setup_MarkdownReports(OutDir = OutDir, scriptname = "01.Guo_Vertesy_coclustering.Figure.2.D.Monocle.R")

# Parameters ----------------------------
metadata <- read.simple.tsv("metadata_cells_CoClustering.tsv");
attach_w_rownames(metadata)

Genes_of_Interest <- read.simple.tsv("20160309_GermCellDevelopmentGenes.tsv")

FemaleGerms_HQ = which_names(gender & HQ_by_Transcriptome & cell_type_annot_simple != "Soma"); l(FemaleGerms_HQ)
sum(names(gender) !=  names(cell_type_annot_simple))

CellsAnalyzed= table(cell_type_annot_simple[FemaleGerms_HQ])
MarkDown_Table_writer_NamedVector(CellsAnalyzed)

CellTypeCol = read.simple.tsv.named.vector("ColorTheme.final.tsv")

# Paramters ------------------------
minCellExx = 10
minExx = 50

InferGenes = T # use all genes found significantly diff exp for ordering.
DAZL_exp = F
SecondPart = T
RunLastPart = F
PlotRanking = T
filterHE = F # done previously

log_settings_MarkDown(minCellExx, minExx, InferGenes, DAZL_exp, RunLastPart)

# Go Baby ------------------------------------------------------------------------------------------------

xx = colnames(Expression_All_norm)


TPM_matrix <- round(Expression_All_norm )
stopifnot(as.logical(dim(TPM_matrix)))

germ_cells = grep(xx,pattern = "GC$", value = T)
TPM_matrix <- TPM_matrix[, (germ_cells)]
llprint(NCOL(TPM_matrix)," HQ Female germ cells are analyzed")

# ------------------------------------------------------------------------------------------
# phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
CreateMeta = T
if (CreateMeta) {
  sample_sheet = as.data.frame(stringr::str_split_fixed(germ_cells,pattern = " ", n = 4))
  rownames(sample_sheet) = germ_cells; str(sample_sheet)
  colnames(sample_sheet) = c("donor", "CellNR", "Celltype", "Batch")
  sample_sheet$"cell_type_annot_simple" = clusterNames_vec_simplest[germ_cells]
  sample_sheet$"cell_type_annot_simple2" = as.factor.numeric(clusterNames_vec_simplest[germ_cells])
  sample_sheet$"cell_type_annot" = clusterNames_vec[germ_cells]
  sample_sheet$"Batch" = as.numeric(sample_sheet$"donor" %in% c("D1", "D2", "D3", "D4", "D5"))
}

sample_sheet$Sample = 1
pd <- new("AnnotatedDataFrame", data = sample_sheet)


# ------------------------------------------------------------------------------------------
## featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

gene_names= rownames(TPM_matrix)

if (filterHE) {
  llprint("## Filtering lowly expressed genes")
  GenesAboveThr = rowSums	(TPM_matrix>minExx)>minCellExx
  llprint(pc_TRUE(GenesAboveThr), "of the genes are expressed above",minExx,"TPM in at least",minCellExx ,"cells.
          **", sum(GenesAboveThr), "genes kept for analysis**")
  GenesAboveThr = which_names(GenesAboveThr)
  TPM_matrix <- TPM_matrix[GenesAboveThr, ];  dim(TPM_matrix)
}  else { GenesAboveThr = rownames(TPM_matrix)}

TPM_matrix = as.matrix(TPM_matrix); dim(TPM_matrix)
llogit("*It seems that you have to have 10K rows in your data frame, otherwise the 'newCellDataSet' function collapses.*")

gene_annot = as.data.frame(rownames(TPM_matrix)); rownames (gene_annot) = rownames(TPM_matrix); colnames(gene_annot)="gene_annot"
gene_annot <- gene_annot[GenesAboveThr,, drop=F]

if (sum(row.names(TPM_matrix) != row.names(gene_annot))) {print("Error w row names")}
fd <- new(Class = "AnnotatedDataFrame", data =gene_annot)


# Error Check -----------------------------
dim(TPM_matrix); dim(gene_annot); dim(sample_sheet)
stopifnot(NROW(TPM_matrix) == NROW(gene_annot))
stopifnot(NCOL(TPM_matrix) == NROW(sample_sheet))
stopif (sum(colnames(TPM_matrix) != rownames(sample_sheet)), message = "The CELL names are messed up!")
stopif (sum(rownames(TPM_matrix) != rownames(gene_annot)), message = "The GENE names are messed up!")
l(setdiff(row.names(gene_annot), row.names(TPM_matrix))) == l(setdiff(row.names(TPM_matrix), row.names(gene_annot)))

# create the CellDataSet object -----------------------------

  llprint("newCellDataSet requires 1000 rows at least" )
  my_data <- newCellDataSet( cellData = TPM_matrix, phenoData = pd, featureData = fd )
  HSMM = my_data

  GenesExpressedHi = rowSums	(TPM_matrix>minExx)>minCellExx
  llprint(pc_TRUE(GenesExpressedHi), "of the genes are expressed above",minExx,"TPM in at least",minCellExx ,"cells.
          **", sum(GenesExpressedHi), "genes kept for analysis**")

  # QC -----------------------------

  # LogNorm, -----------------------------
  # Once you've excluded cells that do not pass your quality control alters, you should verify that the expression values
  # stored in your CellDataSet follow a distribution that is roughly lognormal:

  llprint("## QC Expression follows Lognormal distribution")
  llprint(" - Why do we test? What does it mean?")
  llprint(" - Nevertheless it looks exaclty as it looks in the Bioconductor Vignette (same spike to the right)")

  # Log-transform each value in the expression matrix.
  L <- log(exprs(HSMM))
  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily"
  melted_dens_df <- melt(t(scale(t(L))))

  # Plot the distribution of the standardized gene expression values.
  plotname = "Genes Expression is roughly lognormal"
  qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') +
    ggtitle(plotname)+
    xlab("Standardized log(TPM)") +
    ylab("Density")
  wplot_save_this(plotname, mdlink = T)

  # 4 Basic differential expression analysis -----------------------------
  llprint("## Differential Gene Expression Analysis")

  marker_genes =as.vector(rownames(Genes_of_Interest)); l(marker_genes)
  llprint("We found", pc_TRUE(marker_genes %in% rownames(TPM_matrix)),"of the marker genes in our data set, we badly miss:", setdiff(marker_genes, rownames(TPM_matrix)), "These were lowly expressed." )

  m = marker_genes[marker_genes %in% rownames(TPM_matrix)]
  GenesExpressedHi = which_names(GenesExpressedHi )


  HSMM_x = (HSMM[GenesExpressedHi,])
  # if(!InferGenes) {  HSMM_x = (HSMM[m,]);dim(HSMM_x)
  # } else {          HSMM_x = (HSMM[GenesExpressedHi,]);dim(HSMM_x)  }

  llprint("### Differntial Expression was calculated on ", NROW(HSMM_x), "highly expressed genes")
  diff_test_res2 <- differentialGeneTest(cds = HSMM_x,  fullModelFormulaStr="~cell_type_annot", cores = 8)
  diff_test_res = diff_test_res2


if (SecondPart) {

  # save.image(file = "Monocle.Rdata")
  # load(file = "~/Google_Drive/X_react_Data/RNA/Expression_analysis/Transcriptome_Analysis/TrRankings_Monocle_RaceID/Monocle.Rdata")

  # Select genes that are significant at an FDR < 10%
  InferGenes=T
  if (InferGenes) {
    sig_genes <- subset(diff_test_res, qval < 0.01); l(sig_genes[,1]); NROW(sig_genes)
    sig_genes <- subset(sig_genes, pval < 0.01); l(sig_genes[,1]); NROW(sig_genes)
  }  else {
    sig_genes <- subset(diff_test_res, qval < 0.1);
  }
  stopif( nrow(sig_genes)==0 , "No genes are significant")

  llprint(dim(sig_genes)[1], "genes are significant at an FDR < 10%. Highest p-value is:", iround(max(sig_genes$pval)))
  # sig_genes[order(sig_genes$pval),]

  # Attach the HUGO symbols and other featureData for these genes
  sig_genes <- merge(fData(HSMM), sig_genes, by="row.names");
  colnames(sig_genes)[1] ="gene_short_name"

  MostSignGenes = head(sig_genes[order(sig_genes$pval),], 18)

  gene_short_name = sig_genes$gene_short_name
  PlotTheseGenes = MostSignGenes$gene_short_name
  Most_Diff_Expr_Genes <- HSMM[PlotTheseGenes,]

  plot_genes_jitter(Most_Diff_Expr_Genes, grouping="cell_type_annot_simple", ncol=3, cell_size = 2, color_by ='cell_type_annot')
  wplot_save_this(plotname = "MostSignGenes",w= 8.3 ,h=11.7)


  # 5 Ordering cells by progress ---------------------------------------------
  llprint("## Ordering cells by progress")

  if(InferGenes) {
    ordering_genes <- row.names(subset(diff_test_res, qval < 0.01)); l(ordering_genes)
  } else {
    ordering_genes <- rownames(Genes_of_Interest)
  }

  # Only use genes are detectably expressed in a sufficient number of cells
  ordering_genes <- 	intersect(ordering_genes, GenesExpressedHi);ordering_genes
  llprint("We order based on", l(ordering_genes), "marker genes that are above the expression criteria (above).")
  llprint(head(ordering_genes))

  GenesNotFound = sort(setdiff(ordering_genes, gene_short_name))
  llprint("The folowing genes were not used for ordering", GenesNotFound)

  HSMM_orig = HSMM
  HSMM <- setOrderingFilter(HSMM, ordering_genes)

  HSMM <- reduceDimension(cds = HSMM, verbose=F)
  try(HSMM <- orderCells(HSMM, num_paths = 1))

  plot_cell_trajectory(HSMM, color_by = "cell_type_annot", show_cell_names = T, cell_name_size = 2)
  wplot_save_this("FigS2.minimal_spanning_tree_names", mdlink = T)

  plot_spanning_tree(HSMM, color_by = "cell_type_annot")+ aes(size =2)
  wplot_save_this("FigS2.minimal_spanning_tree", mdlink = T)

  # plot per donor
  plot_spanning_tree(HSMM, color_by = "donor", backbone_color = "#B99C19")+ aes(size =2)
  wplot_save_this("FigS2.minimal_spanning_tree_per_Donor", mdlink = T)


  # 7 Write out the table ------------------------------------------------------------------------------------------
  colnames(HSMM@phenoData@data)

  ColumnsToWriteOut = c( "Pseudotime", "donor", "cell_type_annot_simple", "cell_type_annot")
  RankTable_Monocle = HSMM@phenoData@data [, ColumnsToWriteOut]
  RankTable_Monocle = RankTable_Monocle[ order(RankTable_Monocle[,1]),] 			# Sort the df

  nr_celllz = dim(RankTable_Monocle)[1]
  RankTable_Monocle = cbind(RankTable_Monocle,
                            "Rank" = 1:nr_celllz 	)

  write.simple.tsv(RankTable_Monocle)

  # Fig.1e Ranking ------------------------------------------------------------------------------------------
  if (PlotRanking) {
    Nr = dim(RankTable_Monocle)[[1]]

    ccc = ColorList$Donor[sample_sheet[rownames(RankTable_Monocle), "donor"]]
    plot(rep(1,Nr), pch=19, cex=1, col=ccc, ylim = c(.9, 1.1), xlab = "", ylab = "", axes = F, main="Monocle ranking order")

    ccc = (ColorList$CellType[clusterNames_vec[rownames(RankTable_Monocle)]])
    points(rep(.95,Nr), pch=15, cex=1, col=ccc)

    ll = ColorList$CellType[-4]
    legend("topright", legend = names(ll),  fill = ll, bty = "n")

    ll = ColorList$Donor
    legend("topleft", legend = names(ll),  fill = ll, bty = "n")
    wplot_save_this(plotname = "Fig1_e.Monocle_Ranking")
  }

  # ------------------------------------------------------------------------------------------


} # SecondPart


