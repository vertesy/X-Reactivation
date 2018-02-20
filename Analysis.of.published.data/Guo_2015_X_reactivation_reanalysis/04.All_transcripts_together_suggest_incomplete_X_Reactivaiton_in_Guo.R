########################################################################################################################
# 04.All_transcripts_together_suggest_incomplete_X_Reactivaiton_in_Guo.R
########################################################################################################################

require(MarkdownReports)

# FUNCITONS ----------------------------------------
symmetrfy <- function(data, threshld =.5) {
  data = as.matrix(data)
  idx = which(unlist(data)<threshld)
  data[idx] = 1-data[idx]
  return(data)
}


# PARAMTERS ----------------------------------------
ExludeAllEsc =T
AverageSNPsInTheSameGene = T
SymmterfyValuesAndSave = F

datadir = "~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Guo_X_reactivation/data/"
InputDir = if(AverageSNPsInTheSameGene) paste0(datadir,"PerGene/") else paste0(datadir,"PerSNP/")

# SETUP ----------------------------------------
OutDir = "~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Guo_X_reactivation/04.X_reactivation_is_incomplete_in_Guo"
geniez = if (ExludeAllEsc) "Proper-X genes" else "Proper-X and Non-confirmed escapee genes"
TTL = paste0("Median Monoallelic Rate in ",geniez," Suggests Incomplete X-Reactivation in the Guo et al. Data")

setup_MarkdownReports(OutDir = OutDir, scriptname = "04.All_transcripts_together_suggest_incomplete_X_Reactivaiton_in_Guo.R", title = TTL)

log_settings_MarkDown( ExludeAllEsc, AverageSNPsInTheSameGene, SymmterfyValuesAndSave)
llprint("We used",geniez,"for this analysis")

# Metdata ----------------------------------------

ESC = read.simple.vec("./metadata/X_Escape_intersect_always_esc.vec")
HET = read.simple.vec("./metadata/X_Escape_intersect_HeterogenousData.vec")
All_esc = c(ESC, HET)

# Read in Files ----------------------------------------

ls_F = sort(list.files(InputDir, pattern = "^X_ratios_w.+.tsv")); ls_F
sort(ls_F)
nrf = l(ls_F)
trunk = stringr::str_split_fixed(ls_F, pattern = "\\.", n=2)[,1] # get the week put of the string
week = substr(trunk, 10,15); print(week)

Ls_emb= NULL


# Reconstruct Haplotypes ----------------------------------------

FractionNonReactivated = numeric(5); names(FractionNonReactivated)=week
NrGenes = FractionNonReactivated
f=4
for (f in 1:nrf) {
  print(f)
  Ls_emb[[f]] = read.simple.tsv(InputDir,ls_F[f])
  X_ratio = Ls_emb[[f]][, -(2:3)]

  PGC_columns = grep("PGC",colnames(X_ratio))

  gz = if (ExludeAllEsc) setdiff(X_ratio$Gene, All_esc) else setdiff(X_ratio$Gene, ESC)

  if (l(gz)) {
    TheseRows = which(X_ratio$Gene %in% gz)
    NrGenes[f] = l(TheseRows)
    X_expression =X_ratio[TheseRows, PGC_columns]
    X_expression = rbind(symmetrfy(X_expression))
    MedianBias = colMedians(X_expression, na.rm = T)
    disp_order = names(sort(MedianBias, decreasing = F))
    X_expression = X_expression[,disp_order, drop=F]
    Inactivated_Cells = (MedianBias >= .95)
      ccc= sort(!Inactivated_Cells, decreasing = T)+2
    FractionNonReactivated[f] = 100*pc_TRUE(na.omit.strip(Inactivated_Cells), percentify = F)
    PC_Of_Inactivated_Cells = pc_TRUE(na.omit.strip(Inactivated_Cells))
      colnames(X_expression) = stringr::str_split_fixed(colnames(X_expression), n=6, pattern = "_")[,5]
    pname = paste("Allelic Expression in", week[f], "PGCs")
    sbb = paste(sum(Inactivated_Cells, na.rm = T), "or", PC_Of_Inactivated_Cells, "of the PGCs show XaXi expression pattern.")
    plotlist = if(min(dim(X_expression))<2) split(100*X_expression, f=1:l(X_expression)) else   splitByCol(100*X_expression, f = 1:ncol(X_expression))
    names(plotlist) = colnames(X_expression)
    wstripchart(plotlist, plotname = pname, main=pname, col=ccc,colorbyColumn=T, ylim=c(0,100), jitter=.3,
                ylb = "Bias Towards Either Allele, per gene (%)", sub =sbb, tilted_text = T, pchcex = 1)
    abline(h=95, lty=2, col=1)
    wplot_save_this(plotname = plotnameLastPlot , mdlink = T, w=14)
  }
}

