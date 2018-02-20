#####################################################################################################
# 02.Methylation_Analysis.Zecher_plus_controls
######################################################################################################
# source("~/x_reactivation/analysis/Methylation_analysis/02.Methylation_Analysis.Zecher_plus_controls.R")
# rm(list=ls(all = TRUE))
# try(dev.off())


# Functions ------------------------------------------------------------------------------------------------------------------------
require (stringr)
try (source ('~/TheCorvinas/R/CodeAndRoll.R'),silent= F)
require(pheatmap)

# Setup ------------------------------------------------------------------------------------------------------------------------
setup_MarkdownReports("~/Google_Drive/X_react_Data/Guo_2015_Methyation/02.Methylation_Analysis.Zecher_plus_controls", fname = "02.Methylation_Analysis.Zecher_plus_controls.R", append = F)

InputDir = "~/Google_Drive/X_react_Data/Guo_2015_Methyation/Old_stuff/Regions_ENSEMBL/FIle_Chunks/"



# Metadata ------------------------------------------------------------------------------------------------------------------------
regions_of_int = c( 'chr11_2021070_2021302_H19_mC', 'chr11_2721173_2721297_LIT1_mC', 'chr14_101292152_101292376_MEG3_mC', 'chr15_25200010_25200249_SNRPN_mC', 'chr19_57351942_57352097_PEG3_mC', 'chr5_112073373_112073568_APC_mC', 'chr5_55029104_55029220_DDX4_mC', 'chr6_74063525_74063669_DPPA5_mC')
ShortNames = c( 'H19', 'KCNQ1', 'MEG3', 'SNRPN', 'PEG3', 'APC', 'DDX4', 'DPPA5')

fPGC =c( 'PGC_10W_embryo1_F_rep1', 'PGC_10W_embryo1_F_rep2', 'PGC_11W_embryo1_F', 'PGC_17W_embryo1_F')
soma =c( 'Brain_5W_embryo1', 'Heart_5W_embryo1', 'Soma_10W_embryo1_M', 'Soma_11W_embryo1_M', 'Soma_17W_embryo1_F', 'Soma_19W_embryo1_M', 'Soma_7W_embryo1_M', 'Soma_7W_embryo2_M')


# Read in the files & Extract the regions of interest ------------------------------------------------------------------------------------------------------------------------
Zechner_plus_Controls = read.simple.tsv ("mC.per_region.Zechner.tsv")
# View(Zechner_plus_Controls)

index1 = seq(from = 1, to = l(Zechner_plus_Controls), by = 3)
index2 = index1 + 1
index3 = index2 + 1

mC = Zechner_plus_Controls[, index1]
C = Zechner_plus_Controls[, index2]
Total = Zechner_plus_Controls[, index3]
pc_mC= mC/ Total



# Barplot ----------------------------------------------------------------------------------------
pc_mC_F_PGCs = pc_mC[fPGC,regions_of_int]
pc_mC_soma = pc_mC[soma,regions_of_int]

Average_Methylation_rate_in_F_PGCs = 100*colMeans(pc_mC_F_PGCs)
Average_Methylation_rate_in_soma = 100*colMeans(pc_mC_soma)
names(Average_Methylation_rate_in_F_PGCs) = names(Average_Methylation_rate_in_soma) = ShortNames
Average_Methylation_rates = unlist(intermingle2lists(Average_Methylation_rate_in_F_PGCs, Average_Methylation_rate_in_soma))


#sem
pgc_sem = apply(100*pc_mC_F_PGCs, 2, sem)
soma_sem = apply(100*pc_mC_soma, 2, sem)
semz =unlist(intermingle2lists(pgc_sem, soma_sem))

wbarplot(Average_Methylation_rates, errorbar = T, upper = semz, tilted_text = T,  ylab = "% mC", col = c("gold1", "chartreuse4"), mdlink = T)

llprint("### Somatic samples:")
llogit(soma)
llprint("### Female PGC samples:")
llogit(fPGC)

llprint("### Regions:")
llogit(regions_of_int)

# Save data ----------------------------------------------------------------------------------------
write.simple.tsv(mC)
write.simple.tsv(C)
write.simple.tsv(Total)
write.simple.tsv(pc_mC)
