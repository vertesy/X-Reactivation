########################################################################
# Baran Annotation cleanup
########################################################################
#source("~/analysis/Baran_2015_TissueSpecificImprinting/00.Baran_Annotation.R")


rm(list=c("Human_notes", "Imprinted_species", "Imprinting_status_Baran", "Mouse_GeneID","Mouse_Genename"))
try (source ('~/Github_repos/TheCorvinas/R/CodeAndRoll.R'),silent= F)

Baran_Table_S6 = read.simple.tsv("~/Google_Drive/X_react_Data/Reanalysis_with_other_data/Imprinting/Baran_2015_TissueSpecificImprinting/metadata_Baran_2015_GenomeRes/Baran_Table_S6_unduplicated.tsv")
attach_w_rownames(Baran_Table_S6)
table(Imprinting_status_Baran)

AssumedImprinted = which_names(c(Imprinting_status_Baran != "biallelic" & Imprinting_status_Baran != "consistent with biallelic" ))
l(AssumedImprinted)
inline_vec.char(AssumedImprinted)

# BaranImprinted =c( 'AMPD3', 'ANO1', 'CALCR', 'CDKN1C', 'COPG2', 'COPG2IT1', 'CPA4', 'DCN', 'DHCR7', 'DIO3', 'DIRAS3', 'DLGAP2', 'DLK1', 'DLX5', 'FAM50B', 'GABRA5', 'GABRG3', 'GLIS3', 'GNAS-AS1', 'GPR1', 'GRB10', 'H19', 'HM13', 'HYMAI', 'IGF2', 'IGF2-AS', 'IMPACT', 'INPP5F_V2', 'KCNQ1', 'KCNQ1OT1', 'KLF14', 'L3MBTL1', 'LRRTM1', 'MAGEL2', 'MAGI2', 'MCTS1', 'MCTS2P', 'MEG3', 'MEG8', 'MEST', 'MESTIT1', 'MIMT1', 'MIR134', 'MIR184', 'MIR296', 'MIR335', 'MIR371A', 'MIR483', 'MIR675', 'MKRN3', 'NAP1L5', 'NDN', 'NLRP2', 'NNAT', 'NPAP1', 'NTM', 'PEG10', 'PEG3', 'PLAGL1', 'PON2', 'PON3', 'PRIM2', 'PWRN1', 'RB1', 'RBP5', 'RTL1', 'SGCE', 'SGK2', 'SLC22A18', 'SNORD107', 'SNORD108', 'SNORD109A', 'SNORD109B', 'SNORD115-48', 'SNORD116-1', 'SNORD64', 'SNRPN', 'SNURF', 'TCEB3C', 'TH', 'TSPAN32', 'TSSC4', 'UBE3A', 'ZDBF2', 'ZFAT-AS1', 'ZIM2', 'ZIM3', 'ZNF331', 'ZNF597')






