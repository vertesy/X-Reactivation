##########################################################################################
### Imprinted Genes per ICR
##########################################################################################
# source("~/analysis/DevCell_analysis/15.y.Imprinting_ICR.R")


create_set_OutDir(OutDirOrig, "/15.y.Imprinting_ICR") # same as before
setup_logging_markdown("15.y.Imprinting_ICR.R", append = F)


# Setup  ------------------------------------------------------------------------------------
thr_min = 10

# Calculations   ------------------------------------------------------------------------------------
CellSets =  c("ControlNamez_HQ", "GermNamez_HQ")
CS = c("Soma", "Germ")

x=1;
for (x in 1:l(CellSets)) {
	nameCS = CS[x]
	SetOfCells = CellSets[x]
	ICR_PP = PP[ BecomesBiallelic, get(SetOfCells) ]

	llprint ( "We do not observe the following ICR-member genes in our data:", setdiff(BecomesBiallelic, gene_names))

	BecomesBiallelic_ICR_class = ICR_annotation[BecomesBiallelic]
	ICR_classes_bi = split(BecomesBiallelic, f = BecomesBiallelic_ICR_class)

	#### Add Independent iGenes

	ICR_classes_bi[["Independent"]] = unlist(non_ICR_members)

	Silenced_ICR_class = ICR_annotation[Silenced_Genes_in_Erasure]
	ICR_classes_silenced = split(Silenced_Genes_in_Erasure, f = Silenced_ICR_class)


	iPE_matrix = rbind(		PP_matrix[     intersect(gene_names,Impr_MAT), ],
							1 - PP_matrix[ intersect(gene_names,Impr_PAT), ]	)


	DP_candyland = EraSure = EraSure_PE = list (NA)
	j=1
	for ( j in 1:l(ICR_classes_bi) ) {
		DP_observed = unlist(DP_SUM_parental[ ICR_classes_bi[[j]], get(SetOfCells) ])
		DP_candyland[[j]] = DP_observed[which(DP_observed>0)]
		PP_observed =	na.omit.strip(unlist(PP[ ICR_classes_bi[[j]], get(SetOfCells) ]))
		# PP_observed[PP_observed < .5] = 1 - PP_observed[PP_observed < .5] 			# Scale them 1-.5
		EraSure[[j]] = jitter(PP_observed, amount = .025)
			stopifnot(l(PP_observed) ==l(DP_candyland[[j]] ) )
		PE_observed =	na.omit.strip( unlist(iPE_matrix[ ICR_classes_bi[[j]], get(SetOfCells) ]))
			stopifnot(l(PP_observed)==l(DP_candyland[[j]] ))
		EraSure_PE[[j]] = jitter(PE_observed, amount = .025)
	} # for
	names(EraSure) = names(EraSure_PE) = names(DP_candyland) = ICR_short_names
	EraSure = lapply(EraSure, function (x) x* 100)
	EraSure_PE = lapply(EraSure_PE, function (x) x* 100)
	PP_observed = lapply(PP_observed, function (x) x* 100)

	ccc = val2col(log(unlist(DP_candyland)))
		stopifnot(l(ccc) ==l(unlist(DP_candyland)))
	# ccc [unlist(DP_candyland)< thr_min] = 0; ccc
	ccc = as.listalike(ccc, DP_candyland)

	# wstripchart_list(EraSure, bg = ccc, plotname = "EraSure_Expression")
	pname = kollapse("EraSure_Expression_PP_", nameCS, print = F)
	wstripchart_list(EraSure, bg = ccc, plotname = pname, ylab ="% Read from the Paternal allele", jitter = .3)

	pname = kollapse("EraSure_Expression_PE_", nameCS, print = F)
	wstripchart_list(EraSure_PE, bg = ccc, plotname = pname, ylab ="% Read from the expected allele", jitter = .3, ylim =c(0,100))

	NM = lapply(DP_candyland, names)
	NM2 = lapply(NM, function (c) substr(c, 1,8))
	Col_List_Donor = as.numeric.wNames(donor[unlist(NM2)])
	NCC = as.listalike(Col_List_Donor, DP_candyland)
	pname2 = kollapse(pname, "_Donor", print = F)
	wstripchart_list(EraSure_PE, bg = NCC, plotname = pname2, ylab ="% Read from the expected allele", jitter = .35)
	legend("bottomleft", legend =unique(donor), fill = c(1,2,5,3,4), bty="n")
	wplot_save_this(plotname = plotnameLastPlot)

	if (nameCS == "Germ") { # cell type colors on boxplot
		Col_List_Celltype = cell_type_col_simple[unlist(NM2)]
		NCC = as.listalike(Col_List_Celltype, DP_candyland)
		pname3 = kollapse(pname, "_CellTypeCol", print = F)
		wstripchart_list(EraSure_PE, bg = NCC, plotname = pname3, ylab ="% Read from the expected allele", jitter = .35)
	}


} # for "ControlNamez_HQ" or "ControlNamez_HQ"

