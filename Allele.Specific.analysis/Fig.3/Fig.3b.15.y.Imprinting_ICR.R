##########################################################################################
### Imprinted Genes per ICR
##########################################################################################
# source("~/analysis/DevCell_analysis/15.y.Imprinting_ICR.R")
try(dev.off())

create_set_OutDir(OutDirOrig, "/15.y.Imprinting_ICR") # same as before
setup_logging_markdown("15.y.Imprinting_ICR.R", append = F)


# Setup  ------------------------------------------------------------------------------------
thr_min = 10
TheOrder = c("H19", "Meg3", "Kcnq1", "snurf", "Peg3", "nonICR") # Pat (2) and Mat(3) methylated ICRs, in chromosomal order.

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
	names(ICR_classes_bi) = str_split_fixed(names(ICR_classes_bi), pattern = "_", n=3)[ ,1]
	ICR_classes_bi = ICR_classes_bi[TheOrder[1:5]]

	#### Add Independent iGenes

	ICR_classes_bi[["Independent"]] = unlist(non_ICR_members)

	Silenced_ICR_class = ICR_annotation[Silenced_Genes_in_Erasure]
	ICR_classes_silenced = split(Silenced_Genes_in_Erasure, f = Silenced_ICR_class)

# 	iPE_matrix = rbind(		PP_matrix[     intersect(gene_names,Impr_MAT), ],
# 							1 - PP_matrix[ intersect(gene_names,Impr_PAT), ]	)

	DP_candyland = EraSure = EraSure_PE = list ("H19" =NA)
	j=2
	names(ICR_classes_bi)
	for ( j in 1:l(ICR_classes_bi) ) {
		ICR_j =TheOrder[j]; any_print(ICR_j)
		DP_observed = unlist(DP_SUM_parental[ ICR_classes_bi[[j]], get(SetOfCells) ])
		DP_candyland[[ICR_j]] = DP_observed[which(DP_observed>0)]
		PP_observed =	na.omit.strip(unlist(PP[ ICR_classes_bi[[j]], get(SetOfCells) ]))
		# PP_observed[PP_observed < .5] = 1 - PP_observed[PP_observed < .5] 			# Scale them 1-.5
		EraSure[[ICR_j]] = jitter(PP_observed, amount = .025)
			stopifnot(l(PP_observed) ==l(DP_candyland[[ICR_j]] ) )
		PE_observed =	na.omit.strip( unlist(iPE_matrix[ ICR_classes_bi[[j]], get(SetOfCells) ]))
			stopifnot(l(PP_observed)==l(DP_candyland[[ICR_j]] ))
		EraSure_PE[[ICR_j]] = jitter(PE_observed, amount = .025)
	} # for
	# names(EraSure) = names(EraSure_PE) = names(DP_candyland) = ICR_short_names
	EraSure = lapply(EraSure, function (x) x* 100)
	EraSure_PE = lapply(EraSure_PE, function (x) x* 100)
	PP_observed = lapply(PP_observed, function (x) x* 100)

	ccc = val2col(log2(unlist(DP_candyland)))
		stopifnot(l(ccc) ==l(unlist(DP_candyland)))
	# ccc [unlist(DP_candyland)< thr_min] = 0; ccc
	ccc = as.listalike(ccc, DP_candyland)

	# wstripchart_list(EraSure, bg = ccc, plotname = "EraSure_Expression")
	pname = kollapse("EraSure_Expression_PP_", nameCS, print = F)
	wstripchart_list(EraSure, bg = ccc, plotname = pname, ylab ="% readss from the paternal allele", jitter = .3)

	if (x==1) { 		llprint("## Clean Monoallelic expression from the *expected* allele in all imprinting clusters (ICs), but not in the non-IC imprinted genes") }
	else if (x==2) { 	llprint("## Some but not all ICs show biallelic expression. Non-IC genes might not even be imprinted given the results in soma.") }

	pname = kollapse("EraSure_Expression_PE_", nameCS, print = F)
	wstripchart_list(EraSure_PE[TheOrder], bg = ccc, plotname = pname, ylab ="% readss from the expected allele", jitter = .3, ylim =c(0,100), mdlink = T) # , pchlwd = 0

	ICR_PE_statistics= cbind(
			"median" = unlist(lapply(EraSure_PE[TheOrder], median)),
			"mean" = unlist(lapply(EraSure_PE[TheOrder], mean)),
			"SD" = unlist(lapply(EraSure_PE[TheOrder], sd)),
			"SEM" = unlist(lapply(EraSure_PE[TheOrder], sem))
			)
	write.simple.tsv(ICR_PE_statistics)


	llprint("## Celltype color")
	NM = lapply(DP_candyland, names)
	NM2 = lapply(NM, function (c) substr(c, 1,8))
	NC = as.numeric.wNames(donor[unlist(NM2)])
	NCC = as.listalike(NC, DP_candyland)
	pname = paste0 (pname, "_DonorColor")
	wstripchart_list(EraSure_PE[TheOrder], bg = NCC, plotname = pname, ylab ="% reads from the expected allele", jitter = .35, mdlink = T)
	if (nameCS == "Soma") {
  	Labelz = c("9.1wk (d1)",
  	           "8wk (d2)",
  	           "16.4wk male (d3)",
  	           "10wk (d4)",
  	           "14.4 (d5)")
	  legend("bottomleft", legend = Labelz, fill = c(1,2,5,3,4), inset = .02, bty = "n")
	}


	if (nameCS == "Germ") { # cell type colors on boxplot
		llprint("## There is no clear correlation between allelic bias and cell state either.")
		NM = lapply(DP_candyland, names)
		NM2 = lapply(NM, function (c) substr(c, 1,8))
		NC = cell_type_col_simple[unlist(NM2)]
		NCC = as.listalike(NC, DP_candyland)
		pname = kollapse("EraSure_Expression_PE_Germ", "_CellTypeCol", print = F)
		wstripchart_list(EraSure_PE[TheOrder], bg = NCC, plotname = pname, ylab ="% reads from the expected allele", jitter = .35, mdlink = T)
		}


} # for "ControlNamez_HQ" or "ControlNamez_HQ"
