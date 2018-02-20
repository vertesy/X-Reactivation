

# Setup -------------------------------------------------------------------------------------------
NrBins_auto = 6
UseErroBars = F
PlottingThrLo =0
# Plot_LowThr_log2 = -1

ChrSplit_DF_Mat= split(DP_mat_clean[HQ_Cellnames], f=chrz_clean)
ChrSplit_DF_Pat= split(DP_pat_clean[HQ_Cellnames], f=chrz_clean)
ChrSplit_DF_Mat$"chrX" = NULL
ChrSplit_DF_Pat$"chrX" = NULL

# Determine active allel -------------------------------------------------------------------------------------------
ChrExpr_Mat = lapply(ChrSplit_DF_Mat, colSums)
ChrExpr_Pat = lapply(ChrSplit_DF_Pat, colSums)

MatIsActive = NULL
for ( d in 1:l(ChrExpr_Pat)) {
	print (names(ChrExpr_Pat)[d])
	MatIsActive[[d]] = unlist(ChrExpr_Mat[[d]]) > unlist(ChrExpr_Pat[[d]])
}; MatIsActive

# convert to active inactive -------------------------------------------------------------------------------------------
ChrSplit_DF_act = ChrSplit_DF_Pat
ChrSplit_DF_inact = ChrSplit_DF_Mat
for ( d in 1:l(ChrExpr_Pat)) {
	ChrSplit_DF_act[[d]][ MatIsActive[[d]] ] = ChrSplit_DF_Mat[[d]][ MatIsActive[[d]] ]   # overwrite Active with Maternal in
	ChrSplit_DF_inact[[d]][ MatIsActive[[d]] ] = ChrSplit_DF_Pat[[d]][ MatIsActive[[d]] ]
}

DP_act = unsplit(ChrSplit_DF_act, f = chrz_clean[chrz_clean != "chrX"])
DP_inact = unsplit(ChrSplit_DF_inact, f = chrz_clean[chrz_clean != "chrX"])
	# plot (unlist(DP_act[chrz_clean == "chr1",]), unlist(DP_inact[chrz_clean == "chr1",]), pch=".")

# Log transform ---------------------------------------------------------------------------------------------------------
DP_act_log2 = log2(DP_act)
DP_inact_log2 = log2(DP_inact)
	# plot (unlist(DP_act_log2[chrz_clean == "chr1",]), unlist(DP_inact_log2[chrz_clean == "chr1",]), pch=".")
	# abline(a=0, b=1)
DP_act_log2[DP_act_log2== -Inf] = -1
DP_inact_log2[DP_inact_log2== -Inf] = -1

# Determine bin boundaries --------------------------------------------------------------------------------------------------
Exp_Auto = unlist(DP_act_log2 + DP_inact_log2)
Exp_Auto = Exp_Auto[Exp_Auto > PlottingThrLo]
slices_auto = 	binner(Exp_Auto, binsize = round(l(Exp_Auto)/NrBins_auto)); 	slices_auto

slices_auto[1] = -1.1
slices_auto[NrBins_auto+1] = ceiling(max(xact+xinact, na.rm = T)); names(slices_auto) = NULL; slices_auto

# sort into bins --------------------------------------------------------------------------------------------------
Exp_Auto = DP_act_log2 + DP_inact_log2
BinnedAuto = NULL
for (i in 1:NrBins_auto) { # i=1
	HP = slices_auto[i]
	LP = slices_auto[i+1]
	BinnedAuto[["Auto Act"]][[i]] = 	DP_act_log2[Exp_Auto > HP & Exp_Auto <= LP];
	BinnedAuto[["Auto Inact"]][[i]] = 	DP_inact_log2[Exp_Auto > HP & Exp_Auto <= LP];
	}
l(BinnedAuto[["Auto Act"]])
dim(Exp_Auto)
dim(DP_inact_log2)

Auto_Act_mean  = unlist(lapply(BinnedAuto[["Auto Act"]], mean))
Auto_Inact_mean  = unlist(lapply(BinnedAuto[["Auto Inact"]], mean))
Auto_Act_sem  = unlist(lapply(BinnedAuto[["Auto Act"]], sem))
Auto_Inact_sem  = unlist(lapply(BinnedAuto[["Auto Inact"]], sem))

if (UseErroBars) {
	xerr_auto = Auto_Act_sem
	yerr_auto = Auto_Inact_sem
	arrows(Auto_Act_mean+xerr_auto, Auto_Inact_mean, Auto_Act_mean-+xerr_auto, Auto_Inact_mean, angle=90, code=3, length=0.1) 		# X dim error bar
	arrows(Auto_Act_mean, Auto_Inact_mean+yerr_auto, Auto_Act_mean,Auto_Inact_mean-yerr_auto, angle=90, code=3, length=0.1) 		# Y dim error bar
} # if (UseErroBars)


