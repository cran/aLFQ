ALF <- function(data, ...) UseMethod("ALF")
# This function implements the workflow used by (1) to select a suitable model for absolute label-free quantification.

# 1.	Ludwig, C., Claassen, M., Schmidt, A. & Aebersold, R. Estimation of Absolute Protein Quantities of Unlabeled Samples by Selected Reaction Monitoring Mass Spectrometry. Molecular \& Cellular Proteomics 11, M111.013987–M111.013987 (2012).

ALF.default <- function(data, report_filename="ALF_report.pdf", prediction_filename="ALF_prediction.csv", peptide_method = "top", peptide_topx = c(1,2,3,4,5,6), peptide_strictness = "loose", peptide_summary = "sum", transition_topx= c(1,2,3,4,5,6), transition_strictness = "loose", transition_summary = "sum", cval_method = "boot" ,cval_mcx = 1000, combine_precursors = TRUE, ...) {
	pdf(file=report_filename)
	
	# nr_peptides nr_transitions tuning
	data.tune <- tune.ALF(data, peptide_method = peptide_method, peptide_topx = peptide_topx, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, cval_method = cval_method, cval_mcx = cval_mcx, transition_topx = transition_topx, transition_strictness = transition_strictness, transition_summary = transition_summary, combine_precursors = combine_precursors)
	transition_topx_min <- which(data.tune == min(data.tune), arr.ind = TRUE)[1]
	peptide_topx_min <- which(data.tune == min(data.tune), arr.ind = TRUE)[2]

	performanceplot.ALF(data.tune)
		
	# calculate optimal model
	optimal.AbsoluteQuantification <- AbsoluteQuantification(data, peptide_method = peptide_method, peptide_topx = peptide_topx_min, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, transition_topx = transition_topx_min, transition_strictness = transition_strictness, transition_summary = transition_summary, combine_precursors = combine_precursors)		
	optimal.AbsoluteQuantification <- predict(optimal.AbsoluteQuantification)
	plot(optimal.AbsoluteQuantification)
	
	optimal.AbsoluteQuantification.cval <- cval.AbsoluteQuantification(optimal.AbsoluteQuantification,method=cval_method, mcx = cval_mcx)
	plot(optimal.AbsoluteQuantification.cval)
	hist.AbsoluteQuantification(optimal.AbsoluteQuantification.cval)
	
	dev.off()
	
	export.AbsoluteQuantification(optimal.AbsoluteQuantification, file = prediction_filename)
}

tune.ALF <- function(data, peptide_method = "top", peptide_topx = c(1,2,3,4,5,6), peptide_strictness = "loose", peptide_summary = "sum", transition_topx= c(1,2,3,4,5,6), transition_strictness = "loose", transition_summary = "sum", cval_method = "boot", cval_mcx = 1000, combine_precursors = "TRUE", ...) {
	cvmfe.mx <- matrix(nrow = length(transition_topx), ncol = length(peptide_topx))
	rownames(cvmfe.mx) <- transition_topx
	colnames(cvmfe.mx) <- peptide_topx

	i <- 1
	while (i <= length(peptide_topx)) {
		j <- 1
		while (j <= length(transition_topx)) {
			cvmfe.mx[j,i] <- cval.AbsoluteQuantification(AbsoluteQuantification(data, peptide_method = peptide_method, peptide_topx = peptide_topx[i], peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, transition_topx = transition_topx[j], transition_strictness = transition_strictness, transition_summary = transition_summary), method = cval_method, mcx = 1000, combine_precursors = combine_precursors)$cv$mfe
			j <- j + 1
		}
		i <- i + 1
	}
		
	return(cvmfe.mx)
}

performanceplot.ALF <- function(x, ...) {
	require(lattice)
	transition_topx_min <- which(x == min(x), arr.ind = TRUE)[1]
	peptide_topx_min <- which(x == min(x), arr.ind = TRUE)[2]
	print(levelplot(t(x),xlab="Peptides",ylab="Transitions",main=paste("Optimal model: ",transition_topx_min," Transitions, ",peptide_topx_min," Peptides",sep="")))
}
