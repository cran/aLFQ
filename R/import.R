# import of mass spectrometry proteomics data analysis software reports.
import <-
function(ms_filenames, ms_filetype, ...) UseMethod("import")

import.default <- function(ms_filenames, ms_filetype, concentration_filename=NA, fasta = NA, targeted_fdr=0.01, spectral_prob=0.95, target_column="ADJNSAF", mergeruns=FALSE, sumruns=FALSE, replace_run_id=FALSE, filtertop=TRUE, removedecoys=TRUE,...) {
	transition_intensity <- peptide_intensity <- protein_intensity <- NULL

	if (!ms_filetype %in% c("skyline","openswath","mprophet","openmslfq","abacus_peptide","abacus_protein")) {
		stop("Please select a valid filetype. Options: \"skyline\", \"openswath\", \"mprophet\", \"openmslfq\", \"abacus_peptide\", \"abacus_protein\"")
	}
	
	# ms_filenames must be provided as vector and are converted to a list
	ms_filenames <- as.list(ms_filenames)
	
	# Skyline import is facilitated by the Transitions Results report
	if (ms_filetype=="skyline") {
		# Skyline specfic adapter
		data.ms <- skyline_converter.import(ms_filenames)
		if (mergeruns==TRUE){
			data.ms <- mergeruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "merged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "summed"
		}
		# Skyline exports all selected transitions. But we only want to use those with an associated MS2 intensity
		data <- subset(data.ms,is.finite(transition_intensity))
	}
	# OpenSWATH import is facilitated by using the output from either mProphet or the TOPPtool OpenSwathFeatureXMLToTSV
	else if (ms_filetype=="openswath") {
		# OpenSWATH specific adapter
		data.ms <- openswath_converter.import(ms_filenames,fdr=targeted_fdr, replace_run_id=replace_run_id, filtertop = filtertop, removedecoys=removedecoys)

		if (mergeruns==TRUE){
			data.ms <- mergeruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "merged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(transition_intensity))
	}
	# mProphet import is facilitated by using the all_peakgroups.xls output
	else if (ms_filetype=="mprophet") {
		# mProphet specific adapter
		data.ms <- mprophet_converter.import(ms_filenames,fdr=targeted_fdr, replace_run_id=replace_run_id, filtertop = filtertop, removedecoys=removedecoys)

		if (mergeruns==TRUE){
			data.ms <- mergeruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "merged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(transition_intensity))
	}
	# OpenMS LFQ import is facilitated from the TOPPtool ProteinQuantifier and the resulting peptides.csv file
	else if (ms_filetype=="openmslfq") {
		# OpenMS LFQ specific adaptor
		data.ms <- openmslfq_converter.import(ms_filenames)
		if (mergeruns==TRUE){
			data.ms <- mergeruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "merged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(peptide_intensity))
	}
	# Abacus peptide import is facilitated from the Abacus with peptide output settings
	else if (ms_filetype=="abacus_peptide") {
		data.ms <- abacus_peptide_converter.import(ms_filenames, prob=spectral_prob, fasta=fasta)
		if (mergeruns==TRUE){
			data.ms <- mergeruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "merged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(peptide_intensity))
	}
	# Abacus protein import is facilitated from the Abacus with default output settings
	else if (ms_filetype=="abacus_protein") {
		data.ms <- abacus_protein_converter.import(ms_filenames, prob=spectral_prob, target_column=target_column)
		if (mergeruns==TRUE){
			data.ms <- mergeruns.import(data.ms,target="protein_intensity")
			data.ms$run_id <- "merged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="protein_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(protein_intensity))
	}

	if ("peptide_sequence" %in% names(data)) {
		# peptide_sequence is stripped and can only contain the 20 natural AA without any modifications
		data$peptide_sequence<-stripsequence.import(data$peptide_sequence)
	}

	# if no concentration file is supplied, "?" is used to indicate absence
	if (!is.na(concentration_filename)) {
		data.conc <- read.csv(concentration_filename, stringsAsFactors=FALSE)
		data.conc$concentration<-as.numeric(data.conc$concentration)
		if ("protein_id" %in% names(data.conc)) {
			# if proteins were manually (e.g. with Skyline quantified)
			if ("run_id" %in% names(data.conc)) {
				data <- merge(data,data.conc, by.x=c("run_id", "protein_id"), by.y=c("run_id", "protein_id"), all.x = TRUE)
			}
			else {
				data <- merge(data,data.conc, by.x="protein_id", by.y="protein_id", all.x = TRUE)
			}
		}
		else if ("peptide_id" %in% names(data.conc)) {
			# if only the spiked-in peptides are known
			if ("run_id" %in% names(data.conc)) {
				data <- merge(data, data.conc, by.x=c("run_id", "peptide_id"), by.y=c("run_id", "peptide_id"), all.x = TRUE)
			}
			else {
				data <- merge(data, data.conc, by.x="peptide_id", by.y="peptide_id", all.x = TRUE)
			}
		}

		data$concentration[ is.na(data$concentration) ] <- "?"
	}
	else {
		data$concentration <- "?"
	}

	return(data)
}

skyline_converter.import <- function(ms_filenames, ...)  {
	data.files = lapply(ms_filenames, read.csv, stringsAsFactors=FALSE)
	data.import <- do.call("rbind", data.files)

	# Skyline headers
	data.ms <- data.import[,c("ReplicateName","ProteinName","PeptideSequence","FragmentIon","PeptideSequence","PrecursorCharge","Area")]
	
	# aLFQ headers
	names(data.ms) <- c("run_id","protein_id","peptide_id","transition_id","peptide_sequence","precursor_charge","transition_intensity")
	
	# replace Skyline NA with R NA values.
	data.ms <- replace(data.ms, data.ms=="#N/A", NA)
	
	data.ms$precursor_charge <- as.numeric(data.ms$precursor_charge)
	data.ms$transition_intensity <- as.numeric(data.ms$transition_intensity)

	return(data.ms)
}

openswath_converter.import <- function(ms_filenames, fdr, replace_run_id, filtertop, removedecoys,...)  {
	transition_code <- transition_code_intensity <- transition_id <- peptide_id <- transition_fragment <- m_score <- peak_group_rank <- decoy <- NULL

	data.files = lapply(ms_filenames, read_table.import, replace_run_id=replace_run_id)
	data.import <- data.table(do.call("rbind", data.files))


	# mProphet headers
	data.ms <- data.import[,c("run_id","ProteinName","FullPeptideName","m_score","peak_group_rank","Sequence","Charge","aggr_Fragment_Annotation","aggr_Peak_Area","decoy"), with = FALSE]

	# aLFQ headers
	setnames(data.ms,c("run_id","ProteinName","FullPeptideName","m_score","peak_group_rank","Sequence","Charge","aggr_Fragment_Annotation","aggr_Peak_Area","decoy"),c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","transition_code","transition_code_intensity","decoy"))
	
	setkeyv(data.ms,c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","decoy"))

	data.ms<-data.ms[, list("transition_fragment"=strsplit(transition_code,";")[[1]],"transition_intensity"=as.numeric(strsplit(transition_code_intensity,";")[[1]])), by=key(data.ms)]

	data.ms[, transition_id:=paste(peptide_id,transition_fragment)]

	# filter FDR
	data.ms <- subset(data.ms, m_score <= fdr)

	# filter top peakgroup
	if (filtertop) {
		data.ms <- subset(data.ms, peak_group_rank == 1)
	}

	# filter decoys
	if (removedecoys) {
		data.ms <- subset(data.ms, decoy == FALSE)
	}

	return(data.frame(data.ms[,c("run_id","protein_id","peptide_id","transition_id","peptide_sequence","precursor_charge","transition_intensity"), with = FALSE]))
}

mprophet_converter.import <- function(ms_filenames, fdr, replace_run_id, filtertop, removedecoys, ...)  {
	transition_code <- transition_code_intensity <- transition_id <- peptide_id <- transition_fragment <- m_score <- peak_group_rank <- decoy <- NULL

	data.files = lapply(ms_filenames, read_table.import)
	data.import <- data.table(do.call("rbind", data.files))

	# mProphet headers
	data.ms <- data.import[,c("run_id","protein","transition_group_record","m_score","peak_group_rank","transition_group_pepseq","transition_group_charge","transition_code_target","abs_area_code_target","decoy"), with = FALSE]

	# aLFQ headers
	setnames(data.ms,c("run_id","protein","transition_group_record","m_score","peak_group_rank","transition_group_pepseq","transition_group_charge","transition_code_target","abs_area_code_target","decoy"),c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","transition_code","transition_code_intensity","decoy"))
	
	setkeyv(data.ms,c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","decoy"))

	data.ms<-data.ms[, list("transition_fragment"=strsplit(transition_code,",")[[1]],"transition_intensity"=as.numeric(strsplit(transition_code_intensity,",")[[1]])), by=key(data.ms)]

	data.ms[, transition_id:=paste(peptide_id,transition_fragment)]

	# filter FDR
	data.ms <- subset(data.ms, m_score <= fdr)

	# filter top peakgroup
	if (filtertop) {
		data.ms <- subset(data.ms, peak_group_rank == 1)
	}

	# filter decoys
	if (removedecoys) {
		data.ms <- subset(data.ms, decoy == FALSE)
	}

	return(data.frame(data.ms[,c("run_id","protein_id","peptide_id","transition_id","peptide_sequence","precursor_charge","transition_intensity"), with = FALSE]))
}

read_table.import <- function(filename, replace_run_id = FALSE, ...) {
	# helper function for mProphet to read runs. filenames are converted to the run_id.
	data <- read.csv(filename, sep = "\t",stringsAsFactors=FALSE)

	if (replace_run_id) {
		data$run_id <- filename
	}

	return(data)
}

openmslfq_converter.import <- function(ms_filenames, ...)  {
	# we need to skip the comments in the header
	data.files = lapply(ms_filenames, read.csv, sep = ",",  comment.char = "#", blank.lines.skip = TRUE, stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,openmslfq_transform.import)

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
openmslfq_transform.import <- function(data, ...) {
	n_proteins <- NULL

	# converts short to long format
	# filter proteotypic peptides
	data <- subset(data, n_proteins == 1)
	
	data.trans <- data.frame(run_id = NA, protein_id = NA, peptide_id = NA, peptide_sequence = NA, precursor_charge = NA, peptide_intensity = NA, stringsAsFactors=FALSE)[0,]

	i <- length(names(data))-1
	while( i > 4 ) {
		data.trans <- rbind(data.trans,data.frame(run_id = names(data)[i], protein_id = data$protein, peptide_id = data$peptide, peptide_sequence = data$peptide, precursor_charge = data$charge, peptide_intensity = data[,i], stringsAsFactors=FALSE))
	
		i <- i-1
	}

	return(data.trans)
}

abacus_peptide_converter.import <- function(ms_filenames, prob, fasta, ...)  {
	if (!is.na(fasta)) {
		peptide_sequences.fasta <- trypsin(fasta)
	}
	else {
    	stop("Specify FASTA file for ABACUS import.")
	}

	data.files = lapply(ms_filenames, read.csv, sep="\t", stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,function(X){abacus_peptide_transform.import(X,prob,fasta,peptide_sequences.fasta)})

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
abacus_peptide_transform.import <- function(data, prob, fasta, peptide_sequences.fasta, ...)  {

	data.re<-ldply(lapply(unlist(strsplit(names(data)[grepl("_MAXPROB", names(data))],"_MAXPROB")),function(X){data.frame("RUN"=X,"MODPEP"=data$MODPEP,"CHARGE"=data$CHARGE,"MAXPROB"=data[,paste(X,"MAXPROB",sep="_")],"NSPECS"=data[,paste(X,"NSPECS",sep="_")])}),data.frame)
	data.filt <- subset(data.re, "MAXPROB" >= prob)

	data.filt$peptide_sequence<-data.filt$MODPEP
	data.filt$peptide_sequence<-stripsequence.import(data.filt$peptide_sequence)
							
	data.ms <- data.filt[,c("RUN","MODPEP","peptide_sequence","CHARGE","NSPECS")]
	names(data.ms) <- c("run_id","peptide_id","peptide_sequence","precursor_charge","peptide_intensity")
	
	names(peptide_sequences.fasta) <- paste(names(peptide_sequences.fasta),"<<remove>>",sep="")
	peptide_sequences.list <- unlist(peptide_sequences.fasta)
	
	peptide_sequences.df <- as.data.frame(sapply(peptide_sequences.list, rbind))
	
	peptide_sequences.df$protein_id <- sapply(strsplit(row.names(peptide_sequences.df), "\\|"), "[[", 1)
	peptide_sequences.df$protein_id_long <- sapply(strsplit(row.names(peptide_sequences.df), "<<remove>>"), "[[", 1)
	names(peptide_sequences.df) <- c("peptide_sequence","protein_id","protein_id_long")
	row.names(peptide_sequences.df) <- NULL
	
	data.ms <- merge(data.ms,peptide_sequences.df)[,c("run_id","protein_id_long","peptide_id","peptide_sequence","precursor_charge","peptide_intensity")]
	names(data.ms) <- c("run_id","protein_id","peptide_id","peptide_sequence","precursor_charge","peptide_intensity")
	return(data.ms)
}

abacus_protein_converter.import <- function(ms_filenames, prob, target_column, ...)  {
	data.files = lapply(ms_filenames, read.csv, sep="\t", stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,function(X){abacus_protein_transform.import(X,prob,target_column)})

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
abacus_protein_transform.import <- function(data, prob, target_column, ...)  {
	PW <- ISFWD <- NULL
	data.re<-ldply(lapply(unlist(strsplit(names(data)[which(substr(names(data),1,3)!="ALL")][grepl("_ID", names(data)[which(substr(names(data),1,3)!="ALL")])],"_ID")),function(X){data.frame("RUN"=X,"PROTID"=data$PROTID,"PROTLEN"=data$PROTLEN,"ISFWD"=data$ISFWD,"PW"=data[,paste(X,"PW",sep="_")],"NUMSPECSTOT"=data[,paste(X,"NUMSPECSTOT",sep="_")],"TOTNSAF"=data[,paste(X,"TOTNSAF",sep="_")],"NUMSPECSUNIQ"=data[,paste(X,"NUMSPECSUNIQ",sep="_")],"UNIQNSAF"=data[,paste(X,"UNIQNSAF",sep="_")],"NUMSPECSADJ"=data[,paste(X,"NUMSPECSADJ",sep="_")],"ADJNSAF"=data[,paste(X,"ADJNSAF",sep="_")])}),data.frame)

	data.filt <- subset(data.re, PW >= prob & ISFWD==1)

	data.ms <- data.filt[,c("RUN","PROTID",target_column)]
							
	names(data.ms) <- c("run_id","protein_id","protein_intensity")
	
	return(data.ms)
}

stripsequence.import <- function(X, ...) {
	return(gsub('[a-z]','',gsub('\\[[^\\[\\]]*\\]','',gsub('\\[[^\\[\\]]*\\]','',gsub('\\([^()]*\\)','',gsub('\\([^()]*\\)','',gsub('\\{[^{}]*\\}','',gsub('\\{[^{}]*\\}','',X)))),perl=TRUE),perl=TRUE),perl=TRUE))


}

mergeruns.import <- function(data,target="transition_intensity", ...)  {
	protein_id <- peptide_id <- transition_id <- peptide_sequence <- precursor_charge <- NULL

	if (target=="transition_intensity") {
		data.merged <- ddply(data, .(protein_id,peptide_id,transition_id,peptide_sequence,precursor_charge), function(X) {data.frame("transition_intensity" = mean(X$transition_intensity, na.rm = TRUE), stringsAsFactors=FALSE)})
	}
	else if (target=="peptide_intensity") {
		data.merged <- ddply(data, .(protein_id,peptide_id,peptide_sequence,precursor_charge), function(X) {data.frame("peptide_intensity" = mean(X$peptide_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	else if (target=="protein_intensity") {
		data.merged <- ddply(data, .(protein_id), function(X) {data.frame("protein_intensity" = mean(X$protein_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	
	return(data.merged)
}

sumruns.import <- function(data,target="transition_intensity", ...)  {
	protein_id <- peptide_id <- transition_id <- peptide_sequence <- precursor_charge <- NULL

	if (target=="transition_intensity") {
		data.summed <- ddply(data, .(protein_id,peptide_id,transition_id,peptide_sequence,precursor_charge), function(X) {data.frame("transition_intensity" = sum(X$transition_intensity, na.rm = TRUE), stringsAsFactors=FALSE)})
	}
	else if (target=="peptide_intensity") {
		data.summed <- ddply(data, .(protein_id,peptide_id,peptide_sequence,precursor_charge), function(X) {data.frame("peptide_intensity" = sum(X$peptide_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	else if (target=="protein_intensity") {
		data.summed <- ddply(data, .(protein_id), function(X) {data.frame("protein_intensity" = sum(X$protein_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	
	return(data.summed)
}