# Protein inference for aLFQ import data frame
ProteinInference <- function(data, ...)  UseMethod("ProteinInference")

ProteinInference.default <- function(data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE, ...) {
	peptide_sequence <- concentration <- peptide_intensity <- NULL

	data <- data.table(data)

	# if the data is on the transition level
	if ("transition_intensity" %in% names(data)) {
		peptide <- peptide_inference.ProteinInference(data, transition_topx = transition_topx, transition_strictness = transition_strictness, transition_summary = transition_summary, consensus_proteins = consensus_proteins, consensus_transitions = consensus_transitions)
	}
	else {
		peptide <- data
	}

	# if the aLFQ import data frame contains any anchor peptides (not proteins!), the concentrations of the according endogenous peptides (and proteins) is inferred
	if (dim(unique(peptide[,c("run_id","protein_id","concentration"), with = FALSE]))[1]!=dim(unique(peptide[,c("run_id","protein_id"), with = FALSE]))[1]) {
		# first step: all peptides with the same sequence as the anchor peptides are selected
	
		setkeyv(peptide,c("run_id","protein_id","peptide_id"))

		conc_peptides<-subset(peptide,peptide_sequence %in% subset(peptide,concentration != "?")$peptide_sequence & peptide_intensity > 0)
		setkeyv(conc_peptides,c("run_id","protein_id","peptide_sequence"))

		# second step: the intensity of the endogenous peptides is divided by the anchor peptide intensity and multiplied by the concentration of the anchor peptide
		conc_peptides<-conc_peptides[, list("concentration"=(peptide_intensity[which(concentration=="?")]/peptide_intensity[which(concentration!="?")])*as.numeric(concentration[which(concentration!="?")])), by=key(conc_peptides)]

		# third step: the endogenous peptides with know assigned concentrations are averaged to compute a mean protein concentration
		conc_proteins<-conc_peptides
		setkeyv(conc_peptides,c("run_id","protein_id"))

		conc_proteins<-conc_proteins[, list("concentration"=mean(concentration)), by=key(conc_proteins)]

		peptide<-merge(subset(peptide,concentration=="?")[,concentration:=NULL],conc_proteins, by=c("run_id","protein_id"), all.x=TRUE)
		peptide$concentration[ is.na(peptide$concentration) ] <- "?"
	}

	# if the data is on the peptide level
	if ("peptide_intensity" %in% names(peptide)) {
		protein <- protein_inference.ProteinInference(peptide, peptide_method = peptide_method, peptide_topx = peptide_topx, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, fasta = fasta, model = model, combine_precursors = combine_precursors, consensus_proteins = consensus_proteins, consensus_peptides = consensus_peptides)
	}
	else {
		protein <- peptide
	}

	# if the data is on the protein level
	if ("protein_intensity" %in% names(protein)) {
		result <- protein
		result$response<-result$protein_intensity
		result<-result[,!(names(result) %in% c("protein_intensity"))]
	}
	else {
		result <- protein
	}

	return(data.frame(result))
}

protein_inference.ProteinInference <- function(data.dt, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, ...) {
	run_id <- protein_id <- peptide_id <- precursor_charge <- omni <- omni_reference <- peptide_intensity <- mean_peptide_intensity <- min_mean_peptide_intensity <- protein_sequence <- protein_sequence_length <- concentration <- NULL

	if (!is.na(fasta)) {
		peptide_sequences.fasta <- trypsin(fasta)
	}

	# consensus filter
	if (consensus_peptides){
		# only use peptide_ids that occur in all run_ids
		setkeyv(data.dt,c("protein_id","peptide_id","precursor_charge"))

		omni.dt <- data.dt[, list("omni" = length(run_id)), by = key(data.dt)]

		data.dt<-data.dt[omni.dt]

		if (consensus_proteins) {
			data.dt<-subset(data.dt, omni == length(unique(data.dt$run_id)))
		}
		else {
			omni_reference.dt <- data.dt[, list("omni_reference" = length(unique(run_id))), by = protein_id]
			data.dt<-data.dt[omni_reference.dt][omni == omni_reference,]
		}
	}

	# should precursors be summed?
	if (combine_precursors) {
		setkeyv(data.dt,c("run_id","protein_id","peptide_id","peptide_sequence","concentration"))
		data.dt<-data.dt[, list("precursor_charge"=0, "peptide_intensity"=sum(peptide_intensity)), by=key(data.dt)]
	}
			
	if (peptide_method == "top") {
		# consensus filter
		if (consensus_peptides){
			setkeyv(data.dt,c("protein_id","peptide_id","precursor_charge"))

			# calculate mean of per peptide_id
			mean_peptides.dt <- data.dt[, list("mean_peptide_intensity"=mean(peptide_intensity)), by=key(data.dt)]

			data.dt<-data.dt[mean_peptides.dt]
			setkeyv(data.dt,c("run_id","protein_id","concentration"))

			# select top consensus peptides
			data.dt<-data.dt[data.dt[, list("min_mean_peptide_intensity"=min(strictnessfilter.ProteinInference(sort(mean_peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]]

			data.dt<-subset(data.dt, mean_peptide_intensity >= min_mean_peptide_intensity)
		}

		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		if (peptide_summary == "mean") {
			data.dt<-data.dt[, list("response"=mean(strictnessfilter.ProteinInference(sort(peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]

		} else if (peptide_summary == "median") {
			data.dt<-data.dt[, list("response"=median(strictnessfilter.ProteinInference(sort(peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]
		} else if (peptide_summary == "sum") {
			data.dt<-data.dt[, list("response"=sum(strictnessfilter.ProteinInference(sort(peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]
		}
	}
	else if (peptide_method == "all") {
		setkeyv(data.dt,c("run_id","protein_id","concentration"))
	
		if (peptide_summary == "mean") {
			data.dt<-data.dt[, list("response"=mean(peptide_intensity)), by=key(data.dt)]
		}
		else if (peptide_summary == "median") {
			data.dt<-data.dt[, list("response"=median(peptide_intensity)), by=key(data.dt)]
		}
		else if (peptide_summary == "sum") {
			data.dt<-data.dt[, list("response"=sum(peptide_intensity)), by=key(data.dt)]
		}
	}
	else if (peptide_method == "iBAQ") {
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-data.dt[, list("response"=sum(peptide_intensity)/length(unlist((lapply(as.list(peptide_sequences.fasta[[protein_id[1]]]),function(X){if(nchar(X)>=6 && nchar(X)<=30){return(X)}else{return(NA)}})))[!is.na(unlist((lapply(as.list(peptide_sequences.fasta[[protein_id[1]]]),function(X){if(nchar(X)>=6 && nchar(X)<=30){return(X)}else{return(NA)}}))))])), by=key(data.dt)]
	}
	else if (peptide_method == "APEX") {
		peptide_sequences.af <- apexFeatures(data.frame("peptide_sequence" = unique(as.vector(unlist(peptide_sequences.fasta[data.dt$protein_id]))), "apex"=NA, stringsAsFactors=FALSE))

        peptide_sequences.apex <- predict(model,peptide_sequences.af)$prediction[,c("peptide_sequence","apex")]

        setkeyv(data.dt,c("run_id","protein_id","concentration"))
        data.dt<-data.dt[, list("response"=sum(peptide_intensity)/sum(peptide_sequences.apex$apex[which(peptide_sequences.apex$peptide_sequence %in% unlist(unname(peptide_sequences.fasta[protein_id[1]])))])), by=key(data.dt)]
	}
	else if (peptide_method == "NSAF") {
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		proteins <- ldply(read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE), function(x) ldply(x))
		names(proteins) <- c("protein_id","protein_sequence")
		proteins<-data.table(proteins)
		setkey(proteins,protein_id)
		proteins<-proteins[,list("protein_sequence_length"=nchar(protein_sequence)), by=key(proteins)]

		data.dt<-merge(data.dt,proteins, by=c("protein_id"))
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-data.dt[, list("peptide_intensity"=sum(peptide_intensity)/protein_sequence_length), by=key(data.dt)]

		setkey(data.dt,"run_id")
		data.dt<-data.dt[, list("protein_id"=protein_id,"concentration"=concentration,"response"=peptide_intensity/sum(peptide_intensity)), by=key(data.dt)]

		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-unique(data.dt)
	}
	
	return(data.dt)
}

peptide_inference.ProteinInference <- function(data.dt, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", consensus_proteins = TRUE, consensus_transitions = TRUE, ...) {
	run_id <- protein_id <- peptide_id <- precursor_charge <- transition_id <- omni <- omni_reference <- transition_intensity <- mean_transition_intensity <- min_mean_transition_intensity <- NULL


	# consensus filter
	if (consensus_transitions){
		# only use transition_ids that occur in all run_ids
		setkeyv(data.dt,c("protein_id","peptide_id","precursor_charge","transition_id"))

		omni.dt <- data.dt[, list("omni" = length(run_id)), by = key(data.dt)]

		data.dt<-data.dt[omni.dt]

		if (consensus_proteins) {
			data.dt<-subset(data.dt, omni == length(unique(data.dt$run_id)))
		}
		else {
			omni_reference.dt <- data.dt[, list("omni_reference" = length(unique(run_id))), by = protein_id]
			data.dt<-data.dt[omni_reference.dt][omni == omni_reference,]
		}

		# calculate mean of per transition_id
		mean_transitions.dt <- data.dt[, list("mean_transition_intensity"=mean(transition_intensity)), by=key(data.dt)]

		data.dt<-data.dt[mean_transitions.dt]
		setkeyv(data.dt,c("run_id","protein_id","peptide_id","precursor_charge"))

		# select top consensus transitions
		data.dt<-data.dt[data.dt[, list("min_mean_transition_intensity"=min(strictnessfilter.ProteinInference(sort(mean_transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]]

		data.dt<-subset(data.dt, mean_transition_intensity >= min_mean_transition_intensity)
	}

	setkeyv(data.dt,c("run_id","protein_id","peptide_id","peptide_sequence","precursor_charge","concentration"))

	if (transition_summary == "mean") {
		data.dt<-data.dt[, list("peptide_intensity"=mean(strictnessfilter.ProteinInference(sort(transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]

	} else if (transition_summary == "median") {
		data.dt<-data.dt[, list("peptide_intensity"=median(strictnessfilter.ProteinInference(sort(transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]
	} else if (transition_summary == "sum") {
		data.dt<-data.dt[, list("peptide_intensity"=sum(strictnessfilter.ProteinInference(sort(transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]
	}

	return(data.dt)
}

strictnessfilter.ProteinInference <- function(data, strictness="loose") {
	if (NA %in% data && strictness == "loose") {
		return(as.vector(na.omit(data)))
	}
	else if (NA %in% data && strictness == "strict") {
		return(as.vector(0))
	}
	else {
		return(as.vector(data))
	}
}

trypsin <- function(fasta, ...)  UseMethod("trypsin")

trypsin.default <- function(fasta, ...) {
	proteins <- read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)

	sequences <- sapply(proteins, strsplit, "(?!P)(?<=[RK])", perl = TRUE)
	return(sequences)
}