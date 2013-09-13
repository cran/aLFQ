# Protein inference for aLFQ import data frame
ProteinInference <- function(data, ...)  UseMethod("ProteinInference")

ProteinInference.default <- function(data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, ...) {
	peptide_sequence <- concentration <- peptide_intensity <- run_id <- protein_id <- concentration <- NULL

	# if the data is on the transition level
	if ("transition_intensity" %in% names(data)) {
		peptide <- peptide_inference.ProteinInference(data, transition_topx = transition_topx, transition_strictness = transition_strictness, transition_summary = transition_summary)
	}
	else {
		peptide <- data
	}

	# if the aLFQ import data frame contains any anchor peptides (not proteins!), the concentrations of the according endogenous peptides (and proteins) is inferred
	if (dim(unique(peptide[,c("protein_id","concentration")]))[1]!=length(unique(peptide[,c("protein_id")]))) {
		# first step: all peptides with the same sequence as the anchor peptides are selected
		# second step: the intensity of the endogenous peptides is divided by the anchor peptide intensity and multiplied by the concentration of the anchor peptide

		conc_peptides<-ddply(subset(peptide,peptide_sequence %in% subset(peptide,concentration != "?")$peptide_sequence & peptide_intensity > 0),.(run_id,protein_id,peptide_sequence),function(X){if(dim(subset(X,concentration=="?"))[1] == dim(subset(X,concentration!="?"))[1]){unique(data.frame("run_id"=X$run_id,"protein_id"=X$protein_id,"peptide_sequence"=X$peptide_sequence,"concentration"=(subset(X,concentration=="?")$peptide_intensity/subset(X,concentration!="?")$peptide_intensity)*as.numeric(subset(X,concentration!="?")$concentration)),stringsAsFactors=FALSE)}})

		# third step: the endogenous peptides with know assigned concentrations are averaged to compute a mean protein concentration
		conc_proteins<-ddply(conc_peptides,.(run_id,protein_id),function(X){data.frame("run_id"=X$run_id,"protein_id"=X$protein_id,"concentration"=mean(X$concentration),stringsAsFactors=FALSE)})

		peptide<-merge(subset(peptide,concentration=="?", select=-c(concentration)), conc_proteins, by=c("run_id","protein_id"), all.x=TRUE)
		peptide$concentration[ is.na(peptide$concentration) ] <- "?"

	}

	# if the data is on the peptide level
	if ("peptide_intensity" %in% names(peptide)) {
		protein <- protein_inference.ProteinInference(peptide, peptide_method = peptide_method, peptide_topx = peptide_topx, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, fasta = fasta, model = model, combine_precursors = combine_precursors, ...)
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

	return(result)
}

protein_inference.ProteinInference <- function(data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, ...) {
	run_id <- protein_id <- peptide_id <- peptide_sequence <- concentration <- response <- peptide_intensity <- protein_sequence_length <- protein_intensity <- NULL

	if (!is.na(fasta)) {
		peptide_sequences.fasta <- trypsin(fasta)
	}

	# should precursors be summed?
	if (combine_precursors) {
		data <- ddply(data, .(run_id,protein_id,peptide_id,peptide_sequence,concentration), function(X) data.frame("peptide_intensity" = sum(X$peptide_intensity)))
	}
			
	if (peptide_method == "top") {
		data.preselection <- ddply(data, .(run_id,protein_id), function(X) top(X, target = "peptide_intensity", topx = peptide_topx, strictness = peptide_strictness, decreasing = TRUE))[c("run_id","protein_id","peptide_id","concentration","peptide_intensity")]
				
		if (peptide_summary == "mean") {
			data.selection <- ddply(data.preselection, .(run_id,protein_id,concentration), function(X) data.frame("response" = mean(X$peptide_intensity)))
		}
		else if (peptide_summary == "median") {
			data.selection <- ddply(data.preselection, .(run_id,protein_id,concentration), function(X) data.frame("response" = median(X$peptide_intensity)))
		}
		else if (peptide_summary == "sum") {
			data.selection <- ddply(data.preselection, .(run_id,protein_id,concentration), function(X) data.frame("response" = sum(X$peptide_intensity)))
		}
	}
	else if (peptide_method == "atop") {
		peptide_sequences.af <- apexFeatures(data.frame("peptide_sequence" = data$peptide_sequence, "apex"=NA, stringsAsFactors=FALSE))
        peptide_sequences.apex <- predict(model,peptide_sequences.af)
		peptide_sequences.response <- peptide_sequences.apex$prediction[,c("peptide_sequence","apex")]
		names(peptide_sequences.response) <- c("peptide_sequence","response")
		
		data <- merge(data,peptide_sequences.response, all.x = TRUE)

		data.preselection <- ddply(data, .(run_id,protein_id,concentration), function(X) top(X, target = "response", topx = peptide_topx, strictness = peptide_strictness))[c("run_id","protein_id","concentration","peptide_intensity")]
		
		if (peptide_summary == "mean") {
			data.selection <- ddply(data.preselection, .(run_id,protein_id,concentration), function(X) data.frame(mean(X$peptide_intensity)))
		}
		else if (peptide_summary == "median") {
			data.selection <- ddply(data.preselection, .(run_id,protein_id,concentration), function(X) data.frame(median(X$peptide_intensity)))
		}
		else if (peptide_summary == "sum") {
			data.selection <- ddply(data.preselection, .(run_id,protein_id,concentration), function(X) data.frame(sum(X$peptide_intensity)))
		}
		
		names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	else if (peptide_method == "all") {
	
		if (peptide_summary == "mean") {
			data.selection <- ddply(data, .(run_id,protein_id,concentration), function(X) data.frame(mean(X$peptide_intensity)))
		}
		else if (peptide_summary == "median") {
			data.selection <- ddply(data, .(run_id,protein_id,concentration), function(X) data.frame(median(X$peptide_intensity)))
		}
		else if (peptide_summary == "sum") {
			data.selection <- ddply(data, .(run_id,protein_id,concentration), function(X) data.frame(sum(X$peptide_intensity)))
		}
		
		names(data.selection) <- c("run_id","protein_id","concentration","response")	
	}
	else if (peptide_method == "iBAQ") {
		data.selection <- ddply(data, .(run_id,protein_id,concentration), function(X) data.frame(sum(X$peptide_intensity)/length(unlist((lapply(as.list(peptide_sequences.fasta[[X$protein_id[1]]]),function(X){if(nchar(X)>=6 && nchar(X)<=30){return(X)}else{return(NA)}})))[!is.na(unlist((lapply(as.list(peptide_sequences.fasta[[X$protein_id[1]]]),function(X){if(nchar(X)>=6 && nchar(X)<=30){return(X)}else{return(NA)}}))))])))
		
		names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	else if (peptide_method == "APEX") {
		peptide_sequences.af <- apexFeatures(data.frame("peptide_sequence" = data$peptide_sequence, "apex"=NA, stringsAsFactors=FALSE))
        peptide_sequences.apex <- predict(model,peptide_sequences.af)
    	peptide_sequences.response <- peptide_sequences.apex$prediction[,c("peptide_sequence","apex")]
        names(peptide_sequences.response) <- c("peptide_sequence","response")

        data <- merge(data,peptide_sequences.response, all.x = TRUE)
        
        data <- ddply(data,.(run_id,protein_id,peptide_id,peptide_sequence,concentration,response),function(X) {sum(X$peptide_intensity)})
        names(data) <- c("run_id","protein_id","peptide_id","peptide_sequence","concentration","response","peptide_intensity")
                			
		data.selection <- ddply(data, .(run_id,protein_id,concentration), function(X) sum(X$peptide_intensity)/sum(data$response[which(data$peptide_sequence %in% unlist(unname(peptide_sequences.fasta[X$protein_id[1]])))]))

		names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	else if (peptide_method == "TSC") {
		data.dt<-data.table(data)
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-data.dt[, list("peptide_intensity"=sum(peptide_intensity)), by=key(data.dt)]

		data.selection<-data.frame(data.dt)

        names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	else if (peptide_method == "NSAF") {
		proteins.fasta <- read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)

		proteins <- ldply(proteins.fasta, function(x) ldply(x))
		names(proteins) <- c("protein_id","protein_sequence")
		proteins$protein_sequence_length<-nchar(proteins$protein_sequence)
		proteinlength<-proteins[,c("protein_id","protein_sequence_length")]

		data<-merge(data,proteinlength,by="protein_id")

		data.dt<-data.table(data)
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data2.dt<-data.dt[, list("peptide_intensity"=sum(peptide_intensity)/protein_sequence_length), by=key(data.dt)]

		setkey(data2.dt,"run_id")
		data2.dt<-data2.dt[, list("protein_id"=protein_id,"concentration"=concentration,"peptide_intensity"=peptide_intensity/sum(peptide_intensity)), by=key(data2.dt)]

		data.selection<-data.frame(data2.dt)

		nasf_factor<-10^(round(log10(length(unique(data.selection$protein_id)))))

		data.selection$peptide_intensity<-data.selection$peptide_intensity*nasf_factor

        names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	else if (peptide_method == "MiST") {
		proteins.fasta <- read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)

		proteins <- ldply(proteins.fasta, function(x) ldply(x))
		names(proteins) <- c("protein_id","protein_sequence")
		proteins$protein_sequence_length<-nchar(proteins$protein_sequence)
		proteinlength<-proteins[,c("protein_id","protein_sequence_length")]

		data<-merge(data,proteinlength,by="protein_id")

		data.dt<-data.table(data)
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		# summary of peptides to proteins
		data2.dt<-data.dt[, list("protein_intensity"=sum(peptide_intensity)/protein_sequence_length), by=key(data.dt)]

		setkey(data2.dt,"run_id")
		data2.dt<-data2.dt[, list("protein_id"=protein_id,"concentration"=concentration,"protein_intensity"=protein_intensity/sum(protein_intensity)), by=key(data2.dt)]

		data.selection<-data.frame(data2.dt)

        names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	else if (peptide_method == "alm") {
		peptide_sequences.af <- apexFeatures(data.frame("peptide_sequence" = data$peptide_sequence, "apex"=NA, stringsAsFactors=FALSE))
		peptide_sequences.apex <- predict(model,peptide_sequences.af)
		peptide_sequences.response <- peptide_sequences.apex$prediction[,c("peptide_sequence","apex")]
		names(peptide_sequences.response) <- c("peptide_sequence","response")
		
		data <- merge(data,peptide_sequences.response, all.x = TRUE)
		
		data$response <- log(data$response)
			
		data.selection <- ddply(data[with(data, order(run_id,protein_id, response)), ], .(run_id,protein_id,concentration), function(X) data.frame(lm(peptide_intensity ~ response,X)$coefficients[[1]]))
		names(data.selection) <- c("run_id","protein_id","concentration","response")
	}
	
	return(data.selection)
}

peptide_inference.ProteinInference <- function(data, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", ...) {
	run_id <- protein_id <- peptide_id <- peptide_sequence <- precursor_charge <- concentration <- NULL

	if (transition_summary == "mean") {
		return(ddply(data, .(run_id,protein_id,peptide_id,peptide_sequence,precursor_charge,concentration), function(X) {data.frame("peptide_intensity"=mean(top(X,target = "transition_intensity", topx = transition_topx, strictness = transition_strictness, decreasing = TRUE)$transition_intensity))}))
	} else if (transition_summary == "median") {
		return(ddply(data, .(run_id,protein_id,peptide_id,peptide_sequence,precursor_charge,concentration), function(X) {data.frame("peptide_intensity"=median(top(X,target = "transition_intensity", topx = transition_topx, strictness = transition_strictness, decreasing = TRUE)$transition_intensity))}))
	} else if (transition_summary == "sum") {
		return(ddply(data, .(run_id,protein_id,peptide_id,peptide_sequence,precursor_charge,concentration), function(X) {data.frame("peptide_intensity"=sum(top(X,target = "transition_intensity", topx = transition_topx, strictness = transition_strictness, decreasing = TRUE)$transition_intensity))}))
	}
}

trypsin <- function(fasta, ...)  UseMethod("trypsin")

trypsin.default <- function(fasta, ...) {
	proteins <- read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)

	sequences <- sapply(proteins, strsplit, "(?!P)(?<=[RK])", perl = TRUE)
	return(sequences)
}

top <- function(data, ...)  UseMethod("top")

top.default <- function(data, target = names(data)[dim(data)[1]], topx = 3, strictness = "loose", decreasing = TRUE, ...) {
	if (strictness == "strict") {
		if (dim(data)[1] >= topx) {
			return(data[with(data, order(data[,target], decreasing = decreasing)), ][1:topx,])
		}
	}
	else if (strictness == "loose")  {
		return(na.omit(data[with(data, order(data[,target], decreasing = decreasing)), ][1:topx,]))
	}
}
