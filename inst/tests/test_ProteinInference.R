library(testthat)
# aqua_data
aqua_data_single_peptide<-data.table("run_id"=c(1,1,1,1,2,2,2,2),"protein_id"=c("PROTEINA","PROTEINA","PROTEINB","PROTEINB","PROTEINA","PROTEINA","PROTEINB","PROTEINB"),"peptide_id"=c("PEPTIDEK","PEPTIDEK(UniMod:259)","PEPTIDEC","PEPTIDEC","PEPTIDEK","PEPTIDEK(UniMod:259)","PEPTIDEC","PEPTIDEC"),"peptide_sequence"=c("PEPTIDEK","PEPTIDEK","PEPTIDEC","PEPTIDEC","PEPTIDEK","PEPTIDEK","PEPTIDEC","PEPTIDEC"),"precursor_charge"=c(2,2,2,3,2,3,2,3),"peptide_intensity"=c(100,80,60,40,95,75,55,35),concentration=c("?",10,"?","?","?",20,"?","?"),stringsAsFactors=FALSE)

aqua_data_multi_peptide<-data.table("run_id"=c(1,1,1,1,2,2,2,2),"protein_id"=c("PROTEINA","PROTEINA","PROTEINA","PROTEINA","PROTEINA","PROTEINA","PROTEINA","PROTEINA"),"peptide_id"=c("PEPTIDEK","PEPTIDEK(UniMod:259)","PEPTIDEC","PEPTIDEC(UniMod:259)","PEPTIDEK","PEPTIDEK(UniMod:259)","PEPTIDEC","PEPTIDEC(UniMod:259)"),"peptide_sequence"=c("PEPTIDEK","PEPTIDEK","PEPTIDEC","PEPTIDEC","PEPTIDEK","PEPTIDEK","PEPTIDEC","PEPTIDEC"),"precursor_charge"=c(2,2,2,3,2,3,2,3),"peptide_intensity"=c(100,80,60,40,95,75,55,35),concentration=c("?",10,"?",10,"?",20,"?",20),stringsAsFactors=FALSE)

# ProteinInference.default
test_that("ProteinInference.default: Absolute abundance estimation of endogenous proteins with spiked-in SIS peptides: Single peptides / Concentration", {
	expect_that(ProteinInference.default(aqua_data_single_peptide, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE, consensus_transitions = TRUE)$concentration,equals(c("12.5","?","25.3333333333333","?")))
	expect_that(ProteinInference.default(aqua_data_single_peptide, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE, consensus_transitions = TRUE)$response,equals(c(100,60,95,55)))
	expect_that(ProteinInference.default(aqua_data_multi_peptide, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE, consensus_transitions = TRUE)$concentration,equals(c("13.75","28.3809523809524")))
	expect_that(ProteinInference.default(aqua_data_multi_peptide, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE, consensus_transitions = TRUE)$response,equals(c(80,75)))
})

# peptide_data
peptide_data<-data.table("run_id"=c(1,1,1,1,2,2,2,2),"protein_id"=c("PROTEINA","PROTEINA","PROTEINB","PROTEINB","PROTEINA","PROTEINA","PROTEINB","PROTEINB"),"peptide_id"=c("ANDPEPTIDEA","PEPTIDEK","PEPTIDEC","PEPTIDEC","PEPTIDEK","PEPTIDEK","PEPTIDEC","PEPTIDEC"),"peptide_sequence"=c("ANDPEPTIDEA","PEPTIDEK","PEPTIDEC","PEPTIDEC","PEPTIDEK","PEPTIDEK","PEPTIDEC","PEPTIDEC"),"precursor_charge"=c(2,2,2,3,2,3,2,3),"peptide_intensity"=c(100,80,60,40,95,75,55,35),concentration=c("?","?","?","?","?","?","?","?"),stringsAsFactors=FALSE)

# protein_inference.ProteinInference
test_that("protein_inference.ProteinInference: Consensus peptide selection", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 3, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(80,50,95,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(80,50,95,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(80,60,95,55)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 3, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(100,60,95,55)))
})

test_that("protein_inference.ProteinInference: Combine precursors", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(100,60,95,55)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = TRUE, consensus_peptides = FALSE)$response,equals(c(100,100,170,90)))
})

test_that("protein_inference.ProteinInference: Summarization of peptides", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(180,100,170,90)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "median", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
})

test_that("protein_inference.ProteinInference: topx selection criterion", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(100,60,95,55)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 3, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
})

test_that("protein_inference.ProteinInference: loose / strict selection criterion", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 3, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "top", peptide_topx = 3, peptide_strictness = "strict", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(0,0,0,0)))
})

test_that("protein_inference.ProteinInference: all peptide_method", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "all", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "all", peptide_summary = "median", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(90,50,85,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "all", peptide_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(180,100,170,90)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "all", peptide_summary = "mean", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(80,50,95,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "all", peptide_summary = "median", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(80,50,95,45)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "all", peptide_summary = "sum", fasta = NA, model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(80,100,95,90)))
})

test_that("protein_inference.ProteinInference: iBAQ peptide_method", {
	data(UPS2MS)

	expect_that(protein_inference.ProteinInference(data.table(UPS2_LFQ), peptide_method = "iBAQ", peptide_summary = "sum", fasta = system.file("extdata","UPS2.fasta",package="aLFQ"), model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(25106081.31,125856278.00,87587312.93,95134656.33,17871715.00,6529159.17,28423.72,208145399.64,14028046.12,171435.00,347060.15,3673181.00,431140.67,156229765.07,1045140.00,12268096.30,60432526.25,4489209.67,67995597.05,37791831.53,16983291.54)))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "iBAQ", peptide_summary = "sum", fasta = system.file("extdata","example.fasta",package="aLFQ"), model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(40,100,47.5,90)))
})

test_that("protein_inference.ProteinInference: APEX peptide_method", {
	set.seed(131)

	data(APEXMS)

	APEX_ORBI<-head(APEX_ORBI,50) # Remove this line for real applications
	APEX_ORBI.af <- apexFeatures(APEX_ORBI)
	APEX_ORBI.apex <- APEX(data=APEX_ORBI.af)

	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "APEX", peptide_summary = "sum", fasta = system.file("extdata","example.fasta",package="aLFQ"), model = APEX_ORBI.apex, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(136.0544,257.7320,161.5646,231.9588), tolerance = .001))
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "APEX", peptide_summary = "sum", fasta = system.file("extdata","example.fasta",package="aLFQ"), model = APEX_ORBI.apex, combine_precursors = FALSE, consensus_peptides = FALSE)$response,equals(c(306.1224,257.7320,289.1156,231.9588), tolerance = .001))
})

test_that("protein_inference.ProteinInference: NSAF peptide_method", {
	expect_that(protein_inference.ProteinInference(peptide_data, peptide_method = "NSAF", peptide_summary = "sum", fasta = system.file("extdata","example.fasta",package="aLFQ"), model = NA, combine_precursors = FALSE, consensus_peptides = TRUE)$response,equals(c(0.1441441,0.4279279,0.1818182,0.4090909), tolerance = .001))
})

# transition_data
transition_data<-data.table("run_id"=c(1,1,1,1,2,2,2,2),"protein_id"=c("prot1","prot1","prot1","prot1","prot1","prot1","prot1","prot1"),"peptide_id"=c("pep1","pep1","pep2","pep2","pep1","pep1","pep2","pep2"),"peptide_sequence"=c("PEPTIDEA","PEPTIDEA","PEPTIDEC","PEPTIDEC","PEPTIDEA","PEPTIDEA","PEPTIDEC","PEPTIDEC"),"precursor_charge"=c(2,2,3,3,2,2,3,3),"transition_id"=c(1,2,3,4,5,2,3,4),"transition_intensity"=c(100,80,60,40,95,75,55,35),concentration=c("?","?","?","?","?","?","?","?"),stringsAsFactors=FALSE)

# peptide_inference.ProteinInference
test_that("peptide_inference.ProteinInference: Consensus peptide selection", {
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = TRUE)$peptide_intensity,equals(c(80,100,75,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 2, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = TRUE)$peptide_intensity,equals(c(80,100,75,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 1, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = TRUE)$peptide_intensity,equals(c(80,60,75,55)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(180,100,170,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 2, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(180,100,170,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 1, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(100,60,95,55)))
})

test_that("peptide_inference.ProteinInference: Summarization of peptides", {
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(180,100,170,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "mean", consensus_transitions = FALSE)$peptide_intensity,equals(c(90,50,85,45)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "median", consensus_transitions = FALSE)$peptide_intensity,equals(c(90,50,85,45)))
})

test_that("peptide_inference.ProteinInference: topx selection criterion", {
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(180,100,170,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 1, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(100,60,95,55)))
})

test_that("peptide_inference.ProteinInference: loose / strict selection criterion", {
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(180,100,170,90)))
	expect_that(peptide_inference.ProteinInference(transition_data, transition_topx = 3, transition_strictness = "strict", transition_summary = "sum", consensus_transitions = FALSE)$peptide_intensity,equals(c(0,0,0,0)))
})

# strictnessfilter.ProteinInference
test_that("trictnessfilter.ProteinInference", {
	expect_that(strictnessfilter.ProteinInference(c(1,2,3,4,NA), strictness="loose"),equals(c(1,2,3,4)))
	expect_that(strictnessfilter.ProteinInference(c(1,2,3,4,NA), strictness="strict"),equals(c(0)))
})

# trypsin.default
test_that("trypsin.default", {
	expect_that(trypsin.default(fasta=system.file("extdata","example.fasta",package="aLFQ"))$PROTEINA,equals(c("PEPTIDEK","ANDPEPTIDEA")))
})
