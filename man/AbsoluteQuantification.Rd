\name{AbsoluteQuantification}
\alias{AbsoluteQuantification}
\alias{AbsoluteQuantification.default}
\alias{predict.AbsoluteQuantification}
\alias{cval}
\alias{cval.default}
\alias{cval.AbsoluteQuantification}
\alias{print.AbsoluteQuantification}
\alias{plot.AbsoluteQuantification}
\alias{hist.AbsoluteQuantification}
\alias{pivot}
\alias{pivot.default}
\alias{pivot.AbsoluteQuantification}
\alias{export}
\alias{export.default}
\alias{export.AbsoluteQuantification}

\title{Absolute label-free quantification of mass spectrometry proteomics experiments}
\description{Absolute label-free quantification of mass spectrometry proteomics experiments.}
\usage{
\method{AbsoluteQuantification}{default}(data, peptide_method = "top", peptide_topx = 1,
  peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3,
  transition_strictness = "loose", transition_summary = "sum", fasta = NA,
  model = NA, total_protein_concentration = 1, combine_precursors = FALSE,
  consensus_peptides = TRUE, consensus_transitions = TRUE, ...)
\method{cval}{AbsoluteQuantification}(object, cval_method = "mc", mcx = 1000, ...)
\method{print}{AbsoluteQuantification}(x, ...)
\method{plot}{AbsoluteQuantification}(x, ...)
\method{hist}{AbsoluteQuantification}(x, ...)
\method{pivot}{AbsoluteQuantification}(x, ...)
\method{export}{AbsoluteQuantification}(x, file, ...)
}
\arguments{
  \item{data}{a mandatory data frame containing the columns \code{"run_id"}, 
    \code{"protein_id"}, \code{"protein_intensity"}, and \code{"concentration"} 
    for quantification on the protein level. For quantification on the peptide
    level, the columns \code{"run_id"}, \code{"protein_id"}, \code{"peptide_id"},
    \code{"peptide_sequence"}, \code{"precursor_charge"}, \code{"peptide_intensity"},
    and \code{"concentration"} are required. For quantification on the transition level,
    the columns \code{"protein_id"}, \code{"peptide_id"}, \code{"transition_id"},
    \code{"peptide_sequence"}, \code{"precursor_charge"}, \code{"transition_intensity"}
    and \code{"concentration"} are required. The id columns can be defined in any format,
    while the \code{"_intensity"} and \code{"concentration"} columns need to be numeric
    and in non-log form. The data may contain calibration data (with numeric
    \code{"concentration"} and test data (with \code{"concentration"} = "?"))}
  \item{peptide_method}{one of \code{"top"}, \code{"all"}, \code{"iBAQ"}, \code{"APEX"} or \code{"NSAF"} peptide to protein intensity estimation methods.}
  \item{peptide_topx}{a positive integer value of the top x peptides to consider for \code{"top"} methods.}
  \item{peptide_strictness}{whether \code{peptide_topx} should only consider proteins with the minimal peptide
    number (\code{"strict"}) or all (\code{"loose"}).}
  \item{peptide_summary}{how to summarize the peptide intensities for the \code{"top"} methods: \code{"mean"}, 
    \code{"median"}, \code{"sum"}.}
  \item{transition_topx}{a positive integer value of the top x transitions to consider for transition to 
    peptide intensity estimation methods.}
  \item{transition_strictness}{whether \code{transition_topx} should only consider peptides with the minimal transition
    number (\code{"strict"}) or all (\code{"loose"}).}
  \item{transition_summary}{how to summarize the transition intensities: \code{"mean"}, 
    \code{"median"}, \code{"sum"}.}
  \item{fasta}{the path and filename to an amino acid fasta file containing the proteins of interest, which
    is used for the \code{"iBAQ"}, \code{"APEX"} and \code{"NSAF"} peptide to protein intensity estimation methods.}
  \item{model}{the name of the \code{"APEX"} model to use}
  \item{total_protein_concentration}{the total protein concentration in the sample in any unit. This will be used for the normalized protein and concentration columns.}
  \item{combine_precursors}{whether to sum all precursors of the same peptide.}
  \item{object}{an \code{AbsoluteQuantification} object.}
  \item{cval_method}{a method for doing crossvalidation: \code{"boot"} (bootstrapping), \code{"mc"} (monte carlo cross-validation), \code{"loo"} (leaving-one-out).}
  \item{consensus_peptides}{if multiple runs are provided, select identical peptides among all runs.}
  \item{consensus_transitions}{if multiple runs are provided, select identical transitions among all runs.}
  \item{mcx}{a positive integer value of the number of folds for cross-validation.}
  
  \item{file}{the location of the output csv file.}
  
  \item{x}{an \code{AbsoluteQuantification} object.}
  \item{...}{future extensions.}
}

\details{
If absolute quantity estimation based on anchor peptides or proteins is demanded, the calibration peptide or protein abundance must be provided. Both estimated calibration protein intensities and separately determined calibration protein concentrations are log transformed and a first order linear least-squares regression of this log-log data is calculated. The abundance of the target proteins is predicted based on this regression. The error of the regression arises from biological and technical variation as well from the protein and peptide intensity estimators. To perform model selection and to estimate the error of the predicted protein concentrations, bootstrapping and Monte Carlo cross-validation as suggested (Malmstrom et al., 2009; Ludwig et al., 2012) were implemented. For both methods, the objective function is the minimization of the mean fold-error. 

If, on the other hand, the total protein concentration per cell is supplied in proteome-wide experiments, the absolute protein concentrations are estimated by normalization of the MS intensities or spectral counts to this number (Lu et al., 2006).
}
\value{
  An object of class \code{AbsoluteQuantification}.
}
\references{
Malmstrom, J. et al. \emph{Proteome-wide cellular protein concentrations of the human pathogen Leptospira interrogans.} Nature 460, 762-765 (2009).

Ludwig, C., Claassen, M., Schmidt, A. & Aebersold, R. \emph{Estimation of Absolute Protein Quantities of Unlabeled Samples by Selected Reaction Monitoring Mass Spectrometry.} Molecular & Cellular Proteomics 11, M111.013987-M111.013987 (2012).

Lu, P., Vogel, C., Wang, R., Yao, X. & Marcotte, E. M. \emph{Absolute protein expression profiling estimates the relative contributions of transcriptional and translational regulation.} Nat Biotech 25, 117-124 (2006).
}
\author{George Rosenberger \email{rosenberger@imsb.biol.ethz.ch}}

\seealso{\code{\link{import}}, \code{\link{ProteinInference}}, \code{\link{ALF}}, \code{\link{APEX}}, \code{\link{apexFeatures}}, \code{\link{proteotypic}}}

\examples{
data(UPS2MS)

UPS2_SRM<-head(UPS2_SRM,100) # Remove this line for real applications

data_AbsoluteQuantification <- cval(AbsoluteQuantification(UPS2_SRM),mcx=2)
print(data_AbsoluteQuantification)
plot(data_AbsoluteQuantification)
hist(data_AbsoluteQuantification)
pivot(data_AbsoluteQuantification)
}

\keyword{ AbsoluteQuantification }
\keyword{ absolute }
\keyword{ label-free }
\keyword{ quantification }
