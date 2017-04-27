#' find reference run/sample of each group
#'
#' @param x exprs matrix
#' @param group grouping variable, i.e. protein accession
#' @return column index of reference fraction per group
#' @noRd
.referenceFraction <- function(x, group) {
  notNA <- !is.na(x)
  mode(notNA) <- "numeric"
  rs <- rowsum(x, group=group, reorder=FALSE, na.rm=TRUE)
  rsNA <- rowsum(notNA, group=group, reorder=FALSE, na.rm=TRUE)
  setNames(max.col((.rowMaxs(rsNA) == rsNA) * rs), rownames(rs))
}

#' find reference run/sample value of each group
#'
#' @param x exprs matrix
#' @param group grouping variable, i.e. protein accession
#' @return double, vector of length group && nrow(x) with reference values
#' @noRd
.referenceFractionValues <- function(x, group) {
  stopifnot(is.matrix(x))
  stopifnot(nrow(x) == length(group))
  ## if group is a factor the order would be different
  group <- as.character(group)
  ref <- .referenceFraction(x, group)
  nr <- nrow(x)
  x[(ref[group] - 1L) * nr + 1L:nr]
}

#' Combine peptides into proteins.
#'
#' This function combines peptides into their proteins by normalising the
#' intensity values to a reference run/sample for each protein.
#'
#' This function is not intented to be used directly (that's why it is not
#' exported via \code{NAMESPACE}). Instead the user should use
#' \code{\link[MSnbase]{combineFeatures}}.
#'
#' The algorithm is described in Nikolovski et al., briefly it works as
#' follows:
#'
#' \enumerate{
#'  \item Find reference run (column) for each protein (grouped rows).
#'  We use the run (column) with the lowest number of \code{NA}.
#'  If multiple candidates are available we use the one with the highest
#'  intensity. This step is skipped if the user use his own \code{reference}
#'  vector.
#'  \item For each protein (grouped rows) and each run (column):
#'  \enumerate{
#'    \item Find peptides (grouped rows) shared by the current run (column) and
#'    the reference run (column).
#'    \item Sum the shared peptides (grouped rows) for the current run (column)
#'    and the reference run (column).
#'    \item The ratio of the shared peptides (grouped rows) of the current run
#'    (column) and the reference run (column) is the new intensity for the
#'    current protein for the current run.
#'  }
#' }
#'
#' @param x \code{matrix}, \code{\link{exprs}} matrix of an
#'   \linkS4class{MSnSet} object.
#' @param group \code{double} or \code{factor}, grouping variable,
#'   i.e. protein accession; has to be of length equal \code{nrow(x)}.
#' @param reference \code{double}, vector of reference values, has to be of the
#' same length as \code{group} and \code{nrow(x)}.
#' @author Sebastian Gibb <mail@@sebastiangibb.de>, Pavel Shliaha
#' @return a matrix with one row per protein.
#' @references
#' Nikolovski N, Shliaha PV, Gatto L, Dupree P, Lilley KS. Label-free protein
#' quantification for plant Golgi protein localization and abundance. Plant Physiol.
#' 2014 Oct;166(2):1033-43. DOI: 10.1104/pp.114.245589. PubMed PMID: 25122472.
#' @seealso \code{\link[MSnbase]{combineFeatures}}
#' @aliases NTR
#' @examples
#' library("MSnbase")
#' data(msnset)
#'
#' # choose the reference run automatically
#' combineFeatures(msnset, groupBy=fData(msnset)$ProteinAccession)
#'
#' # use a user-given reference
#' combineFeatures(msnset, groupBy=fData(msnset)$ProteinAccession,
#'  reference=rep(2, 55))
#'
normToReference <- function(x, group,
                            reference=.referenceFractionValues(x=x,
                                                               group=group)) {
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(reference))
  stopifnot(nrow(x) == length(group))
  stopifnot(length(reference) == length(group))

  mRef <- matrix(reference, nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
  x[is.na(mRef)] <- NA_real_ # mark all values as NA where reference is NA
  mRef[is.na(x)] <- NA_real_
  rowsum(x, group=group, reorder=FALSE, na.rm=TRUE) /
    rowsum(mRef, group=group, reorder=FALSE, na.rm=TRUE)
}
