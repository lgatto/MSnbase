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
  ref <- .referenceFraction(x, group)
  nr <- nrow(x)
  x[(ref[group] - 1L) * nr + 1L:nr]
}

#' norm to reference run/sample for each protein
#'
#' described in: https://doi.org/10.1104/pp.114.245589
#' The peptide data were converted to protein intensities as follows. For each
#' protein, a fraction with the highest number of peptides quantified was
#' nominated as a reference fraction. The protein abundance in that reference
#' fraction was taken as one. Then, other fractions were quantified against the
#' reference fraction by finding the ratio of the combined intensity for
#' peptides shared between the interrogated and reference fractions. When the
#' quantities in all fractions for a given protein were computed, the values
#' were renormalized to give a sum = 1 across all 10 fractions.
#'
#' @param x exprs matrix
#' @param group grouping variable, i.e. protein accession
#' @param reference double, vector of reference values (has to be of length
#'  group && nrow(x))
#' @param norm normalise proteins to sum 1?
#' @return column index of reference fraction per group
#' @noRd
.normToReference <- function(x, group,
                             reference=.referenceFractionValues(x=x, group=group),
                             norm=TRUE) {
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(reference))
  stopifnot(nrow(x) == length(group))
  stopifnot(length(reference) == length(group))

  mRef <- matrix(reference, nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
  x[is.na(mRef)] <- NA_real_ # mark all values as NA where reference is NA
  mRef[is.na(x)] <- NA_real_
  r <- rowsum(x, group=group, reorder=FALSE, na.rm=TRUE) /
        rowsum(mRef, group=group, reorder=FALSE, na.rm=TRUE)
  if (norm) {
    r <- r/rowSums(r, na.rm=TRUE)
  }
  r
}


