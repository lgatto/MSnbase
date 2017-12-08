#' SummarizedExperiment to MSnSet object conversion
#'
#' \code{se2msnset} generates an MSnSet object from an SummarizedExperiment object.
#'
#' @param se SummarizedExperiment,
#' Object which will be turned into an MSnSet object.
#' @return MSnSet object.
#' @export
se2msnset <- function(se) {
  # Check input
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Extract expression, feature and pheno data
  raw <- SummarizedExperiment::assay(se)
  feat_data <- data.frame(SummarizedExperiment::rowData(se))
  rownames(feat_data) <- se@NAMES
  pheno_data <- data.frame(SummarizedExperiment::colData(se))

  # Generate MSnSet
  msnset <- MSnbase::MSnSet(exprs = as.matrix(raw),
                            pData = Biobase::AnnotatedDataFrame(pheno_data),
                            fData = Biobase::AnnotatedDataFrame(feat_data))
  return(msnset)
}

#' MSnSet to SummarizedExperiment object conversion
#'
#' \code{msnset2se} generates an SummarizedExperiment object from an MSnSet object.
#'
#' @param msnset MSnSet object,
#' Object which will be turned into an SummarizedExperiment object.
#' @return SummarizedExperiment object.
#' @export
msnset2se <- function(msnset) {
  # Check input
  assertthat::assert_that(inherits(msnset, "MSnSet"))

  # Extract assay, rowData and colData
  raw <- Biobase::exprs(msnset)
  row_data <- Biobase::fData(msnset)
  col_data <- Biobase::pData(msnset)

  # Generate SE
  se <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(raw),
                                                   rowData = row_data,
                                                   colData = col_data)
  return(se)
}
