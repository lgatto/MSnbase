#' SummarizedExperiment to MSnSet object conversion
#'
#' \code{se2msnset} generates an MSnSet object from an SummarizedExperiment object.
#'
#' @param se SummarizedExperiment,
#' Object which will be turned into an MSnSet object.
#' @return MSnSet object.
se2msnset <- function(se) {
  # Check input
  stopifnot(inherits(se, "SummarizedExperiment"))
	
	if(!requireNamespace("SummarizedExperiment")) {
		stop("The SummarizedExperiment package is required for SummarizedExperiment to MSnSet conversion.")
	}

  # Extract expression, feature and pheno data
  raw <- SummarizedExperiment::assay(se)
  featData <- data.frame(SummarizedExperiment::rowData(se), row.names = names(se))
  phenoData <- data.frame(SummarizedExperiment::colData(se))

  # Generate MSnSet
  MSnSet(exprs = as.matrix(raw),
  			 pData = AnnotatedDataFrame(phenoData),
  			 fData = AnnotatedDataFrame(featData))
}

#' MSnSet to SummarizedExperiment object conversion
#'
#' \code{msnset2se} generates an SummarizedExperiment object from an MSnSet object.
#'
#' @param msnset MSnSet object,
#' Object which will be turned into an SummarizedExperiment object.
#' @return SummarizedExperiment object.
msnset2se <- function(msnset) {
  # Check input
  stopifnot(inherits(msnset, "MSnSet"))

	if(!requireNamespace("SummarizedExperiment")) {
		stop("The SummarizedExperiment package is required for MSnSet to SummarizedExperiment conversion.")
	}
	
  # Extract assay, rowData and colData
  raw <- exprs(msnset)
  rowData <- fData(msnset)
  colData <- pData(msnset)

  # Generate SE
  SummarizedExperiment::SummarizedExperiment(assays = as.matrix(raw),
  																					 rowData = rowData,
  																					 colData = colData)
}

setAs("MSnSet", 
			"SummarizedExperiment",
			function(from) {
				msnset2se(from)
			})

setAs("SummarizedExperiment", 
			"MSnSet",
			function(from) {
				se2msnset(from)
			})