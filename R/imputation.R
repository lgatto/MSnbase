##' @title Quantitative proteomics data imputation
##'
##' @description
##'
##' The `impute` method performs data imputation on `MSnSet` instances
##' using a variety of methods.
##'
##' Users should proceed with care when imputing data and take
##' precautions to assure that the imputation produce valid results,
##' in particular with naive imputations such as replacing missing
##' values with 0.
##'
##' See [MsCoreUtils::impute_matrix()] for details on the different
##' imputation methods available and strategies.
##'
##' @param object An `MSnSet` object with missing values to be
##'     imputed.
##' 
##' @param method `character(1)` defining the imputation method. See
##'     `MsCoreUtils::imputeMethods()` for available ones. See
##'     [MsCoreUtils::impute_matrix()] for details.
##' 
##' @param ... Additional parameters passed to the inner imputation
##'     function. See [MsCoreUtils::impute_matrix()] for details.
##' 
##' @rdname impute
##'
##' @aliases impute,MSnSet-method naset
##'
##' @examples
##'
##' data(naset)
##' 
##' ## table of missing values along the rows
##' table(fData(naset)$nNA)
##'
##' ## table of missing values along the columns
##' pData(naset)$nNA
##'
##' ## non-random missing values
##' notna <- which(!fData(naset)$randna)
##' length(notna)
##' notna
##'
##' impute(naset, method = "min")
##' 
##' if (require("imputeLCMD")) {
##'     impute(naset, method = "QRILC")
##'     impute(naset, method = "MinDet")
##' }
##' 
##' if (require("norm"))
##'     impute(naset, method = "MLE")
##' 
##' impute(naset, "mixed",
##'        randna = fData(naset)$randna,
##'        mar = "knn", mnar = "QRILC")
##'
##'
##' ## neighbour averaging
##' x <- naset[1:4, 1:6]
##' 
##' exprs(x)[1, 1] <- NA ## min value
##' exprs(x)[2, 3] <- NA ## average
##' exprs(x)[3, 1:2] <- NA ## min value and average
##' ## 4th row: no imputation
##' exprs(x)
##'
##' exprs(impute(x, "nbavg"))
setMethod("impute", "MSnSet",
          function(object, method, ...) {
              res <- MsCoreUtils::impute_matrix(exprs(object), method, ...)
              exprs(object) <- res
              if (missing(method))
                  method <- "user-defined function"
              object@processingData@processing <-
                  c(object@processingData@processing,
                    paste("Data imputation using",
                          method, date()))              
              if (validObject(object))
                  return(object)
          })
