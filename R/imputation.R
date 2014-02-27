##' Imputation of missing values
##'
##' @title Missing values imputation
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param method The method of imputation to use. Default value is bpca.
##' @references Oba et al., A Bayesian missing value estimation method for 
##' gene expression profile data, Bioinformatics (2003) 19 (16): 2088-2096.
##' @return An instance of class \code{"\linkS4class{MSnSet}"}.
##' @author Laurent Gatto, Samuel Wieczorek, Vasile-Cosmin Lazar
setMethod("impute", "MSnSet",
          function(object, method = c("bpca","knn"), ...) {
            method <- match.arg(method)
            ## so far, only knn imputation
            if (method == "knn"){
              .eset <- impute.knn(exprs(object), ...)$data
              exprs(object) <- .eset
            }
            else if (method=="bpca"){
             nSamples <- dim(exprs(object))[2]
             .resultBPCA = pca(exprs(object), method="bpca", nPcs=(nSamples-1), verbose=F)
             exprs(object) = completeObs(.resultBPCA)
            }
            impargs <- pairlist(...)
            object@processingData@processing <-
              c(object@processingData@processing,
                paste("Data imputation using",
                       method, date()))
            if (!is.null(impargs)) {
              impargs <- unlist(impargs)
              impargs <- paste(names(impargs), impargs, sep = "=")
              impargs <- paste0("  Using parameter(s) ", impargs)
            } else {
              impargs <- "  Using default parameters"
            }
            object@processingData@processing <-
              c(object@processingData@processing,
                impargs)            
            if (validObject(object))
              return(object)           
          })



