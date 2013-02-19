setMethod("impute", "MSnSet",
          function(object, method = "knn", ...) {
            method <- match.arg(method)
            ## so far, only knn imputation
            .eset <- impute.knn(exprs(object), ...)$data
            exprs(object) <- .eset
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

