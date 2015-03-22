imputeMethods <- function()
    c("bpca","knn", "QRILC", "MLE",
      "MinDet", "MinProb", "min", "zero",
      "mixed", "nbavg")


setMethod("impute", "MSnSet",
          function(object,
                   method,
                   randna,
                   mar,
                   mnar,
                   ...) {
              if (missing(method))
                  stop("Please specify an imputation method. ",
                       "See '?impute' for details.")              
              method <- match.arg(method,
                                  choices = imputeMethods(),
                                  several.ok = FALSE)
              
              if (method %in% c("CRILC", "MinDet", "MinProb"))
                  if (!require("imputeLCMD"))
                      stop("Method ", method,
                           "requires the imputeLCMD package.")              
              ## 
              ## imputaton methods
              ##
              if (method == "knn") {
                  .eset <- impute.knn(exprs(object), ...)$data
                  exprs(object) <- .eset
              } else if (method == "nbavg") {
                  message("Assuming values are ordered.")
                  impargs <- pairlist(...)                  
                  if (is.null(impargs$k)) k <- min(exprs(object), na.rm = TRUE)
                  else k <- impargs$k
                  exprs(object) <- imp_neighbour_avg(exprs(object),
                                                     k = k)
              } else if (method == "MLE") {
                  require("norm") || stop("Package 'norm' is required.")
                  x <- exprs(object)
                  s <- prelim.norm(x)  ## preliminary manipulations
                  th <- em.norm(s, ...) ## find the MLE
                  seed <- sample(.Machine$integer.max, 1)
                  rngseed(seed) ## set random number generator seed
                  exprs(object) <-
                      imp.norm(s, th, x)  ## impute missing data under the MLE
              } else if (method == "bpca"){
                  nSamples <- dim(exprs(object))[2]
                  .resultBPCA <- pca(exprs(object), method = "bpca",
                                     nPcs = (nSamples-1), verbose = FALSE, ...)
                  exprs(object) <- completeObs(.resultBPCA)
              } else if (method == "QRILC") {
                  exprs(object) <-
                      imputeLCMD::impute.QRILC(exprs(object), ...)[[1]]
              } else if (method == "MinDet") {
                  exprs(object) <-
                      imputeLCMD::impute.MinDet(exprs(object), ...)
              } else if (method == "MinProp") {
                  exprs(object) <-
                      imputeLCMD::impute.MinProb(exprs(object), ...)
              } else if (method == "min") {
                  val <- min(exprs(object), na.rm = TRUE)
                  exprs(object)[is.na(exprs(object))] <- val
              } else if (method == "mixed") {
                  if (missing(randna))
                      stop("Mixed imputation requires 'randna' argument. See ?impute.",
                           call. = FALSE)
                  stopifnot(is.logical(randna))
                  if (missing(mar))
                      stop("Mixed imputation requires 'mar' argument. See ?impute.",
                           call. = FALSE)
                  if (missing(mnar))
                      stop("Mixed imputation requires 'mnar' argument. See ?impute.",
                           call. = FALSE)
                  if (length(randna) != nrow(object))
                      stop("Number of proteins and length of randna must be equal.",
                           call. = FALSE)
                  exprs(object)[randna, ] <- exprs(impute(object[randna, ],
                                                          mar, ...))
                  exprs(object)[!randna, ] <- exprs(impute(object[!randna, ],
                                                           mnar, ...))
                  method <- paste(method, collapse = "/") ## for logging                  
              } else { ## method == "zero"
                  exprs(object)[is.na(exprs(object))] <- 0
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
                  if (!method %in% c("min", "zero"))
                      impargs <- "  Using default parameters"
              }
              object@processingData@processing <-
                  c(object@processingData@processing,
                    impargs)            
              if (validObject(object))
                  return(object)           
          })





