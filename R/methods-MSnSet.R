#######################################################################
## NOTE: MSnSet is based in the ExpressionSet class defined in the
##       Biobase package and should have the exact same functionnality.
##       Some of the code is heavily inspired and sometimes
##       directly copies from the respective ExpressionSet method.

setMethod("initialize", "MSnSet",
          function(.Object,
                   assayData,
                   phenoData,
                   featureData,
                   experimentData,
                   exprs = new("matrix"),
                   ... ) {
            if (missing(assayData)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(exprs, byrow = FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(exprs, byrow = TRUE)
              if (missing(experimentData))
                experimentData <- new("MIAPE")
              .Object <- callNextMethod(.Object,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        exprs = exprs,
                                        experimentData = experimentData,
                                        ...)
            } else if (missing(exprs)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(assayData, byrow = FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(assayData, byrow = TRUE)
              if (missing(experimentData))
                experimentData <- new("MIAPE")
              .Object <- callNextMethod(.Object,
                                        assayData = assayData,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        experimentData = experimentData,
                                        ...)
            } else stop("provide at most one of 'assayData' or 'exprs' to initialize MSnSet",
                     call. = FALSE)
            .Object@processingData <- new("MSnProcess")
            if (validObject(.Object))
              Biobase:::.harmonizeDimnames(.Object)
          })



setValidity("MSnSet", function(object) {
  msg <- validMsg(NULL, Biobase:::isValidVersion(object, "MSnSet"))
  msg <- validMsg(msg, Biobase::assayDataValidMembers(assayData(object),
                                                      c("exprs")))
  if ( nrow(qual(object)) != 0 ) {
    nrow.obs <- nrow(qual(object))
    nrow.exp <- nrow(object) * length(object$reporters)
    if (nrow.obs != nrow.exp)
      msg <- validMsg(msg,
                      "number of rows in assayData and qual slots do not match.")
  }
  if (!inherits(experimentData(object), "MIAPE"))
    msg <- validMsg(msg,
                    "experimentData slot in MSnSet must be 'MIAPE' object")
  if (!inherits(exprs(object), "matrix"))
    msg <- validMsg(msg,
                    "exprs(.) must be a matrix")
  if ("" %in% featureNames(object))
      msg <- validMsg(msg, "Empty string is not a valid feature name.")
  if (is.null(msg)) TRUE else msg
})

setMethod("exprs", signature(object = "MSnSet"),
          function(object) assayDataElement(object, "exprs"))

setReplaceMethod("exprs", signature(object = "MSnSet", value = "matrix"),
                 function(object, value)
                     assayDataElementReplace(object, "exprs", value))

setMethod("show","MSnSet",
          function(object) {
            callNextMethod()
            show(processingData(object))
            invisible(NULL)
          })


setMethod("normalize", "MSnSet",
          function(object,
                   method = c("sum", "max", "center.mean",
                              "center.median", "diff.median",
                              "quantiles", "quantiles.robust", "vsn"), ...)
              normalise_MSnSet(object, match.arg(method), ...)
          )

normalise <- normalize

setMethod("scale", "MSnSet",
          function(x, center = TRUE, scale = TRUE) {
            e <- scale(exprs(x), center = center, scale = scale)
            proc <- NULL
            if (!is.null(attr(e, "scaled:scale"))) {
              escale <- attr(e, "scaled:scale")
              attr(e, "scaled:scale") <- NULL
              escale <- paste(round(escale, ), 3)
              proc <- c(proc,
                        paste("Scaled",
                              paste0("[", escale, collapse = ", "), "]:",
                              date()))
            }
            if (!is.null(attr(e, "scaled:center"))) {
              ecenter <- attr(e, "scaled:center")
              attr(e, "scaled:center") <- NULL
              ecenter <- paste(round(ecenter, ), 3)
              proc <- c(proc,
                        paste("Centered",
                              paste0("[", ecenter, collapse = ", "), "]:",
                              date()))
            }
            exprs(x) <- e
            x@processingData@processing <-
              c(x@processingData@processing,
                proc)
            if (validObject(x))
              return(x)
          })

setMethod("purityCorrect",
          signature = signature("MSnSet", "matrix"),
          function(object, impurities) {
            if (ncol(impurities) != nrow(impurities))
                stop("Impurity matrix must be a square matrix")
            if (ncol(object) != ncol(impurities))
                stop("Impurity matrix should be ",
                     ncol(object), " by ", ncol(object))
            impurities <- t(impurities)
            .purcor <- function(x, .impurities = impurities) {
                keep <- !is.na(x)
                if (sum(keep) > 1)
                    x[keep] <- solve(.impurities[keep, keep], x[keep])
                x[x<0] <- NA
                return(x)
            }
            corr.exprs <- apply(exprs(object), 1, .purcor)
            exprs(object) <- t(corr.exprs)
            object@processingData@processing <-
              c(object@processingData@processing,
                paste0("Purity corrected: ", date()))
            if (validObject(object))
              return(object)
          })


setMethod("dim", "MSnSet",function(x) dim(exprs(x)))
setMethod("qual", "MSnSet", function(object) object@qual)

setMethod("fileNames",
          signature(object = "MSnSet"),
          function(object) processingData(object)@files)

## setReplaceMethod("fileNames",
##           signature(object="MSnSet", value="character"),
##           function(object, value) {
##             fileNames(object@processingData) <- value
##             return(object)
##           })

setMethod("processingData",
          signature(object="MSnSet"),
          function(object) object@processingData)

setMethod("msInfo","MSnSet",
          function(object) msInfo(experimentData(object)))

setMethod("expinfo","MSnSet",
          function(object) expinfo(experimentData(object)))

setMethod("exptitle","MSnSet",
          function(object) exptitle(experimentData(object)))

setMethod("expemail","MSnSet",
          function(object) expemail(experimentData(object)))

setMethod("ionSource","MSnSet",
          function(object) ionSource(experimentData(object)))

setMethod("analyser","MSnSet",
          function(object) analyser(experimentData(object)))

setMethod("analyzer","MSnSet",
          function(object) analyzer(experimentData(object)))

setMethod("detectorType","MSnSet",
          function(object) detectorType(experimentData(object)))

setMethod("meanSdPlot",
          signature="MSnSet",
          definition =
          function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
                   ylab = "sd", pch  = ".", plot = TRUE, ...)
          vsn::meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, plot=plot, ...))

t.MSnSet <- function(x) {
  ans <- new("MSnSet",
             exprs = t(exprs(x)),
             phenoData = featureData(x),
             featureData = phenoData(x),
             experimentData = experimentData(x),
             annotation = annotation(x))
  ans@processingData@processing <-
      x@processingData@processing
  ans <- logging(ans, "MSnSet transposed")
  if (validObject(ans))
      return(ans)
}


setMethod("[", "MSnSet", function(x, i, j, ...) {
  dim0 <- dim(x)
  .Object <- callNextMethod(...)
  dim1 <- dim(.Object)
  ## subsetting qual - requires pData(x)$mz!
  ## fn <- featureNames(.Object)
  ## reps <- match(.Object$mz,x$mz)
  ## qrows <- paste(rep(fn,each=length(reps)),reps,sep=".")
  ## .Object@qual <- .Object@qual[qrows,]
  .Object@qual <- data.frame()
  dim0 <- paste0("[", paste0(dim0, collapse = ","), "]")
  dim1 <- paste0("[", paste0(dim1, collapse = ","), "]")
  .Object@processingData@processing <-
      c(.Object@processingData@processing,
        paste0("Subset ", dim0, dim1, " ", date()))
  if (validObject(.Object))
      return(.Object)
})


setAs("MSnSet", "ExpressionSet",
      function (from)
      new("ExpressionSet",
          exprs = exprs(from),
          phenoData = phenoData(from),
          featureData = featureData(from),
          annotation = annotation(from),
          experimentData = as(experimentData(from), "MIAME"),
          protocolData = protocolData(from))
      )

as.ExpressionSet.MSnSet <- function(x) as(x,"ExpressionSet")

setAs("ExpressionSet", "MSnSet",
      function (from)
      new("MSnSet",
          exprs = exprs(from),
          phenoData = phenoData(from),
          featureData = featureData(from),
          annotation = annotation(from),
          protocolData = protocolData(from))
      )

as.MSnSet.ExpressionSet <- function(x) as(x, "MSnSet")

setAs("MSnSet", "data.frame",
      function (from) {
          ## MSnSet -> ExpressionSet -> data.frame
          from <- as(from, "ExpressionSet")
          as(from, "data.frame")
      })

as.data.frame.MSnSet <-
    function(x, row.names = NULL, optional = FALSE, ...) as(x, "data.frame")

ms2df <- function(x, fcols = fvarLabels(x)) {
    if (is.null(fcols)) {
        res <- data.frame(exprs(x))
    } else {
        sel <- fvarLabels(x) %in% fcols
        res <- data.frame(exprs(x),
                          fData(x)[, sel])
    }
    return(res)
}

setMethod("write.exprs",
          signature(x = "MSnSet"),
          function(x,
                   fDataCols = NULL,
                   fcol,
                   file = "tmp.txt", quote = FALSE,
                   sep = "\t", col.names = NA, ...) {
            res <- exprs(x)
            if (!missing(fcol))
                fDataCols <- fcol
            if (!is.null(fDataCols))
              res <- cbind(res, fData(x)[, fDataCols])
            write.table(res, file = file, quote = quote, sep = sep,
                        col.names = col.names, ...)
          })

setReplaceMethod("experimentData",
                 signature = signature(
                   object = "MSnSet",
                   value = "MIAPE"),
                 function(object, value) {
                   if (!validObject(value))
                     stop("Not a valid MIAPE instance.")
                   object@experimentData <- value
                   object
                 })



setMethod("combine",
          signature = signature(
            x = "MSnSet", y = "MSnSet"),
          function(x, y, ...) {
            if (class(x) != class(y))
              stop(paste("objects must be the same class, but are ",
                         class(x), ", ", class(y), sep=""))
            ## if (!isCurrent(x)[["MSnSet"]])
            ##     x <- updateObject(x)
            assayData(x) <- combine(assayData(x), assayData(y))
            phenoData(x) <- combine(phenoData(x), phenoData(y))
            featureData(x) <- combine(featureData(x), featureData(y))
            experimentData(x) <- combine(experimentData(x), experimentData(y))
            protocolData(x) <- combine(protocolData(x), protocolData(y))
            x@processingData <- combine(processingData(x), processingData(y))
            x@processingData@processing <- paste("Combined [",
                                                 paste(dim(x), collapse = ","),
                                                 "] and [",
                                                 paste(dim(y), collapse = ","),
                                                 "] MSnSets ", date(), sep = "")
            x@qual <- data.frame() ## dropping qual slot
            ## annotation -- constant / not used
            if (validObject(x))
              return(x)
          })


setMethod("topN", signature(object = "matrix"),
          function(object, groupBy, n=3, fun, ...) {
            if (missing(groupBy))
              stop("Specify how to group features to select top ", n, ".")
            if (missing(fun)) {
              fun <- sum
              if (ncol(object) > 1)
                message("Ranking features using their sum.")
            }
            rn <- rownames(object)
            idx <- by(object, groupBy, getTopIdx, n, fun, ...)
            object <- subsetBy(object, groupBy, idx)
            if (!is.null(rn)) {
              rownames(object) <- subsetBy(rn, groupBy, idx)
            } else {
              rownames(object) <- NULL
            }
            return(object)
          })

setMethod("topN", signature(object = "MSnSet"),
          function(object, groupBy, n=3, fun, ...) {
            if (missing(groupBy))
              stop("Specify how to group features to select top ", n, ".")
            if (missing(fun)) {
              fun <- sum
              if (ncol(object) > 1)
                message("Ranking features using their sum.")
            }
            idx <- by(exprs(object), groupBy, getTopIdx, n, fun, ...)
            fn <- subsetBy(featureNames(object), groupBy, idx)
            .eset <- subsetBy(exprs(object), groupBy, idx)
            if (!is.matrix(.eset))
              .eset <- matrix(.eset, ncol = 1)
            rownames(.eset) <- fn
            .proc <- processingData(object)
            .proc@processing <- c(.proc@processing,
                                  paste0("Selected top ", n,
                                         " features: ", date()))
            .fdata <- subsetBy(fData(object), groupBy, idx)
            ans <- new("MSnSet",
                       experimentData = experimentData(object),
                       exprs = .eset,
                       phenoData = phenoData(object),
                       featureData = new("AnnotatedDataFrame", data = .fdata),
                       annotation = object@annotation,
                       protocolData = protocolData(object))
            ans@processingData <- .proc
            featureNames(ans) <- fn
            if (validObject(ans))
              return(ans)
          })


getRatios <- function(x, log = FALSE) {
  ## x: a vector of numerics
  ## returns a vector of all xi/xj ratios
    x <- as.numeric(x)
    cmb <- combn(length(x), 2)
    r <- numeric(ncol(cmb))
    for (i in 1:ncol(cmb)) {
        j <- cmb[1, i]
        k <- cmb[2, i]
        ifelse(log,
               r[i] <- x[j] - x[k],
               r[i] <- x[j] / x[k])
    }
    return(r)
}


setMethod("exprsToRatios",
          "MSnSet",
          function(object, log = FALSE) {
            if (ncol(object) == 2) {
              ifelse(log,
                     r <- exprs(object)[, 1] - exprs(object)[, 2],
                     r <- exprs(object)[, 1] / exprs(object)[, 2])
              dim(r) <- c(length(r), 1)
            } else {
              r <- apply(exprs(object), 1, getRatios, log)
              r <- t(r)
            }
            rownames(r) <- featureNames(object)
            cmb <- combn(ncol(object), 2)
            ratio.description <-
                apply(cmb, 2,
                      function(x)
                          paste(sampleNames(object)[x[1]],
                                sampleNames(object)[x[2]],
                                sep = "/"))
            phenodata <- new("AnnotatedDataFrame",
                             data = data.frame(ratio.description))
            processingdata <- processingData(object)
            processingdata@processing <- c(processingdata@processing,
                                           paste("Intensities to ratios: ",
                                                 date(), sep = ""))
            message("Dropping protocolData.")
            res <- new("MSnSet",
                       exprs = r,
                       featureData = featureData(object),
                       phenoData = phenodata,
                       processingData = processingdata,
                       experimentData = experimentData(object))
            if (validObject(res))
                return(res)
          })


setMethod("exprsToRatios",
          "matrix",
          function(object, log = FALSE) {
            if (ncol(object) == 2) {
              ifelse(log,
                     r <- object[, 1] - object[, 2],
                     r <- object[, 1] / object[, 2])
              dim(r) <- c(length(r), 1)
            } else {
              r <- apply(object, 1, getRatios, log)
              r <- t(r)
              colnames(r) <-
                apply(combn(ncol(object), 2), 2,
                      paste, collapse = ".")
            }
            r
          })

setMethod("image", "MSnSet",
          function(x,
                   facetBy = NULL,
                   sOrderBy = NULL,
                   legend = "",
                   low, high,
                   fnames,
                   nmax = 50,
                   plot = TRUE) {
              ## get rid of 'no visible global function definition' note
              sample.name <- feature.id <- Expression <- NULL
              isFC <- any(exprs(x) < 0, na.rm = TRUE)
              xlong <- melt(exprs(x))
              colnames(xlong) <- c("feature.id", "sample.name", "Expression")
              xlong[['feature.id']] <- as.character(xlong[['feature.id']])
              xlong[['sample.name']] <- as.character(xlong[['sample.name']])
              xlong <- merge(xlong, fData(x), by.x = "feature.id", by.y = 0)
              xlong <- merge(xlong, pData(x), by.x = "sample.name", by.y = 0)
              x <- xlong
              if (!is.null(sOrderBy))
                  x[['sample.name']] <- reorder(x[['sample.name']], x[[sOrderBy]])
    
              if (!is.null(facetBy)) x$facetBy <- x[[facetBy]]
              p <- ggplot(x, aes(x = `sample.name`, y = `feature.id`,
                                 fill = `Expression`)) +
                  geom_raster() +
                  theme(
                      axis.text.x = element_text(angle = +90),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                  xlab("Sample names") +
                  ylab("Features")
              if (isFC) {
                  if (missing(low)) low <- "gold1"
                  if (missing(high)) high <- "#08306B"
                  p <- p + scale_fill_gradient2(legend, low = low,
                                                high = high,
                                                mid = "white")
              } else {
                  if (missing(low)) low <- "#F7FBFF"
                  if (missing(high)) high <- "#08306B"
                  p <- p + scale_fill_gradient(legend, low = low, high = high)
              }
              if (missing(fnames)) {
                  if (length(unique(x$`feature.id`)) > nmax) {
                      p <- p + theme(axis.text.y = element_blank())
                      p <- p + theme(axis.ticks.y = element_blank())
                  }
              } else {
                  if (!fnames) {
                      p <- p + theme(axis.text.y = element_blank())
                      p <- p + theme(axis.ticks.y = element_blank())
                  }
              }
              if (!is.null(facetBy))
                  p <- p + facet_grid( . ~ facetBy, scales='free', space='free')
              if (plot) plot(p)
              invisible(p)
          })


image2 <- function(x,
                   yticks = 10,
                   x.cex.axis = .75,
                   y.cex.axis = .75,
                   xlab = "Samples",
                   ylab = "Features",
                   ...) {
    if (inherits(x, "MSnSet"))
        x <- exprs(x)
    nc <- ncol(x)
    nr <- nrow(x)
    lab <- colnames(x)
    if (is.null(lab))
        lab <- 1:nc 
    graphics::image(t(x),
                    xlab = xlab, ylab = ylab,
                    xaxt = "n", yaxt = "n", ...)
    axis(1, seq(0, 1 , 1 / (nc - 1)),
         labels = lab,
         cex.axis = x.cex.axis)
    yticks <- seq(0, 1, 1 / (yticks - 1)) * nr
    axis(2, seq(0,1, 1 / (length(yticks) - 1)),
         labels = round(yticks, 0),
         cex.axis = y.cex.axis)
    invisible(NULL)
}

setMethod("plotNA", signature(object = "MSnSet"),
          function(object, pNA = .5) {
              if (pNA > 1)
                  pNA <- 1
              if (pNA < 0)
                  pNA <- 0
              X <- exprs(object)
              p <- plotNA_matrix(X, pNA)
              invisible(p)
          })


setMethod("plotNA", signature(object = "matrix"),
          function(object, pNA = .5) {
              if (pNA > 1)
                  pNA <- 1
              if (pNA < 0)
                  pNA <- 0
              p <- plotNA_matrix(object, pNA)
              invisible(p)
          })

setMethod("filterNA", signature(object = "matrix"),
          function(object, pNA = 0, pattern) {
              if (missing(pattern)) { ## using pNA
                  if (pNA > 1)
                      pNA <- 1
                  if (pNA < 0)
                      pNA <- 0
                  k <- apply(object, 1,
                             function(x) sum(is.na(x)) / length(x))
                  accept <- k <= pNA
                  if (sum(accept) == 1) {
                      ans <- matrix(object[accept, ], nrow = 1)
                      rownames(ans) <- rownames(object)[accept]
                  } else {
                      ans <- object[accept, ]
                  }
              } else { ## using pattern
                  accept <- getRowsFromPattern(object, pattern)
                  ans <- object[accept, ]
              }
              return(ans)
          })

setMethod("filterZero", "matrix",
          function(object, ...) {
              object[object == 0] <- NA
              object <- filterNA(object, ...)
              object[is.na(object)] <- 0
              object
          })


setMethod("filterNA", signature(object = "MSnSet"),
          function(object, pNA = 0, pattern, droplevels = TRUE) {
              if (missing(pattern)) { ## using pNA
                  if (pNA > 1)
                      pNA <- 1
                  if (pNA < 0)
                      pNA <- 0
                  k <- apply(exprs(object), 1,
                             function(x) sum(is.na(x)) / length(x))
                  accept <- k <= pNA
                  ans <- object[accept, ]
                  ans@processingData@processing <-
                      c(processingData(ans)@processing,
                        paste0("Removed features with more than ",
                               round(pNA, 3), " NAs: ", date()))
              } else { ## using pattern
                  accept <- getRowsFromPattern(exprs(object), pattern)
                  ans <- object[accept, ]
                  ans@processingData@processing <-
                      c(processingData(ans)@processing,
                        paste0("Removed features with according to pattern ",
                               pattern, " ", date()))
              }
              if (droplevels)
                  ans <- droplevels(ans)
              if (validObject(ans))
                  return(ans)
          })

setMethod("filterZero", signature = "MSnSet",
          function(object, ...) {
              exprs(object)[exprs(object) == 0] <- NA
              object <- filterNA(object, ...)
              ## updating processing log
              px <- object@processingData@processing
              lx <- length(px)
              ## if last one is 'Dropped ...', then update previous
              i <- ifelse(grep("Dropped featureData's levels", px[lx]),
                          lx - 1,
                          lx)
              px[i] <- sub("NAs", "zeros", px[i])
              exprs(object)[is.na(object)] <- 0
              object@processingData@processing <- px
              if (validObject(object))
                  return(object)
          })

is.na.MSnSet <- function(x) is.na(exprs(x))

droplevels.MSnSet <- function(x, ...) {
    fData(x) <- droplevels(fData(x), ...)
    x@processingData@processing <-
        c(x@processingData@processing,
          paste("Dropped featureData's levels", date()))
    if (validObject(x))
        return(x)
}

setMethod("log",
          signature = "MSnSet",
          function(x, base, ...) {
              exprs(x) <-  log(exprs(x), base)
              x@processingData@processing <-
                  c(x@processingData@processing,
                    paste0("Log transformed (base ", base, ") ", date()))
              if (validObject(x))
                  return(x)
          })


setMethod("MAplot",
          signature = "MSnSet",
          function(object, log.it = TRUE, base = 2, ...) {
              if (ncol(object) < 2)
                  stop("Need at least 2 samples to produce an MA plot.")
              if (log.it)
                  object <- log(object, base)
              if (ncol(object) == 2) {
                  x <- exprs(object[, 1])
                  y <- exprs(object[, 2])
                  sel <- !is.na(x) & !is.na(y)
                  x <- x[sel]
                  y <- y[sel]
                  M <- x - y
                  A <- (x + y)/2
                  ma.plot(A, M, ...)
              } else {
                  mva.pairs(exprs(object), log.it = FALSE, ...)
              }
          })

setMethod("addIdentificationData", c("MSnSet", "character"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol = c("spectrumFile", "acquisitionnum"),
                   verbose = isMSnbaseVerbose()) {
              addIdentificationData(object,
                                    id = mzID(id, verbose = verbose),
                                    fcol = fcol, icol = icol)
          })

setMethod("addIdentificationData", c("MSnSet", "mzIDClasses"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol = c("spectrumFile", "acquisitionnum"), ...) {
                       addIdentificationData(object, id = flatten(id),
                                             fcol = fcol, icol = icol)
                   })

setMethod("addIdentificationData", c("MSnSet", "data.frame"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol = c("spectrumFile", "acquisitionnum"), ...) {
                       fd <- fData(object)

                       if (!nrow(fd))
                           stop("No feature data found.")

                       fd$spectrum.file <- basename(fileNames(object)[fd$file])

                       fd <- utils.mergeSpectraAndIdentificationData(fd, id,
                                                                     fcol = fcol,
                                                                     icol = icol)

                       ## after adding the identification data we remove the
                       ## temporary data to avoid duplication and problems in quantify
                       cn <- colnames(fd)
                       keep <- cn[!(cn %in% c("spectrum.file"))]
                       fData(object)[, keep] <- fd[, keep, drop = FALSE]
                       if (validObject(object))
                           return(object)
          })

setMethod("removeNoId", "MSnSet",
          function(object, fcol = "pepseq", keep = NULL)
              utils.removeNoId(object, fcol, keep))

setMethod("removeMultipleAssignment", "MSnSet",
          function(object, fcol = "nprot")
              utils.removeMultipleAssignment(object, fcol))

setMethod("idSummary",
          signature = "MSnSet",
          function(object) {
              ## we temporaly add the spectrumFile information
              ## to our fData data.frame because utils.idSummary
              ## needs this information for matching
              fd <- fData(object)
              fd$spectrumFile <- basename(fileNames(object)[fd$file])
              return(utils.idSummary(fd))
          })

##############################################
## This should also be implemented for pSet!

## o MSnSet $ and [[ now dispatch on featureData,
##   instead of phenoData (as was inherited from
##   ExpressionSet via eSet)

## setMethod("$", "MSnSet", function(x, name) {
##   eval(substitute(featureData(x)$NAME_ARG, list(NAME_ARG=name)))
## })


## setReplaceMethod("$", "MSnSet", function(x, name, value) {
##   featureData(x)[[name]] = value
##   x
## })

## setMethod("[[", "MSnSet", function(x, i, j, ...) featureData(x)[[i]])

## setReplaceMethod("[[", "MSnSet",
##                  function(x, i, j, ..., value) {
##                      featureData(x)[[i, ...]] <- value
##                      x
##                  })


setAs("IBSpectra", "MSnSet",
      function (from, to = "MSnSet") {
          ans <- MSnSet(exprs = assayData(from)$ions,
                        fData = fData(from),
                        pData = pData(from))
          exp <- experimentData(from)
          ## the example data in isobar has MIAME
          ## experimental data ?!?!
          if (inherits(exp, "MIAPE"))
              ans@experimentData <- exp
          ans@protocolData <- protocolData(from)
          if (validObject(ans))
              return(ans)
      })

as.IBSpectra.MSnSet <- function(x) {
    ans <- MSnSet(exprs = assayData(x)$ions,
                  fData = fData(x),
                  pData = pData(x))
    exp <- experimentData(x)
    ## the example data in isobar has MIAME
    ## experimental data ?!?!
    if (inherits(exp, "MIAPE"))
        ans@experimentData <- exp
    ans@protocolData <- protocolData(x)
    if (validObject(ans))
        return(ans)
}


## setAs("MSnSet", "IBSpectra",
##       function (from, to = "IBSpectra") {
##           ## see IBSpectraTypes() for possible types
##           ## if (ncol(from)) == 2) ...
##           ## if (ncol(from)) == 2) ...      
##       })
