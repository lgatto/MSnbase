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
        stopifnot(all(fcols %in% fvarLabels(x)))
        res <- data.frame(exprs(x),
                          fData(x)[, fcols])
        colnames(res)[-seq_len(ncol(x))] <- fcols
    }
    return(res)
}

## suggested use: keep only non-empty slots - see .reduce_list in
## utils.R

setAs("AnnotatedDataFrame", "list",
      function (from) as.list(from@data))

setAs("MIAxE", "list",
      function (from) {
          nms <- slotNames(from)
          nms <- setdiff(nms, ".__classVersion__")
          ans <- vector("list", length = length(nms))
          names(ans) <- nms
          for (k in nms)
              ans[[k]] <- slot(from, k)
          ans
      })

setAs("MSnProcess", "list",
      function (from) {
          nms <- slotNames(from)
          nms <- setdiff(nms, ".__classVersion__")
          ans <- vector("list", length = length(nms))
          names(ans) <- nms
          for (k in nms)
              ans[[k]] <- slot(from, k)
          ans
      })

setAs("MSnSet",
    "SummarizedExperiment",
    function (from) {
        if (!requireNamespace("SummarizedExperiment"))
            stop("The SummarizedExperiment package is required",
                 "to coerce an MSnSet to a SummarizedExperiment.")
        from <- logging(from, "Coerced to SummarizedExperiment.")

        raw <- exprs(from)
        rowData <- fData(from)
        colData <- pData(from)
        proc <- processingData(from)
        metaData <- list(MSnbaseFiles = proc@files,
                         MSnbaseProcessing = proc@processing,
                         MSnbaseVersion = proc@MSnbaseVersion)

        SummarizedExperiment::SummarizedExperiment(
            assays = as.matrix(raw),
            rowData = rowData,
            colData = colData,
            metadata = metaData)
})

setAs("SummarizedExperiment",
    "MSnSet",
    function(from) {

        if (!requireNamespace("SummarizedExperiment"))
            stop("The SummarizedExperiment package is required",
                 "to coerce a SummarizedExperiment to an MSnSet.")

        raw <- SummarizedExperiment::assay(from)
        featData <- data.frame(
            SummarizedExperiment::rowData(from),
            row.names = names(from))
        phenoData <- data.frame(SummarizedExperiment::colData(from))

        ## Extract metadata based solely on MSnSet slot names
        .processingData <- metadata(from)$processingData
        .experimentData <- metadata(from)$experimentData
        .protocolData <- metadata(from)$protocolData
        .qual <- metadata(from)$qual

        msnset <- MSnSet(exprs = as.matrix(raw),
                         pData = AnnotatedDataFrame(phenoData),
                         fData = AnnotatedDataFrame(featData))

        ## Assign metadata based on class (which fails if the proper
        ## name above was missing and the variable was assign NULL)
        if (inherits(.experimentData, "MIAPE"))
            msnset@experimentData <- .experimentData
        if (inherits(.processingData, "MSnProcess"))
            msnset@processingData <- .processingData
        if (inherits(.protocolData, "AnnotatedDataFrame"))
            msnset@protocolData <- .protocolData
        if (inherits(.qual, "data.frame"))
            msnset@qual <- .qual
        return(msnset)
})

##' A helper function to add all metadata from an MSnSet to a
##' SummarizedExperiment object.
##'
##' @title Add all of an MSnSet's metadata
##' @param x An instance of class `SummarizedExperiment` created with
#          `x <- as(y, "MSnSet")`.
##' @param y An instance of class `MSnSet`
##' @return The `x` object with ``
##' @author Laurent Gatto
##' @noRd
addMSnSetMetadata <- function(x, y) {
    stopifnot(inherits(x, "SummarizedExperiment"))
    stopifnot(inherits(y, "MSnSet"))
    x@metadata$processingData <- processingData(y)
    x@metadata$experimentData <- experimentData(y)
    x@metadata$protocolData <- protocolData(y)
    x@metadata$qual <- qual(y)
    if (validObject(x))
        return(x)
}

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
