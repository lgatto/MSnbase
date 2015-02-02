count_MSnSet <- function(object) {
    object <- removeNoId(object)
    .exprs <- matrix(1, ncol = 1, nrow = length(object))
    rownames(.exprs) <- featureNames(object)
    colnames(.exprs) <- sampleNames(object)
    .qual <- data.frame()

    fd <- header(object)
    if (nrow(fData(object)) > 0) {
        if (nrow(fData(object)) == length(object)) {
            fd <- combine(fData(object), fd)
        } else {
            warning("Unexpected number of features in featureData slot. Dropping it.")
        }
    }
    ## featureData rows must be reordered to match assayData rows
    .featureData <- new("AnnotatedDataFrame",
                        data = fd[rownames(.exprs), ])
    msnset <- new("MSnSet",
                  qual = .qual,
                  exprs = .exprs,
                  experimentData = experimentData(object),
                  phenoData = phenoData(object),
                  featureData = .featureData,
                  annotation = "No annotation")

    msnset@processingData <- object@processingData
    if (validObject(msnset))
        return(msnset)
}

tic_MSnSet <- function(object) {
    .exprs <- matrix(tic(object), ncol = 1)
    rownames(.exprs) <- featureNames(object)
    colnames(.exprs) <- sampleNames(object)
    .qual <- data.frame()

    fd <- header(object)
    if (nrow(fData(object)) > 0) {
        if (nrow(fData(object)) == length(object)) {
            fd <- combine(fData(object), fd)
        } else {
            warning("Unexpected number of features in featureData slot. Dropping it.")
        }
    }
    ## featureData rows must be reordered to match assayData rows
    .featureData <- new("AnnotatedDataFrame", data=fd[rownames(.exprs), ])
    msnset <- new("MSnSet",
                  qual = .qual,
                  exprs = .exprs,
                  experimentData = experimentData(object),
                  phenoData = phenoData(object),
                  featureData = .featureData,
                  annotation = "No annotation")

    msnset@processingData <- object@processingData
    if (validObject(msnset))
        return(msnset)
}



## SI Spectra Index
## Griffin et al. 2010, PMID: 20010810
SI <- function(object,
               method = c("SI", "SIgi", "SIn"),
               groupBy = "accession",
               plength = "length",
               verbose = TRUE) {
    method <- match.arg(method)
    if (!plength %in% fvarLabels(object))
        stop(plength,
             " not found in fvarLabel(.). 'plength' must a feature variable")
    if (!groupBy %in% fvarLabels(object))
        stop(groupBy,
             " not found in fvarLabel(.). 'groupBy' must a feature variable")

    object <- tic_MSnSet(object)
    groupBy <- as.factor(fData(object)[, groupBy])

    ## SI: protein-wise summed tic
    object <- combineFeatures(object, groupBy = groupBy,
                              fun = "sum", verbose = verbose)

    if (method %in% c("SIgi", "SIn"))
        exprs(object) <- exprs(object)/colSums(exprs(object))

    if (method == "SIn") {
        plength <- fData(object)[, plength]
        exprs(object) <- exprs(object)/plength
    }

    object <- logging(object, paste0("Quantification by ", method))
    if (validObject(object))
        return(object)
}

## (N)SAF (Normalised) Spectral Abundance Factor
## Paoletti et al. 2006, PMID: 17138671
SAF <- function(object,
                method = c("SAF", "NSAF"),
                groupBy = "accession",
                plength = "length",
                verbose = TRUE) {
    method <- match.arg(method)
    if (!plength %in% fvarLabels(object))
        stop(plength,
             " not found in fvarLabel(.). 'plength' must a feature variable")
    if (!groupBy %in% fvarLabels(object))
        stop(groupBy,
             " not found in fvarLabel(.). 'groupBy' must a feature variable")
    object <- count_MSnSet(object)

    groupBy <- as.factor(fData(object)[, groupBy])

    object <- combineFeatures(object,
                              groupBy = groupBy,
                              fun = length, verbose = verbose)

    plength <- fData(object)[, plength]
    exprs(object) <- exprs(object)/plength

    if (method == "NSAF")
        exprs(object) <- exprs(object)/colSums(exprs(object))

    object <- logging(object, paste0("Quantification by ", method))
    if (validObject(object))
        return(object)
}
