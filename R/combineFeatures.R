combineFeatures <- function(object,
                            groupBy,
                            fun = c("mean",
                                "median",
                                "weighted.mean",
                                "sum",
                                "medpolish",
                                "iPQF",
                                "NTR"),
                            fcol,                            
                            redundancy.handler = c("unique", "multiple"),
                            cv = TRUE,
                            cv.norm = "sum",
                            verbose = isMSnbaseVerbose(),
                            ... ## further arguments to fun
                            ) {
    if (is.character(fun))
        fun <- match.arg(fun)
    if (missing(groupBy)) {
        if (missing(fcol))
            stop("Require either 'groupBy' or 'fcol'.")
        stopifnot(fcol %in% fvarLabels(object))
        groupBy <- fData(object)[, fcol]
    }
    if (is.list(groupBy)) {
        if (length(groupBy) != nrow(object))
            stop("'length(groupBy)' must be equal to 'nrow(object)': ",
                 length(groupBy), " != ", nrow(object), ".")
        if (!is.null(names(groupBy))) {
            if (!all(names(groupBy) %in% featureNames(object)))
                stop("'groupBy' names and 'featureNames(object)' don't match.")
            if (!all(names(groupBy) == featureNames(object))) {
                warning("Re-ordering groupBy to match feature names.")
                groupBy <- groupBy[featureNames(object)]
            }
        }
        redundancy.handler <- match.arg(redundancy.handler)
        result <- combineFeaturesL(object, groupBy, fun,
                                   redundancy.handler,
                                   cv, cv.norm, verbose, ...)
    } else { ## factor, numeric or character
        result <- combineFeaturesV(object, groupBy, fun,
                                   cv, cv.norm, verbose, ...)
    }
    if (validObject(result))
        return(result)
}

combineFeaturesL <- function(object,   ## MSnSet
                             groupBy,  ## list
                             fun,
                             redundancy.handler,
                             cv = TRUE,
                             cv.norm = "sum",
                             verbose = isMSnbaseVerbose(),
                             ...    ## additional arguments to fun
                             ) {
    ## handling of the redundancy
    if (redundancy.handler == "multiple") {
        expansion.index <- rep(seq_len(nrow(object)), sapply(groupBy, length))
        new.exprs <- exprs(object)[expansion.index, , drop = FALSE]
        rownames(new.exprs) <- NULL
        new.feature.data <- fData(object)[expansion.index, , drop = FALSE]
        rownames(new.feature.data) <- NULL
        ## TODO: check this
        object <- new("MSnSet", exprs = new.exprs,
                      featureData = new(
                          "AnnotatedDataFrame",
                          data = new.feature.data),
                      phenoData = phenoData(object))
        groupBy <- unlist(groupBy)
    } else { ## redundancy.handler == "unique" - checked by match.arg
        idx.unique <- sapply(groupBy, length) < 2
        object <- object[idx.unique, ]
        groupBy <- unlist(groupBy[idx.unique])
    }
    combineFeaturesV(object, groupBy, fun, cv, cv.norm, verbose, ...)
}


combineFeaturesV <- function(object,   ## MSnSet
                             groupBy,  ## factor, character or numeric
                             fun,
                             cv = TRUE,
                             cv.norm = "sum",
                             verbose = isMSnbaseVerbose(),
                             ...    ## additional arguments to fun
                             ) {
    groupBy <- as.character(groupBy)
    if (cv) {
        cv.mat <- featureCV(object, groupBy = groupBy,
                            norm = cv.norm)
    }
    n1 <- nrow(object)
    ## !! order of features in matRes is defined by the groupBy factor !!
    if (is.character(fun) && fun == "iPQF") {
        ## NB: here, we pass the object, not only assay data,
        ##     because iPGF also needs the feature data, otherwise
        ##     not passed and used in combineFeatureMatrix
        ##     iPQF still returns a matrix, though.
        matRes <- iPQF(object, groupBy, ...)
    } else if (is.character(fun) && fun == "NTR") {
        matRes <- normToReference(exprs(object), group=groupBy, ...)
        ## order matrix according to groupBy (is important because rownames are
        ## overwritten a few lines below
        matRes <- matRes[order(unique(groupBy)), , drop=FALSE]
    } else {
        matRes <- as.matrix(combineMatrixFeatures(exprs(object),
                                                  groupBy, fun,
                                                  verbose = verbose,
                                                  ...))
    }
    fdata <- fData(object)[!duplicated(groupBy), , drop = FALSE] ## takes the first occurences
    fdata <- fdata[order(unique(groupBy)), , drop = FALSE] ## ordering fdata according to groupBy factor
    rownames(matRes) <- rownames(fdata)
    colnames(matRes) <- sampleNames(object)
    if (cv)
        fdata <- cbind(fdata, cv.mat)
    res <- new("MSnSet", exprs = matRes,
               featureData = new("AnnotatedDataFrame",
                                 data = fdata))
    res@processingData@merged <- TRUE
    res@qual <- object@qual[0, ]
    pData(res) <- pData(object)
    if (is.character(fun)) {
        msg <- paste("Combined ", n1, " features into ",
                     nrow(res), " using ", fun, sep = "")
    } else {
        msg <- paste("Combined ", n1, " features into ",
                     nrow(res), " using user-defined function",
                     sep = "")
    }
    if (verbose)
        message(msg)
    res@processingData@processing <- c(object@processingData@processing,
                                       paste(msg, ": ", date(), sep = ""))
    ## update feature names according to the groupBy argument
    ## new in version 1.5.9
    fn <- sort(unique(groupBy))
    featureNames(res) <- fn
    if (validObject(res))
        return(res)
}

combineMatrixFeatures <- function(matr,    ## matrix
                                  groupBy, ## char/factor
                                  fun,
                                  verbose = isMSnbaseVerbose(),
                                  ...) { ## additional arguments to fun
    if (is.character(fun)) {
        ## Using a predefined function
        if (fun == "medpolish") {
            summarisedFeatures <- by(matr,
                                     groupBy,
                                     function(x) {
                                         medpol <- medpolish(as.numeric(x),
                                                             trace.iter = verbose, ...)
                                         return(medpol$overall + medpol$col)
                                     })
        } else if (fun == "weighted.mean") {
            ## Expecting 'w' argument
            args <- list(...)
            if (is.null(args$w))
                stop("Expecting a weight parameter 'w' for 'weigthed mean'.")
            w <- args$w
            if (is.null(colnames(matr)))
                colnames(matr) <- paste0("X", 1:ncol(matr))
            summarisedFeatures <- apply(matr,2,
                                        function(x) {
                                            .data <- data.frame(x=x, groupBy, w=w)
                                            ddply(.data,
                                                  "groupBy",
                                                  summarise,
                                                  wmn = weighted.mean(x,w))
                                        })
            summarisedFeatures <- do.call(cbind, as.list(summarisedFeatures))
            rn <- summarisedFeatures[,1]
            summarisedFeatures <-
                summarisedFeatures[, grep("wmn", colnames(summarisedFeatures))]
            colnames(summarisedFeatures) <- colnames(matr)
            rownames(summarisedFeatures) <- rn
            return(summarisedFeatures)
        } else {
            ## using either 'sum', 'mean', 'median'
            summarisedFeatures <- by(matr,
                                     groupBy,
                                     function(x) apply(x, 2, eval(parse(text = fun)), ...))
        }
    } else {
        ## using user-defined function
        summarisedFeatures <- by(matr,
                                 groupBy,
                                 function(x) apply(x, 2, fun, ...))
    }
    return(do.call(rbind, as.list(summarisedFeatures)))
}


##' This function evaluates the variability within all protein group
##' of an \code{MSnSet}. If a protein group is composed only of a
##' single feature, \code{NA} is returned.
##'
##' This function can be used to identify protein groups with
##' incoherent feature (petides or PSMs) expression patterns. Using
##' \code{max} as a function, one can identify protein groups with
##' single extreme outliers, such as, for example, a mis-identified
##' peptide that was erroneously assigned to that protein group. Using
##' \code{mean} identifies more systematic inconsistencies where, for
##' example, the subsets of peptide (or PSM) feautres correspond to
##' proteins with different expression patterns.
##'
##' @title Identify aggregation outliers
##' @param object An object of class \code{MSnSet}.
##' @param groupBy A \code{character} containing the protein grouping
##'     feature variable name.
##' @param fun A function the summarise the distance between features
##'     within protein groups, typically \code{max} or
##'     \code{mean}.\code{median}.
##' @return A \code{matrix} providing the number of features per
##'     protein group (\code{nb_feats} column) and the aggregation
##'     summarising distance (\code{agg_dist} column).
##' @author Laurent Gatto
##' @seealso \code{\link{combineFeatures}} to combine PSMs
##'     quantitation into peptides and/or into proteins.
##' @examples
##' library("pRolocdata")
##' data(hyperLOPIT2015ms3r1psm)
##' groupBy <- "Protein.Group.Accessions"
##' res1 <- aggvar(hyperLOPIT2015ms3r1psm, groupBy, fun = max)
##' res2 <- aggvar(hyperLOPIT2015ms3r1psm, groupBy, fun = mean)
##' par(mfrow = c(1, 3))
##' plot(res1, log = "y", main = "Single outliers (max)")
##' plot(res2, log = "y", main = "Overall inconsistency (mean)")
##' plot(res1[, "agg_dist"], res2[, "agg_dist"],
##'      xlab = "max", ylab = "mean")
aggvar <- function(object, groupBy, fun) {
    stopifnot(inherits(object, "MSnSet"))
    stopifnot(groupBy %in% fvarLabels(object))
    .d <- function(x, fun. = fun) {
        dx <- dist(as.matrix(x))
        if (length(dx) == 0) return(NA)
        do.call(fun., list(dx))
    }
    gb <- fData(object)[, groupBy]
    nr <- by(exprs(object), gb, nrow)
    d <-  by(exprs(object), gb, .d, fun)
    ans <- cbind(agg_dist = d, nb_feats = nr)
    attr(ans, "agg_dist") <- fun
    ans
}
