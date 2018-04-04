##' Reads as set of `mzId` files containing PSMs an generates a
##' `data.frame`.
##'
##' This function uses the functionality provided by the `mzR` package
##' to access data in the `mzId` files. An object of class `mzRident`
##' can also be coerced to a `data.frame` using `as(, "data.frame")`.
##'
##' @title Import peptide-spectrum matches
##' @param files A `character` of `mzid` files.
##' @return A `data.frame` containing the PSMs stored in the `mzId`
##'     files.
##' @md
##' @aliases coerce,mzRident,data.frame-method
##' @author Laurent Gatto
##' @seealso [filterIdentificationDataFrame()] to filter out unreliable PSMs.
##' @examples
##' idf <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid"
##' f <- msdata::ident(full.names = TRUE, pattern = idf)
##' basename(f)
##' readMzIdData(f)
readMzIdData <- function(files) {
    if (length(files) == 1 && file.exists(files)) {
        id <- mzR::openIDfile(files)
        iddf <- as(id, "data.frame")
    } else {
        if (!all(flex <- file.exists(files)))
            stop(paste(files[!flex], collapse = ", "), " not found.")
        iddf <- lapply(files,
                       function(f) as(openIDfile(f), "data.frame"))
        iddf <- do.call(rbind, iddf)
    }
    iddf
}

##' A function to filter out PSMs matching to the decoy database, of
##' rank greater than one and matching non-proteotypic peptides.
##'
##' The PSMs should be stored in a `data.frame` such as those produced
##' by [readMzIdData()]. Note that this function should be called
##' before calling the [reduce][reduce,data.frame-method] method on a
##' PSM `data.frame`.
##'
##' @md
##' @title Filter out unreliable PSMs.
##' @param x A `data.frame` containing PSMs.
##' @param decoy The column name defining whether entries match the
##'     decoy database. Default is `"isDecoy"`. The column should be a
##'     `logical` and only PSMs holding a `FALSE` are
##'     retained. Ignored is set to `NULL`.
##' @param rank The column name holding the rank of the PSM. Default
##'     is `"rank"`. This column should be a `numeric` and only PSMs
##'     having rank equal to 1 are retained. Ignored is set to `NULL`.
##' @param accession The column name holding the protein (groups)
##'     accession. Default is `"DatabaseAccess"`. Ignored is set to
##'     `NULL`.
##' @param spectrumID The name of the spectrum identifier
##'     column. Default is `spectrumID`.
##' @param verbose A `logical` verbosity flag. Default is to take
##'     `isMSnbaseVerbose()`.
##' @return A new `data.frame` with filtered out peptides and with the
##'     same columns as the input `x`.
##' @author Laurent Gatto
filterIdentificationDataFrame <- function(x,
                                          decoy = "isDecoy",
                                          rank = "rank",
                                          accession = "DatabaseAccess",
                                          spectrumID = "spectrumID",
                                          verbose = isMSnbaseVerbose()) {
    n0 <- nrow(x)
    if (verbose)
        message("Starting with ", n0, " PSMs:")
    if (!is.null(decoy)) {
        x <- x[!x[, decoy], ]
        n1 <- nrow(x)
        if (verbose)
            message(" removed ", n0 - n1, " decoy hits")
        n0 <- n1
    }
    if (!is.null(rank)) {
        x <- x[x[, rank] == 1, ]
        n2 <- nrow(x)
        if (verbose)
            message(" removed ", n0 - n2, " PSMs with rank > 1")
        n0 <- n2
    }
    if (!is.null(accession)) {
        keep <- tapply(x[, accession],
                       x[, spectrumID],
                       function(xx) length(unique(xx))) == 1
        x <- x[keep, ]
        n3 <- nrow(x)
        if (verbose)
            message(" removed ", n0 - n3, " non-proteotypic peptides")
    }
    if (verbose)
        message(nrow(x), " PSMs left.")
    x
}


## Each of these function does the same whether object is an MSnExp or
## and MSnSet. The different implementations process the id input
## differently, depending if it's a charater (file name), an mzRident
## object, or a mzID(Collection) object. They end up all calling the
## .addDataFrameIdentificationData.

.addCharacterIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key, decoy,
             rank, accession, verbose, ...) {
        ## The code below is the same as in readMzIdData but filters
        ## and reduces after each mzId file is read and converted to a
        ## data.frame.
        if (length(id) == 1 && file.exists(id)) {
            id <- mzR::openIDfile(id)
            iddf <- as(id, "data.frame")
            iddf <- filterIdentificationDataFrame(iddf, decoy = decoy,
                                                  rank = rank,
                                                  accession = accession,
                                                  verbose = verbose,
                                                  ...)
            iddf <- reduce(iddf, key = key)
        } else {
            if (!all(flex <- file.exists(id)))
                stop(paste(id[!flex], collapse = ", "), " not found.")
            iddf <- lapply(id,
                           function(x) {
                               iddf <- as(openIDfile(x), "data.frame")
                               iddf <- filterIdentificationDataFrame(iddf, decoy = decoy,
                                                                     rank = rank,
                                                                     accession = accession,
                                                                     verbose = verbose,
                                                                     ...)
                               iddf <- reduce(iddf, key = key)
                           })
            iddf <- do.call(rbind, iddf)
        }
        ## Filtering already done - set these args to NULL
        .addDataFrameIdentificationData(object, iddf, fcol, icol, acc,
                                        desc, pepseq, decoy = NULL,
                                        rank = NULL, accession = NULL,
                                        verbose = verbose,
                                        ...)
    }

.addMzRidentIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key, decoy,
             rank, accession, verbose, ...) {
        iddf <- as(id, "data.frame")
        iddf <- filterIdentificationDataFrame(iddf, decoy = decoy,
                                              rank = rank, accession = accession,
                                              verbose = verbose,
                                              ...)
        iddf <- reduce(iddf, key = key)
        ## Filtering already done - set these args to NULL
        .addDataFrameIdentificationData(object, iddf, fcol, icol, acc,
                                        desc, pepseq, decoy = NULL,
                                        rank = NULL, accession = NULL,
                                        verbose = verbose)
    }

.addMzIDIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key, decoy,
             rank, accession, verbose, ...) {
        iddf <- flatten(id)
        names(iddf) <- make.names(names(iddf))
        iddf <- filterIdentificationDataFrame(iddf, decoy = decoy,
                                              rank = rank, accession = accession,
                                              spectrumID = "spectrumid",
                                              verbose = verbose,
                                              ...)
        iddf <- reduce(iddf, key = key)
        ## Filtering already done - set these args to NULL
        .addDataFrameIdentificationData(object, iddf, fcol, icol, acc,
                                        desc, pepseq, decoy = NULL,
                                        rank = NULL, accession = NULL,
                                        verbose = verbose)
    }

.addDataFrameIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key, decoy,
             rank, accession, verbose, ...) {
        if (!missing(key)) { ## otherwise, id is reduced
            id <- reduce(id, key)
        }
        ## Filtering arguments could all be NULL to bypass filtering,
        ## as in the other .addClassIdentificationData function above
        id <- filterIdentificationDataFrame(id, accession = accession,
                                            rank = rank,
                                            decoy = decoy,
                                            verbose = verbose,
                                            ...)
        ## we temporaly add the spectrum.file/acquisition.number information
        ## to our fData data.frame because
        ## utils.mergeSpectraAndIdentificationData needs this information
        ## for matching
        fd <- fData(object)

        if (!nrow(fd))
            stop("No feature data found.")

        fd$spectrum.file <- basename(fileNames(object)[fromFile(object)])
        fd$acquisition.number <- acquisitionNum(object)
        fd <- utils.mergeSpectraAndIdentificationData(fd, id,
                                                      fcol = fcol,
                                                      icol = icol,
                                                      acc = acc,
                                                      desc = desc,
                                                      pepseq = pepseq,
                                                      rank)
        ## after adding the identification data we remove the
        ## temporary data to avoid duplication and problems in quantify
        ## (We don't remove acquisition.number here because the featureData slot
        ## is the only place where this information is stored in an MSnSet
        ## object; see also https://github.com/lgatto/MSnbase/issues/235.)
        cn <- colnames(fd)
        keep <- cn[cn != "spectrum.file"]
        fData(object)[, keep] <- fd[, keep, drop=FALSE]

        if (validObject(object))
            object
    }
