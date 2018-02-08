################################################
## Quantitation parmeter constructors


##' @rdname QuantitationParam-class
##'
##' @param reporters See class slot.
##' @param msLevel See class slot.
##' @param method See class slot.
IsobaricTagging <- function(reporters,
                            msLevel, 
                            method =
                                c("max",
                                  "trapezoidation",
                                  "sum")) {
    if (missing(reporters))
        stop("Please provide the isobaric reporters.")
    if (missing(msLevel)) {
        if (isMSnbaseVerbose())
            message("No MS level provided. Will be determined from data.")
        msLevel <- integer()
    }
    method <- match.arg(method)
    .QuantitationParam(msLevel = as.integer(msLevel),
                       reporters = reporters,
                       method = method,
                       wd = width(reporters),
                       strict = FALSE,
                       name = "IsobaricTagging")
}

##' @rdname QuantitationParam-class
##' @param pepseq See class slot.
##' @param dbaccess See class slot.
##' @param plength See class slot.
##' @param nprot See class slot.
SpectralCounting <- function(method =
                                 c("count", "tic",
                                   "SI", "SIgi",
                                   "SIn", "SAF",
                                   "NSAF"),
                             pepseq = "sequence", 
                             dbaccess = "DatabaseAccess",
                             plength = "DBseqLength",
                             nprot = "nprot") {
    method <- match.arg(method)
    .QuantitationParam(msLevel = 2L,
                       method = method,
                       name = "SpectralCounting",
                       pepseq = pepseq,     ## needed for count
                       dbaccess = dbaccess, ## needed for SI* and *SAF
                       plength = plength,   ## needed for SI* and *SAF
                       nprot = nprot)       ## needed for SI* and *SAF
}

##' @rdname QuantitationParam-class
##'
##' @param object Object of class `QuantiationParam`.
setMethod("show", "QuantitationParam",
          function(object) {
              cat("Quantitation method: '", object@name, "'\n", sep = "")
              cat(" Method:", object@method, "\n")
              if (length(object@msLevel))
                  cat(" MS level:", object@msLevel, "\n")
              if (!is.null(object@reporters))
                  cat(" Reporters:", names(object@reporters), "\n")
          })
