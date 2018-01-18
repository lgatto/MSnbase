setClassUnion("ReporterIonsOrNull", c("ReporterIons", "NULL"))

.QuantitationParam <-
    setClass("QuantitationParam",
         slots = c(msLevel = "integer",
                   reporters = "ReporterIonsOrNull",
                   method = "character",
                   methargs = "list",
                   name = "character"),
         contains = "Versioned",
         prototype = prototype(
             name = "Quantitation method",
             new("Versioned",
                 versions = c(classVersion("ReporterIons"),
                              QuantitationParam = "0.1.0"))
         ))


setMethod("show", "QuantitationParam",
          function(object) {
              cat("Quantitation method: '", object@name, "'\n", sep = "")
              cat(" Method:", object@method, "\n")
              if (length(object@msLevel))
                  cat(" MS level:", object@msLevel, "\n")
              if (!is.null(object@reporters))
                  cat(" Reporters:", names(object@reporters), "\n")
          })

IsobaricQuantitation <- function(reporters,
                                 msLevel, 
                                 method = c("max",
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
    .QuantitationParam(msLevel = msLevel,
                       reporters = reporters,
                       method = method,
                       name = "Isobaric labelling")
}

SpectralCountingQuantitation <- function(method = c("count", "SI",
                                                    "SIgi", "SIn",
                                                    "SAF")) {
    method <- match.arg(method)
    .QuantitationParam(msLevel = 2L,
                       reporters = NULL,
                       method = method,
                       name = "Spectral counting")
}


quantify2 <- function(object, 
                      params,
                      ...) {
    stopifnot(inherits(object, "OnDiskMSnExp"))
    stopifnot(inherits(params, "QuantitationParam"))
    ## hard-coding isobaric quantitation
    if (!length(params@msLevel))
        params@msLevel <- max(msLevel(object))
    obj2 <- filterMsLevel(object, params@msLevel)
    res <- quantify(obj2,
                    method = params@method,
                    reporters = params@reporters,
                    ...)
    res
}
