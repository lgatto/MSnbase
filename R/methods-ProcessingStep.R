############################################################
## Methods and functions for class ProcessingStep

############################################################
## ProcessingStep
##
## Constructor
ProcessingStep <- function(FUN = character(), ARGS = list())  {
    if (missing(FUN))
        FUN <- character()
    return(new("ProcessingStep", FUN = FUN, ARGS = ARGS))
}

############################################################
## show
##
setMethod("show", "ProcessingStep", function(object) {
    cat("Object of class \"", class(object), "\"\n", sep = "")
    cat(" Function: ", object@FUN, "\n", sep = "")
    args <- object@ARGS
    if (length(args) > 0) {
        cat(" Arguments:\n")
        for (i in 1:length(args)) {
            cat("  o ", names(args)[i], " = ",
                paste(args[[i]], collapse = ", "), "\n", sep = "")
        }
    }
})

############################################################
## execute
##
## Eventually make a method execute, ProcessingStep, Spectrum

############################################################
## executeProcessingStep
##
## Execute the processing step.
executeProcessingStep <- function(object, ...) {
    if (!is(object, "ProcessingStep"))
        stop("'object' is supposed to be a 'ProcessingStep' object!")
    ## Eventually switch based on the MSlevel?
    ## Check if we've got an msLevel argument in object@ARGS, if so,
    ## check that msLevel(x) matches any of the msLevel, if not,
    ## return x, otherwise call do.call.
    return(do.call(object@FUN, args = c(list(...), object@ARGS)))
}

#' Internal function to apply the lazy processing queue to each spectrum
#' in the provided list.
#'
#' @param x `list` of `Spectrum` objects.
#'
#' @param queue `list` (or `NULL`) of `ProcessingStep` objects.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.apply_processing_queue <- function(x, queue = NULL) {
    if (length(queue)) {
        x <- lapply(x, function(z, q) {
            for (pStep in q) {
                z <- executeProcessingStep(pStep, z)
            }
            z
        }, q = queue)
    }
    x
}
