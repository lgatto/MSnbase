############################################################
## Methods and functions for class ProcessingStep

############################################################
## ProcessingStep
##
## Constructor
ProcessingStep <- function(FUN=character(), ARGS=list()){
    if(missing(FUN))
        FUN <- character()
    return(new("ProcessingStep", FUN=FUN, ARGS=ARGS))
}

############################################################
## show
##
setMethod("show", "ProcessingStep", function(object){
    cat("Object of class \"",class(object),"\"\n",sep="")
    cat(" Function: ", object@FUN, "\n", sep="")
    args <- object@ARGS
    if(length(args) > 0){
        cat(" Arguments:\n")
        for(i in 1:length(args)){
            cat("  o ", names(args)[i], "= ", paste(args[[i]], collapse=", "), "\n", sep="")
        }
    }
})

############################################################
## execute
##
## Execute the processing step.
setMethod("execute", "ProcessingStep", function(object, ...){
    return(do.call(object@FUN, args=c(list(...), object@ARGS)))
})
