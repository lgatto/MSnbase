############################################################
## Methods for "basic" classes.
##

############################################################
## Methods for matrix.
##
############################################################
## removePeaks
##
## The removePeaks method for a matrix. Thus we can apply it
## to the matrix returned by mzR::peaks
setMethod("removePeaks", "matrix", function(object, t, centroided=FALSE){
    ## Second column is intensity.
    if(t == "min")
        t <- min(object[object[, 2] > 0, 2], na.rm=TRUE)
    if(!is.numeric(t))
        stop("'t' must either be 'min' or numeric.")
    if(centroided){
        object[, 2] <- utils.removePeaks_centroided(object[, 2], t)
    }else{
        object[, 2] <- utils.removePeaks(object[, 2], t)
    }
    return(object)
})
