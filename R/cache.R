########################################################################
## **CACHING HAS BEEN SUPERSEDED BY THE OnDiskMSnExp IMPLEMENTATION** ##
########################################################################

testCacheArg <- function(cache, maxCache = 2) {
    ## Function used to test the value of a 'cache'
    ## parameter in a read*Data function
    ## Parameters:
    ##  cache: value of the cache argument to test
    ##  maxCache: max value allowed. Generally
    ##            3, but could be less, depending
    ##            on the input data. maxCache is 1
    ##            for readMgfData for instance.
    ## Value: valid (possibly updated) cache value
    if (!is.numeric(cache))
        stop("'cache' must be numeric.")
    if (cache < 0 | cache > maxCache) {
        warning("cache must be [0:", maxCache, "]!")
        if (cache < 0) {
            warning("Setting cache to 0.")
            cache <- 0
        } else {
            warning("Setting cache to ", maxCache, ".")
            cache <- maxCache
        }
    }
    return(cache)
}

setCacheEnv <- function(toCache, level = 0, lock = TRUE) {
    ## Set the .cache slot of a pSet object.
    ## Parameters
    ##  toCache a list with
    ##     "assaydata": environment - pSet assaydata slot
    ##     "hd": header dataframe
    ##  level: numeric - cache level
    ##  lock: logical - lock env and bindings (default is TRUE)
    ## Return:
    ##  A new cache environment
    level <- testCacheArg(level)
    cacheEnv <- new.env(parent = emptyenv())
    assaydata <- toCache[["assaydata"]]
    hd <- toCache[["hd"]]
    assign("level", level, cacheEnv)
    if (level >= 1) { ## levels 2 and 3 not yet implemented
        ## precursor MZ
        precMz <- unname(eapply(assaydata, precursorMz))
        assign("rangePrecursorMz", range(precMz), cacheEnv)
        assign("nPrecursorMz", length(precMz), cacheEnv)
        assign("uPrecursorMz", length(unique(precMz)), cacheEnv)
        ## MS2 MS range
        assign("rangeMz", range(unname(eapply(assaydata, mz))),
               cacheEnv)
        ## MS2 retention time
        Rtime <- unname(eapply(assaydata, rtime))
        assign("rangeRtime", range(Rtime), cacheEnv)
        assign("nRtime", length(Rtime), cacheEnv)
        ## MS levels
        assign("msLevels", unique(unlist(eapply(assaydata, msLevel))),
               cacheEnv)
        ## precursor scans
        assign("nPrecursorScans",
               length(unique(eapply(assaydata, precScanNum))), cacheEnv)
        ## assay data size
        assign("size",
               sum(unlist(unname(eapply(assaydata, object.size)))),
               cacheEnv)
        ## full header
        assign("hd", hd, cacheEnv)
    }
    if (lock)
        lockEnvironment(cacheEnv, bindings = TRUE)
    return(cacheEnv)
}

getCacheEnv <- function(object, show = FALSE) {
    if (show)
        print(ls.str(object@.cache))
    invisible(object@.cache)
}

cacheEnvIsLocked <- function(object)
    environmentIsLocked(object@.cache)

lockCacheEnv <- function(object) {
    lockEnvironment(object@.cache, bindings = TRUE)
}
