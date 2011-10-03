testCacheArg <- function(cache, maxCache=3) {
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
    warning("cache must be [0:",maxCache,"]!")
    if (cache < 0) {
      warning("Setting cache to 0.")
      cache <- 0
    } else {
      warning("Setting cache to ",maxCache,".")
      cache <- maxCache
    }
  }
  return(cache)
}

newCacheEnv <- function(assaydata, level=0, lock=TRUE) {
  ## Set the .cache slot of a pSet object.
  ## Parameters
  ##  assaydata: environment - pSet assaydata slot
  ##  level: numeric - cache level
  ##  lock: logical - lock env and bindings (default is TRUE)
  ## Return:
  ##  A new cache environment
  cacheEnv <- new.env()
  assign("level", level, cacheEnv)
  if (level >= 1) {
    assign("rangePrecursorMz", range(unname(eapply(assaydata,precursorMz))), cacheEnv)
    assign("rangeMz", range(unname(eapply(assaydata,mz))), cacheEnv)
    assign("rangeRtime", range(unname(eapply(assaydata,rtime))), cacheEnv)
  }
  ## levels 2 and 3 not yet implemented
  if (lock)
    lockEnvironment(cacheEnv, bindings=TRUE)
  return(cacheEnv)
}

setCacheEnv <- function(object, level, lock=TRUE) {
  ## Set the .cache slot of a pSet object.
  ## Parameters
  ##  object: MSnExp instance
  ##  level: numeric - cache level
  ##  lock: logical - lock env and bindings (default is TRUE)
  ## Return:
  ##  Invisibly returns the updated env.
  ## Side effect:
  ##  The object@.cache is updated
  object@.cache <- new.env()
  assign("level", level, object@.cache) 
  if (level >= 1) {
    assign("rangePrecursorMz", range(precursorMz(object)), object@.cache)
    assign("rangeMz", range(mz(object)), object@.cache)
    assign("rangeRtime", range(rtime(object)), object@.cache)
  }
  ## levels 2 and 3 not yet implemented
  if (lock)
    lockCacheEnv(object)
  invisible(validObject(object))
}

getCacheEnv <- function(object) {
  print(ls.str(object@.cache))
  invisible(object@.cache)
}

cacheEnvIsLocked <- function(object)
  environmentIsLocked(object@.cache)

lockCacheEnv <- function(object) {
  lockEnvironment(object@.cache, bindings=TRUE)  
}

## updateCacheEnv <- function(object) {
##   cacheEnv <- object@.cache
##   cacheLevel <- cacheEnv$level
##   newCacheEnv <- setCacheEnv(level = cacheLevel)
## }
