setMethod("calculateFragments", c("character", "missing"),
          function(sequence, type=c("b", "y"), z=1,
                   modifications=c(C=160.030649), verbose=TRUE) {
            l <- lapply(sequence, .calculateFragments,
                        type=type, z=z, modifications=modifications,
                        verbose=verbose)
            return(do.call(rbind, l))
        })
