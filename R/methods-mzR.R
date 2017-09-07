setMethod("plotMzDelta", "mzRramp",
          function(object, reporters = NULL,
                   subset,
                   percentage = 0.1,
                   precMz = NULL,
                   precMzWidth = 2,
                   bw = 1,
                   xlim = c(40,200),
                   withLabels = TRUE,
                   size = 2.5,
                   plot = TRUE,
                   verbose = isMSnbaseVerbose()) {
              ## keep only MS2 spectra
              hd <- header(object)
              ms2 <- which(hd$msLevel == 2)
              if (!missing(subset)) {
                  if (subset <= 0 | subset >= 1) {
                      warning('subset must be in ]0, 1[. Ignoring ',
                              subset, '.', immediate. = TRUE)
                  } else {
                      n <- length(ms2)
                      .subset <- sample(n, ceiling(n * subset))
                      ms2 <- ms2[.subset]
                      if (verbose)
                          message("Subset to ", length(ms2),
                                  " spectra.")
                  }
              }               
              hd <- hd[ms2, ]
              pl <- peaksAsLists(object, ms2)  
              plotMzDelta_list(pl, reporters, percentage,
                               precMz = hd$precursorMZ,
                               precMzWidth, bw,
                               xlim, withLabels, size,
                               plot, verbose)
          })
