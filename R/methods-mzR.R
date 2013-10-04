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
                   verbose = TRUE) {             
              ## keep only MS2 spectra
              hd <- header(ms)
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
              pl <- peaksAsLists(ms, ms2)  
              plotMzDelta_list(pl, reporters, percentage,
                               precMz = hd$precursorMZ,
                               precMzWidth, bw,
                               xlim, withLabels, size,
                               plot, verbose)
          })


setMethod("chromatogram", "mzRramp",
          function(object,
                   y = c("tic", "bpi"),
                   legend = TRUE,
                   plot = TRUE,
                   ...) {                    
          hd <- header(object)          
          .chromatogram(hd, y, legend, plot, ...)
      })

setMethod("chromatogram", "data.frame",
          function(object,
                   y = c("tic", "bpi"),
                   legend = TRUE,
                   plot = TRUE,
                   ...) {          
          stopifnot("retentionTime" %in% colnames(object))
          y <- match.arg(y)
          chck <- switch(y,
                         tic = stopifnot("totIonCurrent" %in% colnames(object)),
                         bpi = stopifnot("basePeakIntensity" %in% colnames(object)))          
          .chromatogram(object, y, legend, plot, ...)
      })


setMethod("xic", "mzRramp",
          function(object, ...) xic_1(object, ...))
