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


setMethod("chromatogram", "character",
          function(object,
                   y = c("tic", "bpi"),
                   legend = TRUE,
                   plot = TRUE,
                   ms = 1L,
                   ...) {
              object <- openMSfile(object)
              on.exit(close(object))
              hd <- header(object)
              f <- basename(fileName(object))
              chromatogram(hd, y, f, legend, plot, ms, ...)
          })

setMethod("chromatogram", "mzRramp",
          function(object,
                   y = c("tic", "bpi"),
                   legend = TRUE,
                   plot = TRUE,
                   ms = 1L,
                   ...) {                    
          hd <- header(object)
          f <- basename(fileName(object))
          chromatogram(hd, y, f, legend, plot, ms, ...)
      })

setMethod("chromatogram", "data.frame",
          function(object,
                   y = c("tic", "bpi"),
                   f,
                   legend = TRUE,
                   plot = TRUE,
                   ms = 1L,
                   ...) { 
          stopifnot("retentionTime" %in% colnames(object))
          stopifnot("msLevel" %in% colnames(object))
          y <- match.arg(y)
          chck <- switch(y,
                         tic = stopifnot("totIonCurrent" %in% colnames(object)),
                         bpi = stopifnot("basePeakIntensity" %in% colnames(object)))
          object <- object[object$msLevel == ms, ]
          if (nrow(object) == 0)
              stop("No spectra of level ", ms, " found.")
          .chromatogram(object, y, f, legend, plot, ...)
      })


setMethod("xic", "mzRramp",
          function(object, ...) xic_1(object, ...))

setMethod("xic", "character",
          function(object, ...) {
              object <- openMSfile(object)
              xic_1(object, ...)
          })
