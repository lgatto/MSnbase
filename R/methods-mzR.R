setMethod("plotMzDelta",c("mzRramp"),
          function(object, reporters=NULL,
                   percentage=0.1,
                   precMz=NULL,
                   precMzWidth=2,
                   bw=1,
                   xlim=c(40,200),
                   withLabels=TRUE,
                   size=2.5,
                   plot=TRUE,
                   verbose=TRUE) {
              ## keep only MS2 spectra
              hd <- header(ms)
              ms2 <- which(hd$msLevel == 2)
              hd <- hd[ms2, ]
              pl <- peaksAsLists(ms, ms2)  
              plotMzDelta_list(pl, reporters, percentage,
                               precMz = hd$precursorMZ,
                               precMzWidth, bw,
                               xlim, withLabels, size,
                               plot, verbose)
          })


