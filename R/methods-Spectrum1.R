##################################################################
## Initialize for Spectrum1.
setMethod("initialize",
          "Spectrum1",
          function(.Object, msLevel = 1L, peaksCount = length(mz),
                   rt = numeric(), acquisitionNum = NA_integer_,
                   scanIndex = integer(), tic = 0, mz = numeric(),
                   intensity = numeric(), fromFile = numeric(),
                   centroided = NA, smoothed = NA,
                   polarity = NA_integer_) {
              res <- Spectrum1_mz_sorted(peaksCount = peaksCount, rt = rt,
                                         acquisitionNum = acquisitionNum,
                                         scanIndex = scanIndex,
                                         tic = tic, mz = mz,
                                         intensity = intensity,
                                         fromFile = fromFile,
                                         centroided = centroided,
                                         smoothed = smoothed,
                                         polarity = polarity)
              if (validObject(res))
                  res
          })
