##################################################################
## Initialize for Spectrum2. Use the C-constructor.
## setMethod("initialize",
##           "Spectrum2",
##           function(.Object, msLevel = 2L, peaksCount = length(mz),
##                    rt = numeric(), acquisitionNum = NA_integer_,
##                    scanIndex = integer(), tic = 0L, mz = numeric(),
##                    intensity = numeric(), fromFile = numeric(),
##                    centroided = NA, smoothed = NA,
##                    polarity = NA_integer_, merged = 1,
##                    precScanNum = NA_integer_, precursorMz = NA,
##                    precursorIntensity = NA, precursorCharge = NA_integer_,
##                    collisionEnergy = NA) {
##               .Object <- Spectrum2_mz_sorted(msLevel, peaksCount, rt,
##                                              acquisitionNum, scanIndex,
##                                              tic, mz, intensity, fromFile,
##                                              centroided, smoothed, polarity,
##                                              merged, precScanNum, precursorMz,
##                                              precursorIntensity, precursorCharge,
##                                              collisionEnergy)
##               if (validObject(.Object))
##                   .Object
##           })

## Initialize method that calls Spectrum initialize and adds the version.
setMethod("initialize",
          "Spectrum2", function(.Object, ...) {
              classVersion(.Object)["Spectrum"] <- getClassVersion("Spectrum")
              classVersion(.Object)["Spectrum2"] <- getClassVersion("Spectrum2")
              callNextMethod(.Object, ...)
          })
