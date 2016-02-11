## msx0 is msx version pre-1.19.12, where the Spectrum2 objects didn't
## have a polarity slot

library("MSnbase")
quantfile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
identfile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzid$")
msx <- readMSData(quantfile, verbose = FALSE)
msx <- addIdentificationData(msx, identfile)
save(msx, file = "../extdata/msx.rda")
