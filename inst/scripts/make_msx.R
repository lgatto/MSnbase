quantfile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
identfile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzid$")
msx <- readMSData(quantfile, verbose = FALSE)
msx <- addIdentificationData(msx, identfile)
save(msx, file = "../extdata/msx.rda")
