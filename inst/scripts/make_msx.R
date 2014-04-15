file <- dir(system.file(package = "MSnbase",dir = "extdata"),
            full.name = TRUE,pattern = "mzXML$")
msx <- readMSData(file, verbose = FALSE)
save(msx, file = "../extdata/msx.rda")
