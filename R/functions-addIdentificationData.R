## Each of these function does the same whether object is an MSnExp or
## and MSnSet. The different implementations process the id input
## differently, depending if it's a charater (file name), an mzRident
## object, or a mzID(Collection) object. They end up all calling the
## addIdentificationData,[MSnExp|MSnSet],data.frame method.

.addCharacterIdentificationData <-
    function(object, id,
             fcol = c("spectrum.file", "acquisition.number"),
             icol = c("spectrumFile", "acquisitionNum"), 
             acc = "DatabaseAccess",
             desc = "DatabaseDescription",
             pepseq = "sequence",
             key = "spectrumID",
             verbose = isMSnbaseVerbose(),
             ...) {
        if (length(id) == 1 && file.exists(id)) {
            id <- mzR::openIDfile(id)
            iddf <- as(id, "data.frame")
            iddf <- reduce(iddf, key = key)
        } else {
            if (!all(flex <- file.exists(id)))
                stop(paste(id[!flex], collapse = ", "), " not found.")
            iddf <- lapply(id,
                           function(x) {
                               iddf <- as(openIDfile(x), "data.frame")
                               iddf <- reduce(iddf, key = key)
                           })
            iddf <- do.call(rbind, iddf)
        }
        addIdentificationData(object, iddf, fcol, icol, acc, desc,
                              pepseq, verbose)
    }

.addMzRidentIdentificationData <-
    function(object, id,
             fcol = c("spectrum.file", "acquisition.number"),
             icol = c("spectrumFile", "acquisitionNum"),
             acc = "DatabaseAccess",
             desc = "DatabaseDescription",
             pepseq = "sequence",
             key = "spectrumID",
             verbose = isMSnbaseVerbose(),
             ...) {
        iddf <- as(id, "data.frame")
        iddf <- reduce(iddf, key = key)
        addIdentificationData(object, iddf, fcol, icol, acc, desc,
                              pepseq, verbose)
    }

.addMzIDIdentificationData <-
    function(object, id,
             fcol = c("spectrum.file", "acquisition.number"),
             icol = c("spectrumFile", "acquisitionnum"),
             acc = "accession",
             desc = "description",
             pepseq = "pepseq",
             key = "spectrumid",
             key = "spectrumid",
             ...) {
        iddf <- flatten(id)
        iddf <- reduce(iddf, key = key)              
        addIdentificationData(object, iddf, fcol, icol, acc, desc,
                              pepseq, verbosea)
    }

## There is no common .addDataFrameIdentificationData function as he
## implementations for MSnExp and MSnSet objects are different. Might
## need to be revised though -- see FIXME in
## addIdentificationData,MSnSet,data.frame.
