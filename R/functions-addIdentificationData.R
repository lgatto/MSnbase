## Each of these function does the same whether object is an MSnExp or
## and MSnSet. The different implementations process the id input
## differently, depending if it's a charater (file name), an mzRident
## object, or a mzID(Collection) object. They end up all calling the
## addIdentificationData,[MSnExp|MSnSet],data.frame method.

.addCharacterIdentificationData <-
          function(object, id, fcol, icol, acc, desc, pepseq, key,
                   verbose, ...) {
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
              addIdentificationData(object, iddf, fcol, icol, acc,
                                    desc, pepseq, verbose = verbose)
          }

.addMzRidentIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key,
             verbose, ...) {
        iddf <- as(id, "data.frame")
        iddf <- reduce(iddf, key = key)
        addIdentificationData(object, iddf, fcol, icol, acc, desc,
                              pepseq, verbose = verbose)
    }

.addMzIDIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key,
             verbose, ...) {
        iddf <- flatten(id)
        names(iddf) <- make.names(names(iddf))
        iddf <- reduce(iddf, key = key)
        addIdentificationData(object, iddf, fcol, icol, acc, desc,
                              pepseq, verbose = verbosea)
    }

.addDataFrameIdentificationData <-
    function(object, id, fcol, icol, acc, desc, pepseq, key,
             verbose, ...) {
        if (!missing(key)) { ## otherwise, id is reduced
            id <- reduce(id, key)
        }
        ## we temporaly add the spectrum.file/acquisition.number information
        ## to our fData data.frame because
        ## utils.mergeSpectraAndIdentificationData needs this information
        ## for matching
        fd <- fData(object)

        if (!nrow(fd))
            stop("No feature data found.")

        fd$spectrum.file <- basename(fileNames(object)[fromFile(object)])
        fd$acquisition.number <- acquisitionNum(object)

        fd <- utils.mergeSpectraAndIdentificationData(fd, id,
                                                            fcol = fcol,
                                                            icol = icol,
                                                            acc = acc,
                                                            desc = desc,
                                                            pepseq = pepseq)
        ## after adding the identification data we remove the
        ## temporary data to avoid duplication and problems in quantify
        ## (We don't remove acquisition.number here because the featureData slot
        ## is the only place where this information is stored in an MSnSet
        ## object; see also https://github.com/lgatto/MSnbase/issues/235.)
        cn <- colnames(fd)
        keep <- cn[cn != "spectrum.file"]
        fData(object)[, keep] <- fd[, keep, drop=FALSE]

        if (validObject(object))
            object
    }
