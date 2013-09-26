setMethod("writeMgfData",
          signature = signature("Spectrum"),
          function(object,
                   con = NULL,
                   COM = NULL,
                   TITLE = NULL) {
            if (is.null(con))
              con <- "spectrum.mgf"
            if (class(con) == "character") {
              if (file.exists(con)) {
                message("Overwriting ", con, "!")
                unlink(con)
              }
              con <- file(description = con,
                          open = "at",
                          blocking = TRUE)
              on.exit(close(con))
            }
            if (!inherits(con, "connection"))
              stop("'con' is not a proper connection!")
            if (is.null(COM))
              COM <- paste("COM=Spectrum exported by MSnbase on ",
                           date(), "\n", sep = "")
            ## write spectrum
            writeLines(COM, con = con)
            writeMgfContent(object, TITLE = TITLE, con = con)
            ## close(con)
          })

setMethod("writeMgfData",
          signature = signature("MSnExp"),
          function(object,
                   con = NULL,
                   COM = NULL,
                   verbose = TRUE) {
            if (is.null(con))
               con <- "experiment.mgf"           
            if (class(con) == "character") {
              if (file.exists(con)) {
                message("Overwriting ", con, "!")
                unlink(con)
              }
              con <- file(description = con,
                          open = "at",
                          blocking = TRUE)
              on.exit(close(con))
            }
            if (!inherits(con, "connection"))
              stop("'con' is not a proper connection!")
            if (is.null(COM))
              COM <- paste("COM=Experiment exported by MSnbase on ",
                           date(), "\n", sep = "")
            writeLines(COM, con = con)
            splist <- spectra(object)
            ## x <- sapply(splist,
            ##             writeMgfContent,
            ##             TITLE = NULL,
            ##             con = con)
            if (verbose)
                pb <- txtProgressBar(min = 0,
                                     max = length(splist),
                                     style = 3)
            for (i in 1:length(splist)) {
                if (verbose) setTxtProgressBar(pb, i)
                writeMgfContent(splist[[i]],
                                TITLE = NULL,
                                con = con)
            }
            if (verbose) close(pb)               
          })

writeMgfContent <- function(sp, TITLE = NULL, con) {
  buffer <- c("BEGIN IONS")
  buffer <- c(buffer,
              paste("SCANS=", sp@acquisitionNum, sep = ""))
  if (is.null(TITLE)) {
    TITLE <- paste("TITLE=msLevel ", sp@msLevel,
                   "; retentionTime ", sp@rt, 
                   "; scanNum ", sp@acquisitionNum, 
                   sep = "")
    if (length(sp@scanIndex) != 0)
      TITLE <- paste(TITLE,
                     "; scanIndex ", sp@scanIndex,
                     sep = "")
    if (sp@msLevel > 1)
      TITLE <- paste(TITLE,
                     "; precMz ", sp@precursorMz,
                     "; precCharge ", sp@precursorCharge,
                     sep = "")
  }
  buffer <- c(buffer, TITLE)
  buffer <- c(buffer,
              paste("RTINSECONDS=", sp@rt, sep = ""))
  buffer <- c(buffer,
              paste("PEPMASS=", sp@precursorMz, sep = ""))
  if ( !is.na(sp@precursorCharge) )
    buffer <- c(buffer,
                paste("CHARGE=", sp@precursorCharge, "+", sep = ""))
  dfr <- as(sp, "data.frame")
  pks <- apply(dfr,
               1,
               base::paste, collapse = " ")  
  buffer <- c(buffer, pks)
  buffer <- c(buffer,
              paste("END IONS"))
  writeLines(buffer, con = con)
}


# Code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
readMgfData <- function(file,
                        pdata = NULL,
                        centroided = TRUE,
                        smoothed = FALSE,
                        verbose = TRUE,
                        cache = 1) {
  if (verbose)
    cat("Scanning ", file,"...\n", sep="")
  mgf <- scan(file = file, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)
  ## From http://www.matrixscience.com/help/data_file_help.html#GEN
  ## Comment lines beginning with one of the symbols #;!/ can be included,
  ## but only outside of the BEGIN IONS and END IONS statements that delimit an MS/MS dataset.
  commentsymbols <- c("^#", "^;", "^!", "^/")
  cmts <- unlist(sapply(commentsymbols, grep, mgf))
  if (length(cmts) > 0)
    mgf <- mgf[-cmts]
  begin <- grep("BEGIN IONS", mgf)
  end <- grep("END IONS", mgf)
  if (verbose) {
    cnt <- 1
    pb <- txtProgressBar(min = 0, max = length(begin), style = 3)    
  }
  spectra <- vector("list",length=length(begin))
  fdata <- c()
  for (i in 1:length(begin)) {
    if (verbose) {
      setTxtProgressBar(pb, cnt)
      cnt <- cnt + 1
    }
    spec <- mgf[seq(begin[i], end[i])]
    spectra[[i]] <- mgfToSpectrum2(spec,centroided=centroided)
    ## adding the spectrum header as fData
    tmpdesc <- spec[grep("=",spec)]
    tmpdesc <- sapply(tmpdesc,strsplit,"=")
    desc <- unlist(lapply(tmpdesc,"[",2))
    names(desc) <- unlist(lapply(tmpdesc,"[",1))
    fdata <- rbind(fdata,desc)
  }
  if (verbose) 
    close(pb)
  nms <- paste("X",seq_len(length(spectra)),sep="")
  names(spectra) <- nms
  assaydata <- list2env(spectra)
  process <- new("MSnProcess",
                 processing = paste("Data loaded:",date()),
                 files = file,
                 smoothed = smoothed)
  if (is.null(pdata)) {
    pdata <- new("NAnnotatedDataFrame",
                 data = data.frame(sampleNames = file, fileNumbers = 1))
  }
  ## fdata <- new("AnnotatedDataFrame",
  ##              data=data.frame(spectrum=1:length(spectra),
  ##                row.names=nms))
  colnames(fdata) <- names(desc)
  rownames(fdata) <- nms
  fdata <- AnnotatedDataFrame(data = data.frame(fdata))
  fdata <- fdata[ls(assaydata), ] ## reorder features
  ## only levels 0 and 1 for mgf peak lists
  cache <- testCacheArg(cache, maxCache = 1)
  if (cache >= 1) {
    tmp <- new("MSnExp",
               assayData = assaydata,
               phenoData = pdata,
               featureData = fdata,
               processingData = process)
    newhd <- .header(tmp)
  } else {
    newhd <- NULL ## not used anyway
  }
  .cacheEnv <- setCacheEnv(list(assaydata = assaydata,
                                hd = newhd),
                           cache, lock = TRUE)
  toReturn <- new("MSnExp",
                  assayData = assaydata,
                  phenoData = pdata,
                  featureData = fdata,
                  processingData = process,
                  .cache = .cacheEnv)
  if (validObject(toReturn))
    return(toReturn)
}


mgfToSpectrum2 <- function(mgf, centroided) {
    mgf <- mgf[c(-1, -length(mgf))] ## remove "BEGIN IONS" and "END IONS"
    mgf <- sub("\t"," ", mgf) ## expecting mz and int to be separated by 1 space
    desc.idx <- grep("=", mgf)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]
    ## mz.i <- adply(.data=mgf[-desc.idx], .margins=1, .fun=function(i) unlist(strsplit(i, split=" ")))
    ## mz <- as.numeric(mz.i$V1)
    ## int <- as.numeric(mz.i$V2)
    mz.i <- sapply(spec,strsplit," ")
    mz <- as.numeric(unlist(lapply(mz.i,"[",1)))
    int <- as.numeric(unlist(lapply(mz.i,"[",2)))
    ## desc <- adply(.data=mgf[desc.idx], .margins=1, .fun=function(i) unlist(strsplit(i, split="=")))
    ## desc <- desc[,-1]
    tmpdesc <- sapply(desc, strsplit,"=")
    desc <- unlist(lapply(tmpdesc, "[", 2))
    names(desc) <- unlist(lapply(tmpdesc, "[", 1))
    ## sp <- new("Spectrum2",
    ##           rt=as.numeric(desc[desc[,1] == "RTINSECONDS",2]),
    ##           acquisitionNum=as.integer(desc[desc[,1] == "SCANS",2]),
    ##           precursorMz=as.numeric(desc[desc[,1] == "PEPMASS",2]),
    ##           #precursorIntensity="",  ##I don't see any definition in MGF.
    ##           precursorCharge=as.integer(sub("[+,-]","",desc[desc[,1] == "CHARGE",2])),
    ##           #scanIndex="",
    ##           #collisionEnergy="",
    ##           peaksCount=length(int),
    ##           mz=mz,
    ##           intensity=int,
    ##           fromFile=1L,
    ##           centroided=centroided)
    .parsePEPMASS <- function(desc, what = c("MZ", "INTENSITY")) {
      x <- strsplit(desc["PEPMASS"], " ")[[1]]
      what <- match.arg(what)
      if (length(x) == 1) x[2] <- 0
      ifelse(what == "MZ",
             as.numeric(x[1]),
             as.numeric(x[2]))
    }
    sp <- new("Spectrum2",
              rt = ifelse("RTINSECONDS" %in% names(desc),
                as.numeric(desc["RTINSECONDS"]), 0),
              scanIndex = ifelse("SCANS" %in% names(desc),
                as.integer(desc["SCANS"]), 0L),
              precursorMz = ifelse("PEPMASS" %in% names(desc),
                .parsePEPMASS(desc, "MZ"),
                0L),
              precursorIntensity = ifelse("PEPMASS" %in% names(desc),
                .parsePEPMASS(desc, "INTENSITY"),
                0L),
              precursorCharge = ifelse("CHARGE" %in% names(desc),
                as.integer(sub("[+,-]","",desc["CHARGE"])), 0L),
              peaksCount = length(int),
              mz = mz,
              intensity = int,
              fromFile = 1L,
              centroided = centroided)

    if(validObject(sp))
        return(sp)
}
