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

# Based on the code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# Modified by Sebastian Gibb <mail@sebastiangibb.de>
readMgfData <- function(file,
                        pdata = NULL,
                        centroided = TRUE,
                        smoothed = FALSE,
                        verbose = TRUE,
                        cache = 1) {
  if (verbose)
    cat("Scanning", file, "...\n")

  mgf <- scan(file = file, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)

  ## From http://www.matrixscience.com/help/data_file_help.html#GEN
  ## Comment lines beginning with one of the symbols #;!/ can be included,
  ## but only outside of the BEGIN IONS and END IONS statements that delimit an MS/MS dataset.
  cmts <- grep("^[#;!/]", mgf)
  if (length(cmts))
    mgf <- mgf[-cmts]

  begin <- grep("BEGIN IONS", mgf) + 1L
  end <- grep("END IONS", mgf) - 1L

  n <- length(begin)

  if (verbose) {
    cnt <- 1L
    pb <- txtProgressBar(min=0L, max=n, style=3L)
  }

  spectra <- vector("list", length=n)
  fdata <- vector("list", length=n)

  for (i in seq(along=spectra)) {
    if (verbose) {
      setTxtProgressBar(pb, cnt)
      cnt <- cnt + 1L
    }
    specInfo <- extractMgfSpectrum2Info(mgf[begin[i]:end[i]],
                                        centroided=centroided)
    spectra[[i]] <- specInfo$spectrum
    fdata[[i]] <- specInfo$fdata
  }
  if (verbose)
    close(pb)

  fdata <- do.call(rbind, fdata)

  names(spectra) <- paste0("X", seq_along(spectra))
  assaydata <- list2env(spectra)
  process <- new("MSnProcess",
                 processing = paste("Data loaded:", date()),
                 files = file,
                 smoothed = smoothed)
  if (is.null(pdata)) {
    pdata <- new("NAnnotatedDataFrame",
                 data = data.frame(sampleNames = file, fileNumbers = 1))
  }
  rownames(fdata) <- names(spectra)
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

extractMgfSpectrum2Info <- function(mgf, centroided) {
  ## grep description
  desc.idx <- grep("=", mgf)
  desc <- mgf[desc.idx]
  spec <- mgf[-desc.idx]

  ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
  mode(ms) <- "double"

  if (!length(ms))
    ms <- matrix(numeric(), ncol=2L)

  desc <- do.call(rbind, strsplit(desc, "=", fixed=TRUE))
  desc <- setNames(desc[, 2L], desc[, 1L])
  fdata <- desc

  desc[c("PEPMASSMZ", "PEPMASSINT")] <- strsplit(desc["PEPMASS"], "[[:space:]]+")[[1L]][1:2]

  ## select only values of interest and convert to numeric
  desc["CHARGE"] <- sub("[+-]", "", desc["CHARGE"])
  voi <- c("RTINSECONDS", "CHARGE", "SCANS", "PEPMASSMZ", "PEPMASSINT")
  desc <- setNames(as.numeric(desc[voi]), voi)
  desc[is.na(desc[voi])] <- 0L

  sp <- new("Spectrum2",
            rt = desc["RTINSECONDS"],
            scanIndex = as.integer(desc["SCANS"]),
            precursorMz = desc["PEPMASSMZ"],
            precursorIntensity = desc["PEPMASSINT"],
            precursorCharge = as.integer(desc["CHARGE"]),
            peaksCount = nrow(ms),
            mz = ms[, 1L],
            intensity = ms[, 2L],
            fromFile = 1L,
            centroided = centroided)

  if(validObject(sp))
    return(list(spectrum=sp, fdata=fdata))
}

