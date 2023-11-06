# Based on the code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# Modified by Sebastian Gibb <mail@sebastiangibb.de>
setMethod("writeMgfData",
          signature = signature("Spectrum"),
          function(object,
                   con = "spectrum.mgf",
                   COM = NULL,
                   TITLE = NULL) {
            writeMgfDataFile(list(object), con = con, COM = COM, TITLE = TITLE,
                             verbose = isMSnbaseVerbose())
          })

setMethod("writeMgfData",
          signature = signature("MSnExp"),
          function(object,
                   con = "experiment.mgf",
                   COM = NULL,
                   verbose = isMSnbaseVerbose()) {
            writeMgfDataFile(spectra(object), con = con, COM = COM,
                             verbose = verbose)
          })

#' @param addFields `data.frame` or `matrix` with optional additional
#'     fields to be added to each spectrum (each column one field).
#'
#' @noRd
writeMgfDataFile <- function(splist, con, COM = NULL, TITLE = NULL,
                             verbose = isMSnbaseVerbose(), addFields = NULL) {
  if (class(con) == "character" && file.exists(con)) {
    message("Overwriting ", con, "!")
    unlink(con)
  }

  if (length(addFields)) {
      if (length(dim(addFields)) != 2)
          stop("'addFields' has to be a matrix or data.frame.")
      if (!is.matrix(addFields))
          addFields <- do.call(cbind, lapply(addFields, as.character))
      if (is.null(colnames(addFields)))
          stop("Column names required on 'addFields'.")
      if (nrow(addFields) != length(splist))
          stop("nrow of 'addFields' has to match length of 'splist'")
  }

  if (class(con)[1] == "character") {
    con <- file(description = con, open = "at")
    on.exit(close(con))
  }

  if (is.null(COM)) {
    COM <- paste0(ifelse(length(splist) <= 1, "Spectrum", "Experiment"),
                  "exported by MSnbase on ", date())
  }
  cat(paste0("COM=",COM), file = con, sep = "")


  verbose <- verbose & length(splist) > 1

  if (verbose)
    pb <- txtProgressBar(min = 0, max = length(splist), style = 3)

  for (i in seq(along=splist)) {
    if (verbose)
      setTxtProgressBar(pb, i)

    writeMgfContent(splist[[i]], TITLE = TITLE, con = con,
                        addFields = addFields[i, ])
  }
  if (verbose)
    close(pb)
}

writeMgfContent <- function(sp, TITLE = NULL, con, addFields = NULL) {
    .cat <- function(..., file = con, sep = "", append = TRUE) {
        cat(..., file = file, sep = sep, append = append)
    }

    .cat("\nBEGIN IONS\n",
       "SCANS=", acquisitionNum(sp))

    if (is.null(TITLE)) {
        .cat("\nTITLE=msLevel ", msLevel(sp),
             "; retentionTime ", rtime(sp),
             "; scanNum ", acquisitionNum(sp))

        if (length(scanIndex(sp))) {
            .cat("; scanIndex ", scanIndex(sp))
        }

        if (msLevel(sp) > 1) {
            .cat("; precMz ", precursorMz(sp),
                 "; precCharge ", precursorCharge(sp))
        }
    } else {
        .cat("\nTITLE=", TITLE)
    }

    if (msLevel(sp) > 1) {
        .cat("\nRTINSECONDS=", rtime(sp), "\nPEPMASS=", precursorMz(sp))
        if (length(precursorCharge(sp)) && !is.na(precursorCharge(sp))) {
            .cat("\nCHARGE=", precursorCharge(sp), "+")
        }
    } else .cat("\nRTINSECONDS=", rtime(sp))

    if (length(addFields) && !is.null(names(addFields)))
        .cat("\n", paste(toupper(names(addFields)),
                         addFields, sep = "=", collapse = "\n"))

    .cat("\n", paste(mz(sp), intensity(sp), collapse = "\n"))
    .cat("\nEND IONS\n")
}

                                        # Based on the code contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# Modified by Sebastian Gibb <mail@sebastiangibb.de>
readMgfData <- function(filename,
                        pdata = NULL,
                        centroided = TRUE,
                        smoothed = FALSE,
                        verbose = isMSnbaseVerbose(),
                        cache = 1) {
  if (verbose)
    cat("Scanning", filename, "...\n")

  mgf <- scan(file = filename, what = "",
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
    pb <- txtProgressBar(min = 0L, max = n, style = 3L)
  }

  spectra <- vector("list", length = n)
  fdata <- vector("list", length = n)

  for (i in seq(along = spectra)) {
    if (verbose) {
      setTxtProgressBar(pb, cnt)
      cnt <- cnt + 1L
    }
    specInfo <- extractMgfSpectrum2Info(mgf[begin[i]:end[i]],
                                        centroided = centroided)
    spectra[[i]] <- specInfo$spectrum
    fdata[[i]] <- specInfo$fdata
  }
  if (verbose)
    close(pb)

  ## This checks that the headers are all the same
  fdn <- sapply(fdata, names, simplify = FALSE)
  fdn1 <- fdn[[1]]
  if (!all(sapply(fdn, function(x) identical(x, fdn1)))) {
      stop("Ion headers identical.")
  }

  fdata <- do.call(rbind, fdata)

  names(spectra) <- paste0("X", seq_along(spectra))
  assaydata <- list2env(spectra)
  process <- new("MSnProcess",
                 processing = paste("Data loaded:", date()),
                 files = filename,
                 smoothed = smoothed)
  if (is.null(pdata)) {
    pdata <- new("AnnotatedDataFrame",
                 data = data.frame(sampleNames = filename,
                                   fileNumbers = 1))
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
        ms <- matrix(numeric(), ncol = 2L)

    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)), substring(desc, 1L, r - 1L))
    fdata <- desc

    desc[c("PEPMASSMZ", "PEPMASSINT")] <-
        strsplit(desc["PEPMASS"], "[[:space:]]+")[[1L]][1:2]

    ## select only values of interest and convert to numeric
    desc["CHARGE"] <- sub("[+-]", "", desc["CHARGE"])
    voi <- c("RTINSECONDS", "CHARGE", "SCANS", "PEPMASSMZ", "PEPMASSINT")
    desc <- setNames(as.numeric(desc[voi]), voi)
    desc[is.na(desc[voi])] <- 0L
    sp <- Spectrum2_mz_sorted(rt = unname(desc["RTINSECONDS"]),
                              scanIndex = unname(as.integer(desc["SCANS"])),
                              precursorMz = unname(desc["PEPMASSMZ"]),
                              precursorIntensity = unname(desc["PEPMASSINT"]),
                              precursorCharge = unname(as.integer(desc["CHARGE"])),
                              mz = ms[, 1L],
                              intensity = ms[, 2L],
                              fromFile = 1L,
                              centroided = centroided)
    return(list(spectrum = sp, fdata = fdata))
}
