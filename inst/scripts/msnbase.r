#! /usr/local/bin/Rscript

#######################################################
## Laurent Gatto <lg390@cam.ac.uk>
##
## This script automates a simple MSnbase
## processing pipeline for isobaric reporter
## tagging data. It is provided as is, whithout
## any guarantee. Permission to use, copy, modify,
## and/or distribute this software for any purpose
## with or without fee is hereby granted, provided
## that the original software is cited.
## 
## ChangeLog:
## * 24-Dec-2012 Initial version v0.1.0
## * 08-Jan-2013 help arg/output v0.1.1
## * 20-Feb-2013 -keepna to retain features
##               with missing values v0.1.2
## * 26-Feb-2013 exact match of all 4/6/8 reporter 
##               ions and extract/return column indices
##               Also considering old TMT6 ispy tags.
#######################################################

version <- "0.1.3"

args <- commandArgs()
args <- args[-(1:match("--args", args))]

msg0 <-
  paste("-------------------------------------------------------------",
        "This free open-source software implements academic research.",
        "If you use it, please cite:  Gatto L. and Lilley KS., ",
        "MSnbase - an R/Bioconductor package for isobaric tagged mass",
        "spectrometry data visualisation, processing and quantitation.",
        "Bioinformatics, 28(2), 288-289, 2012 (PMID: 22113085).",
        "-------------------------------------------------------------\n\n",
        sep = "\n")

cat(msg0)

message("Analysis settings")
message(" - Running msnbase.r version ", version, ".")
message(" - Using MSnbase ", packageVersion("MSnbase"), ".")

if ("-h" %in%  args) {
  message(" -h: prints this help.")
  message(" -i: impurity matrix, generally defined as an csv file.")
  message("     See ?makeImpuritiesMatrix for details.")
  message(" -b: minimum reporter intensity filtering (default is 20000).")
  message(" -keepna: retain features with missing data.\n")
  message(" Other arguments are the files to be analysed.")
  message(" Non-existent file names are ignored.")
  q()
}


iarg <- match("-i", args)
if (!is.na(iarg)) {
  imp <- args[iarg + 1]
  args <- args[-(iarg:(iarg + 1))]
  if (file.exists(imp)) message(" - Using user-defined impurity matrix: ", imp, ".")
  else message(" - matrix ", imp, " not found.")
} else {
  message(" - No purity correction.")
  imp <- NULL
}

barg <- match("-b", args)
if (!is.na(barg)) {
  minint <- as.numeric(args[barg + 1])
  args <- args[-(barg:(barg + 1))]
  message(" - Using user-defined minimum intensity filtering: ", minint, ".")
} else {
  message(" - Minimum intensity 20000.")
  minint <- 20000
}

naarg <- match("-keepna", args)
if (!is.na(naarg)) {
  .na.rm <- FALSE
  args <- args[-naarg]
  message(" - Retaining features with missing data.")
} else {
  message(" - Removing features with missing data.")
  .na.rm <- TRUE
}


files <- unique(args)
if (length(files) == 0) {
  ## display help
  message("No arguments provided - done.\n")
  q()
}

found <- sapply(files, file.exists)
if (!all(found))
  warning('File(s) ', paste(files[!found], collapse = ", "), ' not found - ignoring.\n')

files <- files[found]
if (length(files) == 0) {
  message("No valid files provided - done.\n")
  q()
} else {
  message('Processing ', length(files), ' file(s): ',
          paste(files, collapse = ", "),"\n")
}

suppressPackageStartupMessages(library("MSnbase"))

whatIspyReporters <- function(f) {
  require("MSnbase")
  tmt6 <- c("126.128", "127.132", "128.14", "129.138",
            "130.142", "131.139")
  tmt6old <- c("126.128", "127.131", "128.135", "129.138",
               "130.141", "131.139")
  itraq4 <- c("114.1", "115.1", "116.1", "117.1")
  itraq8 <- c("113.1", "114.1", "115.1", "116.1",
              "117.1", "118.1", "119.1","120.08")
  hd <- scan(what = "character", file = f, nlines = 1, quiet = TRUE)
  if (all(tmt6 %in% hd)) {
    return(list(index = which(hd %in% tmt6),
                reporters = TMT6))
  }
  if (all(tmt6old %in% hd)) {
    return(list(index = which(hd %in% tmt6old),
                reporters = TMT6))
  }
  if (all(itraq4 %in% hd)) {
    return(list(index = which(hd %in% itraq4),
                reporters = iTRAQ4))
  }
  if (all(itraq8 %in% hd)) {
    return(list(index = which(hd %in% itraq8),
                reporters = iTRAQ8))
  }
  stop("ispy header does not match any expected patter.")
}

do <- function(f) {
  message("Processing ", f, "...", appendLF = FALSE)
  require("MSnbase")
  reps <- whatIspyReporters(f)
  ri <- reps$index
  if (is.null(imp)) {
    impurities <- matrix(0, length(ri), length(ri))
    diag(impurities) <- 1
  } else {
    impurities <- makeImpuritiesMatrix(filename = imp, edit = FALSE)
  } 
  ## processing
  x0 <- readIspyData(f, reporters = ri,
                     min.int = minint,
                     na.rm = .na.rm, 
                     verbose = FALSE)
  x2 <- purityCorrect(x0, impurities)
  if (.na.rm)
    x2 <- filterNA(x2, pNA = 0)
  xx <- combineFeatures(x2, groupBy = fData(x2)$ProteinAccession,
                        fun = "median", verbose = FALSE)
  ## additional meta-data
  tmp <- nQuants(x2, "ProteinAccession")
  stopifnot(all.equal(featureNames(xx), rownames(tmp)))
  if (.na.rm) {
    ## we expect having used the same number of features to perform
    ## quantitation if all feature with NAs were removed
    stopifnot(all(apply(tmp, 1, function(x) length(unique(x)) == 1)))
    nQuantSpectra <- tmp[, 1, drop = FALSE]
    colnames(nQuantSpectra) <- "nQuantSpectra"
  } else {
    nQuantSpectra <- tmp
    colnames(nQuantSpectra) <- paste("nQuantSpectra",
                                     colnames(nQuantSpectra),
                                     sep = ".")
  }
  fData(xx) <- cbind(fData(xx), nQuantSpectra)

  nQuantPeps <- tapply(fData(x2)$PeptideSequence,
                       fData(x2)$ProteinAccession,
                       function(x) length(unique(x)))
  
  stopifnot(all.equal(names(nQuantPeps), featureNames(xx)))

  fData(xx)$nQuantPeps <- nQuantPeps
  ## exporting
  cvcols <- grep("CV", fvarLabels(xx), value = TRUE)
  nqntsp <- grep("nQuantSpectra", fvarLabels(xx), value = TRUE)  
  cols <- c("ProteinAccession",                     
            "ProteinDescription",                   
            "ProteinQ.Value.uniquepeptidesonly.",
            "ProteinType1Error.uniquepeptidesonly.",
            "NumberOfUniquePeptides",
            "nQuantPeps",
            nqntsp,
            cvcols)
  f2 <- sub("\\.tsv", "_prot.tsv", f)
  write.exprs(xx, file = f2, fDataCols = cols)
  ## return
  message(" saved in ", f2, ".")
  invisible(xx)
}

res <- lapply(files, do)



