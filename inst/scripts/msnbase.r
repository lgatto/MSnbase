#! /usr/local/bin/Rscript

###################################################
## Laurent Gatto, 2012-2013, <lg390@cam.ac.uk>
## This script automates a simple MSnbase
## processing pipeline for isobaric reporter
## tagging data. It is provided as is, whithout
## any guarantee. Permission to use, copy, modify,
## and/or distribute this software for any purpose
## with or without fee is hereby granted, provided
## that the original software is cited.
## 
## ChangeLog:
## * 24-Dec-2012 Initial version 0.1.0
###################################################

version <- "0.1.0"

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

iarg <- match("-i", args)
if (!is.na(iarg)) {
  imp <- args[iarg + 1]
  args <- args[-(iarg:(iarg + 1))]
  message(" - Using user-defined impurity matrix: ", imp, ".")
} else {
  message(" - No purity correction.")
  imp <- NULL
}

barg <- match("-b", args)
if (!is.na(barg)) {
  minint <- args[barg + 1]
  args <- args[-(barg:(barg + 1))]
  message(" - Using user-defined minimum intensity filtering: ", minint, ".")
} else {
  message(" - Minimum intensity 20000.")
  minint <- 20000
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
  tmt6 <- "126.128"
  itraq4 <- c("114.1", "115.1", "116.1", "117.1" ,"145.1")
  itraq8 <- c("113.1", "114.1", "115.1", "116.1",
              "117.1", "118.1", "119.1","120.08" ,"121.1")
  hd <- scan(what = "character", file = f, nlines = 1, quiet = TRUE)
  if (all(tmt6 %in% hd)) return(TMT6)
  if (all(itraq4 %in% hd)) return(iTRAQ4)
  if (all(itraq8 %in% hd)) return(iTRAQ8)
  stop("ispy header does not match any expected patters.")
}

do <- function(f) {
  message("Processing ", f, "...", appendLF = FALSE)
  require("MSnbase")
  ri <- whatIspyReporters(f)
  r <- 19:(19 + length(ri) - 1)    
  if (is.null(imp)) {
    impurities <- matrix(0, length(ri), length(ri))
    diag(impurities) <- 1
  } else if (file.exists(imp)) {
    impurities <- makeImpuritiesMatrix(filename = imp, edit = FALSE)
  } else {
    impurities <- makeImpuritiesMatrix(as.numeric(imp), edit = FALSE)
  }  
  ## processing
  x0 <- readIspyData(f, reporters = r, min.int = minint, verbose = FALSE)
  x2 <- filterNA(purityCorrect(x0, impurities), pNA = 0)
  xx <- combineFeatures(x2, groupBy = fData(x2)$ProteinAccession,
                        fun = "median", verbose = FALSE)
  ## additional meta-data
  tmp <- nQuants(x2, "ProteinAccession")
  stopifnot(all(apply(tmp, 1, function(x) length(unique(x)) == 1)))
  stopifnot(all.equal(featureNames(xx), rownames(tmp)))
  nQuantSpectra <- as.numeric(tmp[, 1])
  nQuantPeps <- tapply(fData(x2)$PeptideSequence,
                   fData(x2)$ProteinAccession,
                   function(x) length(unique(x)))
  all.equal(names(nQuantPeps), featureNames(xx))
  fData(xx)$nQuantSpectra <- nQuantSpectra
  fData(xx)$nQuantPeps <- nQuantPeps
  ## exporting
  cvcols <- grep("CV", fvarLabels(xx), value = TRUE)
  cols <- c("ProteinAccession",                     
            "ProteinDescription",                   
            "ProteinQ.Value.uniquepeptidesonly.",
            "ProteinType1Error.uniquepeptidesonly.",
            "NumberOfUniquePeptides",
            "nQuantPeps",
            "nQuantSpectra",
            cvcols)
  f2 <- sub("\\.tsv", "_prot.tsv", f)
  write.exprs(xx, file = f2, fDataCols = cols)
  ## return
  message(" saved in ", f2, ".")
  invisible(xx)
}

res <- lapply(files, do)



