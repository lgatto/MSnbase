## Tests to run on a local set of (real) SRM data.
library(MSnbase)
library(testthat)

srm_fls <- dir("~/data/SRM/ss2017_11_28/", pattern = "mzML", full.names = TRUE)

srm <- readSRMData(srm_fls)
expect_equal(ncol(srm), length(srm_fls))
expect_equal(nrow(srm), 71)             # 71 analytes
expect_true(all(intensity(srm[1, 1]) != rtime(srm[1, 1])))

## Check the feature data.
## positive polarity
expect_true(all(fData(srm)$polarity == 1))
expect_true(all(polarity(srm) == 1))

pmz <- precursorMz(srm)
expect_true(is.matrix(pmz))
expect_equal(colnames(pmz), c("mzmin", "mzmax"))
expect_true(all(!is.na(pmz[, 1])))
expect_equal(pmz[, 1], fData(srm)$precursorIsolationWindowTargetMZ)

pmz <- productMz(srm)
expect_true(is.matrix(pmz))
expect_equal(colnames(pmz), c("mzmin", "mzmax"))
expect_true(all(!is.na(pmz[, 1])))
expect_equal(pmz[, 1], fData(srm)$productIsolationWindowTargetMZ)

## With phenoData
expect_error(readSRMData(srm_fls, pdata = "a"))
pheno <- data.frame(a = 1:3, b = 4:6)
expect_error(readSRMData(srm_fls, pdata = pheno))
pheno <- data.frame(sample = 1:length(srm_fls), batch = "A")
srm <- readSRMData(srm_fls, pdata = pheno)

expect_equal(colnames(pData(srm)), c("sample", "batch", "file"))
expect_equal(srm$sample, pheno$sample)
expect_equal(srm$batch, pheno$batch)

## Just visually inspect
plot(srm[1, ])
plot(srm[2, ])

