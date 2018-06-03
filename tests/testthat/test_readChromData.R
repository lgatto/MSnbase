test_that("readSRMData works", {
    library(msdata)
    fl <- proteomics(full.names = TRUE, pattern = "MRM")
    files <- c(fl, fl)

    ## Errors
    expect_error(readSRMData())
    expect_error(readSRMData(files, pdata = data.frame(files = "a")))
    expect_error(readSRMData(files, pdata = "a"))

    ## Read the data
    expect_warning(mrm <- readSRMData(files))
    expect_true(colnames(pData(mrm)) == "file")
    expect_equal(mrm$file, pData(mrm)$file)

    ## Test that the precursorMz from each chromatogram matches the one from
    ## the feature data.
    precMz <- do.call(rbind, lapply(mrm[, 1], precursorMz))
    colnames(precMz) <- c("mzmin", "mzmax")
    expect_equal(precursorMz(mrm), precMz)
    expect_equal(fData(mrm)$precursorIsolationWindowTargetMZ, precMz[, 1])
    ## Same for the productMz
    prodMz <- do.call(rbind, lapply(mrm[, 1], productMz))
    colnames(prodMz) <- c("mzmin", "mzmax")
    expect_equal(productMz(mrm), prodMz)
    expect_equal(fData(mrm)$productIsolationWindowTargetMZ, prodMz[, 1])
    
    ## Read with pheno data
    expect_warning(mrm <- readSRMData(files,
                                      pdata = data.frame(sample = 1:2)))
    expect_equal(mrm$sample, 1:2)
})

test_that(".combine_data.frame works", {
    A <- data.frame(a = c(1, 1, 2, 3, 4, 5, 6), b = c(2, 2, 3, 4, 5, 5, 6))
    B <- data.frame(a = c(1, 3, 3, 3, 3, 4, 5), b = c(2, 4, 4, 4, 5, 5, 5))
    C <- data.frame(a = c(3, 3, 4, 4), b = c(4, 5, 5, 5))
    
    exp <- data.frame(a = c(1, 1, 2, 3, 3, 3, 3, 4, 4, 5, 6),
                      b = c(2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 6))
    expect_equal(exp, .combine_data.frame(list(A, B, C)))
    
    expect_error(.combine_data.frame())
    expect_error(.combine_data.frame(list(A, B, C), cols = c("z")))

    ## When I do have factors?
    D <- data.frame(a = c("z", "b", "a", "g", "g"),
                    b = c(1, 2, 3, 4, 4), stringsAsFactors = FALSE)
    E <- data.frame(a = c("g", "a", "d", "g"),
                    b = c(4, 3, 1, 4), stringsAsFactors = FALSE)
    exp <- data.frame(a = c("a", "b", "d", "g", "g", "z"),
                      b = c(3, 2, 1, 4, 4, 1), stringsAsFactors = FALSE)
    expect_equal(.combine_data.frame(list(D, E)), exp)
    expect_equal(.combine_data.frame(list(E, D)), exp)
})

test_that(".polarity_char works", {
    expect_equal(.polarity_char(1), "+")
    expect_equal(.polarity_char(0), "-")
    expect_equal(.polarity_char(-1), NA)
    expect_equal(.polarity_char(c(1, 1, 0, 1, -1)), c("+", "+", "-", "+", NA))
    expect_error(.polarity_char(c(2, 1, 0)))
})
