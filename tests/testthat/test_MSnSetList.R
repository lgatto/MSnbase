context("MSnSetLists")

test_that("MSnSetList initialisation", {
    expect_true(validObject(MSnSetList()))
    data(msnset)
    n <- 40L
    expect_true(validObject(xl <- split(msnset, "ProteinAccession")))
    expect_identical(nrow(fData(xl)), n)
    expect_identical(ncol(fData(xl)), 0L)
    fData(xl)$foo <- 1:40
    expect_identical(ncol(fData(xl)), 1L)
    expect_true(validObject(lapply(xl, normalise, method = "sum")))
    i <- c(2, 35, 12)
    j <- names(xl)[i]
    k <- rep(FALSE, n)
    k[i] <- TRUE
    expect_true(validObject(xli <- xl[i]))
    expect_true(validObject(xlj <- xl[j]))
    expect_true(validObject(xlk <- xl[k]))
    expect_identical(xli, xlj)
    expect_identical(xli, xlk[j])

    rl <- lapply(xli, nrow)
    rs <- sapply(xli, nrow)
    expect_identical(rs, unlist(rl))
    rs <- sapply(xli, dim)
    expect_identical(dim(rs), c(2L, 3L))
})
