test_that("BackendHdf5 validators work", {
    expect_true(is.null(.valid.BackendHdf5.md5sum(5, 5)))
    expect_true(is.character(.valid.BackendHdf5.md5sum(3, 1:2)))
    expect_true(is.null(.valid.BackendHdf5.hdf5file(5, 5)))
    expect_true(is.character(.valid.BackendHdf5.hdf5file(3, 1:2)))
    expect_true(is.character(.valid.BackendHdf5.hdf5file(c(1, 1), c(1, 2))))
    tmpf <- tempfile()
    write(4, tmpf)
    expect_true(is.null(.valid.BackendHdf5.hdf5file.exist(tmpf)))
    expect_true(is.character(.valid.BackendHdf5.hdf5file.exist("4")))
})

test_that("BackendHdf5 works", {
    tst <- BackendHdf5()
    expect_true(is(tst, "BackendHdf5"))
    expect_true(validObject(tst))
    tst@md5sum <- "r"
    expect_error(validObject(tst))
})

test_that("backendInitialize and backendImportData,BackendHdf5 work", {
    tst <- BackendHdf5()
    dr <- tempdir()
    res <- MSnbase:::backendInitialize(tst, files = sf, path = dr)
    expect_equal(length(res@hdf5file), 2)
    expect_equal(length(res@md5sum), 2)
    expect_true(all(file.exists(res@hdf5file)))
    expect_error(MSnbase:::backendInitialize(tst, files = sf, path = dr))
    expect_true(validObject(res))
    expect_true(all(res@md5sum == c("", "")))
    show(res)

    ## backendImportData
    spdata <- sciex_mzr@spectraData
    res <- MSnbase:::backendImportData(res, spdata)
    cont <- rhdf5::h5ls(res@hdf5file[1])
    expect_equal(sum(spdata$fileIdx == 1) + 3, nrow(cont))
    cont <- rhdf5::h5ls(res@hdf5file[2])
    expect_equal(sum(spdata$fileIdx == 2) + 3, nrow(cont))
})

test_that(".serialize_msfile_to_hdf5 works", {
    h5 <- tempfile()
    md5 <- MSnbase:::.serialize_msfile_to_hdf5(fileNames(sciex)[1], h5)
    cont <- rhdf5::h5ls(h5)
    expect_equal(nrow(cont), sum(fromFile(sciex) == 1) + 3)
})
