test_that("fileNames accessor MSnExp", {
    data(itraqdata, package = "MSnbase")
    expect_identical(fileNames(itraqdata), "dummyiTRAQ.mzXML")
    expect_identical(fileNames(itraqdata),
                     fileNames(processingData(itraqdata)))
})

test_that("fileNames accessor MzTab", {
    f <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/MTBLS2.mztab"
    expect_identical(fileName(MzTab(f)), f)
    expect_identical(fileName(MzTab(f)), fileNames(MzTab(f)))
})

test_that("fileNames accessor MSmap", {
    f <- dir(system.file("threonine", package = "msdata"),
             full.names = TRUE)
    ms <- openMSfile(f)
    map <- MSmap(ms, lowMz = 200, highMz = 500, resMz = 1)
    expect_identical(f, fileName(map))
    expect_identical(fileName(map), fileNames(map))
})

test_that("fileNames accessor MSnSet", {
    data(msnset, package = "MSnbase")
    expect_identical(fileNames(msnset), "dummyiTRAQ.mzXML")
})
