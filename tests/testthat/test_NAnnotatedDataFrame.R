context("NAnnotatedDataFrame")

test_that("NAnnotatedDataFrame validity", {
    na <- new("NAnnotatedDataFrame")
    expect_true(validObject(na))
    dnms <- c(rowNames = 0, columnNames = 0, multiNames = 1)
    expect_true(all.equal(dim(na), dnms))
    pd <- data.frame(pvarA = "A", pvarB = "B")
    na <- new("NAnnotatedDataFrame", data = pd) 
    expect_true(validObject(na))
    dnms <- c(rowNames = 1, columnNames = 2, multiNames = 1)
    expect_true(all.equal(dim(na), dnms))
    na <- new("NAnnotatedDataFrame",
              multiplex = 4, multiLabels = reporterNames(iTRAQ4),
              data = pd)
    expect_true(validObject(na))
    dnms <- c(rowNames = 1, columnNames = 2, multiNames = 4)
    expect_true(all.equal(dim(na), dnms))

    expect_true(multiplex(na) == 4)
    expect_true(all.equal(multiLabels(na), reporterNames(iTRAQ4)))
})
