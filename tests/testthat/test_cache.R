context("cache testing")

test_that("Caching MS2 data", {
  file <- dir(system.file(package = "MSnbase", dir = "extdata"),
              full.names = TRUE,
              pattern = "mzXML$")
  c1 <- readMSData(file, verbose = FALSE, cache = 1)
  expect_true(MSnbase:::checkHeader(c1))
  c0 <- readMSData(file, verbose = FALSE, cache = 0)
  expect_true(MSnbase:::checkHeader(c0))
  expect_true(MSnbase:::getCacheEnv(c1)$level == 1)
  expect_true(MSnbase:::getCacheEnv(c0)$level == 0)
  show(c1)
  show(c0)
  expect_true(all.equal(c0, c1, check.attributes = FALSE))
})
