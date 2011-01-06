context("MSnExp class")

test_that("readMzXMLData and MSnExp instance", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,verbose=FALSE)
  expect_that(class(aa)=="MSnExp",is_true())
  expect_that(length(aa),equals(98))
  ## MS levels
  expect_that(length(msLevel(aa)),equals(98))
  expect_that(unique(msLevel(aa)),equals(2))
  expect_that(length(ms1scan(aa)),equals(98))
  expect_that(length(unique(ms1scan(aa))),equals(7))
  ## Precursor MZ
  expect_that(length(precursorMz(aa)),equals(98))
  expect_that(precursorMz(aa)[1],is_a("numeric"))
  expect_that(length(unique(precursorMz(aa))),equals(15))
  expect_that(range(precursorMz(aa)),
              equals(c(424.76651001,1007.04766850)))
  expect_that(as.numeric(sort(precursorMz(aa))[1]),equals(424.76651001)) ## [*]
  ## Retention time
  expect_that(length(rtime(aa)),equals(98))
  expect_that(rtime(aa)[1],is_a("numeric"))
  expect_that(range(rtime(aa)),
              equals(c(1982.92,6179.02)))
  expect_that(as.numeric(sort(rtime(aa))[1]),equals(1982.92)) ## [*]
  ## [*] using as.numeric because rtime and precursorMz return named numerics
})
