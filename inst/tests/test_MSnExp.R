context("MSnExp class")

test_that("readMzXMLData and MSnExp instance", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMzXMLData(file,verbose=FALSE)
  expect_that(length(aa),equals(98))
  expect_that(length(precursorMz(aa)),equals(98))
  expect_that(length(rtime(aa)),equals(98))
  expect_that(rtime(aa)[1],is_a("numeric"))
  expect_that(as.numeric(sort(rtime(aa))[1]),equals(1982.92))
  expect_that(as.numeric(sort(precursorMz(aa))[1]),equals(424.76651001))
  expect_that(length(unique(precursorMz(aa))),equals(15))
  expect_that(class(aa)=="MSnExp",is_true())
})
