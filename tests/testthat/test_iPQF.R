context("iPQF feature aggregation")

test_that("data is valid", {
              ## test input
              data("msnset2")
              expect_true(validObject(msnset2))
              ## test output
              fl <- list.files(system.file("extdata", package = "MSnbase"),
                               full.names=TRUE, pattern = "ipqfres.rda")
              nms <- load(fl)
              val <- sapply(nms, function(x) validObject(get(x)))
              expect_true(all(val))
          })

test_that("iPQF unit test", {
              data("msnset2")
              fl <- list.files(system.file("extdata", package = "MSnbase"),
                               full.names=TRUE, pattern = "ipqfres.rda")
              load(fl)
              ##
              res1 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = FALSE,
                                      ratio.calc = "sum",
                                      method.combine = FALSE)
              res1@processingData <- processingData(ipqf1)
              expect_true(all.equal(res1, ipqf1))

              dflt <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF")
              dflt@processingData <- processingData(ipqf1)
              expect_true(all.equal(res1, dflt))
              

              res2 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = FALSE,
                                      ratio.calc = "sum",
                                      method.combine = TRUE)
              res2@processingData <- processingData(ipqf2)
              expect_true(all.equal(res2, ipqf2))

              res3 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = TRUE,
                                      ratio.calc = "sum",
                                      method.combine = FALSE)
              res3@processingData <- processingData(ipqf3)
              expect_true(all.equal(res3, ipqf3))

              res4 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = TRUE,
                                      ratio.calc = "none",
                                      method.combine = FALSE)
              res4@processingData <- processingData(ipqf4)
              expect_true(all.equal(res4, ipqf4))

              res5 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = TRUE,
                                      ratio.calc = "X114.ions",
                                      method.combine = FALSE)
              res5@processingData <- processingData(ipqf5)
              expect_true(all.equal(res5, ipqf5))
          })
