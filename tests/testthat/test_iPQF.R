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
                               full.names = TRUE, pattern = "ipqfres.rda")
              load(fl) ## ipqf1 to 5
              ##
              res1 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = FALSE,
                                      ratio.calc = "sum",
                                      method.combine = FALSE)
              res1@processingData <- processingData(ipqf1)
              res1@experimentData <- experimentData(ipqf1)
              res1@.__classVersion__ <- ipqf1@.__classVersion__
              expect_true(all.equal(res1, ipqf1))

              dflt <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF")
              dflt@processingData <- processingData(ipqf1)
              dflt@experimentData <- experimentData(ipqf1)
              res1@.__classVersion__ <- dflt@.__classVersion__
              expect_true(all.equal(res1, dflt))
              

              res2 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = FALSE,
                                      ratio.calc = "sum",
                                      method.combine = TRUE)
              res2@processingData <- processingData(ipqf2)
              res2@experimentData <- experimentData(ipqf2)
              res2@.__classVersion__ <- ipqf2@.__classVersion__
              expect_true(all.equal(res2, ipqf2))

              res3 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = TRUE,
                                      ratio.calc = "sum",
                                      method.combine = FALSE)
              res3@processingData <- processingData(ipqf3)
              res3@experimentData <- experimentData(ipqf3)
              res3@.__classVersion__ <- ipqf3@.__classVersion__
              expect_true(all.equal(res3, ipqf3))

              res4 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = TRUE,
                                      ratio.calc = "none",
                                      method.combine = FALSE)
              res4@processingData <- processingData(ipqf4)
              res4@experimentData <- experimentData(ipqf4)
              res4@.__classVersion__ <- ipqf4@.__classVersion__
              expect_true(all.equal(res4, ipqf4))

              res5 <- combineFeatures(msnset2,
                                      groupBy = fData(msnset2)$accession,
                                      redundancy.handler = "unique",
                                      fun = "iPQF",
                                      low.support.filter = TRUE,
                                      ratio.calc = "X114.ions",
                                      method.combine = FALSE)
              res5@processingData <- processingData(ipqf5)
              res5@experimentData <- experimentData(ipqf5)
              res5@.__classVersion__ <- ipqf5@.__classVersion__
              expect_true(all.equal(res5, ipqf5))
          })
