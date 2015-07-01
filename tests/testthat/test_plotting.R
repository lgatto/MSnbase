context("Plotting functions")

test_that("image,MSnSet-method", {
              data(naset)
              ## only difference between these two plots is that
              ## the former uses Samples/Features as xlab/ylab 
              ## and the latter Columns/Samples
              par(mfrow = c(1, 2))
              expect_null(image(naset))
              expect_null(image2(exprs(naset)))
              dev.off()
          })
