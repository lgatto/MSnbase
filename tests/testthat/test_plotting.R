context("Plotting functions")

test_that("image,MSnSet-method", {
              data(naset)
              par(mfrow = c(1, 2))
              expect_is(image(naset), "gg")
              expect_null(image2(exprs(naset)))
              dev.off()
          })
