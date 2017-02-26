context("missing-data")

test_that(".preparePlotNAData", {
  m <- matrix(c(1, 1, 1, 1, NA,
                1, 1, 1, 1, 1,
                NA, 1, NA, NA, NA,
                NA, NA, NA, 1, 1), ncol=5, nrow=4, byrow=TRUE)
  r <- data.frame(x = 1:4,
                  proteins = c(1, 0.8, 0.4, 0.2),
                  data = c(1, 0.9, 1 - 4/15, 0.6))
  expect_equal(MSnbase:::.preparePlotNAData(m), r)
})
