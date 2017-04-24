context("NTR")

m <- matrix(c(1, 1, 2, NA,
              2, 2, 3, 1,
              1, NA, NA, NA,
              2, 2, NA, 2,
              NA, 3, 3, 4), nrow=5, byrow=TRUE,
            dimnames = list(paste0("P", 1:5), paste0("F", 1:4)))

group <- c("J", "J", "G", "G", "G")

test_that(".referenceFraction", {
  expect_equal(MSnbase:::.referenceFraction(m, group), c(J=3, G=4))
})

test_that(".referenceFractionValues", {
  expect_error(MSnbase:::.referenceFractionValues(1:2, 1:2), "is.matrix\\(x\\)")
  expect_error(MSnbase:::.referenceFractionValues(m, 1:2),
               "nrow\\(x\\) == length\\(group\\)")
  expect_equal(MSnbase:::.referenceFractionValues(m, group), c(2, 3, NA, 2, 4))
  expect_equal(MSnbase:::.referenceFractionValues(m,
                                                  c("J", "G", "J", "G", "J")),
                                                 c(2, 2, NA, 2, 3))
})

test_that(".normToReference", {
  r <- matrix(c(3/5, 3/5, 1, 1/3,
                1, 5/6, 3/4, 1), nrow=2, byrow=TRUE,
              dimnames = list(c("J", "G"), paste0("F", 1:4)))
  r2 <- matrix(c(3/2, 3/2, 5/2, 1,
                 3/2, 5/2, 3, 3), nrow=2, byrow=TRUE,
               dimnames = list(c("J", "G"), paste0("F", 1:4)))
  r3 <- matrix(c(1, 1, 5/3, 1/2,
                 3/7, 5/9, 3/5, 6/9), nrow=2, byrow=TRUE,
               dimnames = list(c("J", "G"), paste0("F", 1:4)))

  expect_equal(MSnbase:::.normToReference(m, group), r/c(38/15, 43/12))
  expect_equal(MSnbase:::.normToReference(m, group, norm=FALSE), r)
  expect_equal(MSnbase:::.normToReference(m, group, reference=rep(1, 5),
                                          norm=FALSE), r2)
  expect_equal(MSnbase:::.normToReference(m, group, reference=1:5,
                                          norm=FALSE), r3)
  expect_equal(as.vector(MSnbase:::.normToReference(m[1,, drop=FALSE], 1,
                                                    norm=FALSE)),
               c(1/2, 1/2, 1, NA))
  expect_equal(as.vector(MSnbase:::.normToReference(m[1,, drop=FALSE], 1)),
               c(1/4, 1/4, 1/2, NA))
})
