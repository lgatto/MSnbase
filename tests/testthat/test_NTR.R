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
  expect_equal(MSnbase:::.referenceFractionValues(m,
                                                  factor(c("J", "G", "J", "G", "J"))),
                                                 c(2, 2, NA, 2, 3))
})

test_that("normToReference", {
  r <- matrix(c(3/5, 3/5, 1, 1/3,
                1, 5/6, 3/4, 1), nrow=2, byrow=TRUE,
              dimnames = list(c("J", "G"), paste0("F", 1:4)))
  r2 <- matrix(c(3/2, 3/2, 5/2, 1,
                 3/2, 5/2, 3, 3), nrow=2, byrow=TRUE,
               dimnames = list(c("J", "G"), paste0("F", 1:4)))
  r3 <- matrix(c(1, 1, 5/3, 1/2,
                 3/7, 5/9, 3/5, 6/9), nrow=2, byrow=TRUE,
               dimnames = list(c("J", "G"), paste0("F", 1:4)))

  expect_equal(MSnbase:::normToReference(m, group), r)
  expect_equal(MSnbase:::normToReference(m, group, reference=rep(1, 5)), r2)
  expect_equal(MSnbase:::normToReference(m, group, reference=1:5), r3)
  expect_equal(as.vector(MSnbase:::normToReference(m[1,, drop=FALSE], 1)),
               c(1/2, 1/2, 1, NA))
})

test_that("combineFeatures(..., method=\"NTR\")", {
  data(msnset)
  pa <- fData(msnset)$ProteinAccession

  ntr <- combineFeatures(msnset, groupBy=pa, method="NTR")
  ntrref <- combineFeatures(msnset, groupBy=pa, method="NTR", reference=
    exprs(msnset)[cbind(1:55, c(4, 4, 2, 3, 2, 3, 2, 4, 2, 2, 3, 4, 4, 4, 2,
                                1, 4, 3, 4, 1, 3, 3, 3, 2, 2, 3, 3, 2, 3, 3,
                                3, 4, 1, 4, 4, 3, 4, 3, 1, 4, 1, 1, 3, 3, 4,
                                3, 4, 4, 4, 4, 1, 3, 2, 3, 3))])
  expect_true(compareMSnSets(ntr, ntrref))

  ## sum and NTR are equal if the reference is
  ## 1/(number of peptides per protein) and contains no NA
  ref <- (1/table(pa)[pa])
  ntr <- combineFeatures(msnset, groupBy=pa, method="NTR", reference=ref)
  s <- combineFeatures(msnset, groupBy=pa, method="sum")
  ## exclude ENO (because of one NA the ref is not correct)
  ntr <- ntr[featureNames(ntr) != "ENO"]
  s <- s[featureNames(s) != "ENO"]
  expect_true(compareMSnSets(ntr, s))
})
