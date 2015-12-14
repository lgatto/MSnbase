context("topN method")

test_that("topN for matrix", {
  M <- matrix(c(1.0, 2.1, 10.9, # group 1
                2.1, 1.1, 2.3,
                3.9, 4.1, 4.1,
                1.0, 1.1, 0.9, # group 2
                1.0, 1.1, 0.9,
                2.1, 2.0, 2.3,
                3.9, 4.1, 4.1,
                1.0, 1.1, 0.9, # group 3
                2.0, 2.1, 2.9),
              9, 3, byrow =TRUE)
  groups <- factor(c(rep("XX", 3),
                     rep("BB", 4),
                     rep("CC", 2)))
  rownames(M) <-
    paste("r", 1:9, groups, sep = "")
  Mtop2sum <- M[c(7, 6,
                  9, 8,
                  1, 3), ]
  Mtop2max <- M[c(7, 6,
                  9, 8,
                  1, 3), ]
  Mtop2min <- M[c(7, 6,
                  9, 8,
                  3, 2), ]
  M2 <- topN(M, groupBy = groups, n = 1)
  expect_equal(dim(M2), c(3, ncol(M)))
  M2 <- topN(M, groupBy = groups, n = 2)
  expect_equal(dim(M2), c(6, ncol(M)))
  expect_equal(Mtop2sum, M2)
  expect_equal(Mtop2max, topN(M, groupBy = groups, n = 2, fun = max))
  expect_equal(Mtop2min, topN(M, groupBy = groups, n = 2, fun = min))
})

test_that("topN for matrix and NAs", {
  M <- matrix(c(10,  NA, 10.9, # group 1
                2.1, 1.1, 2.3,
                3.9, 4.1, 4.1,
                1.1, 1.1, 1.1 ), ## group 2
              4, 3, byrow =TRUE)
  groups <- factor(c(rep("AA", 3),
                     rep("BB", 1)))
  rownames(M) <-
    paste("r", 1:4, groups, sep = "")
  r1 <- rownames(topN(M, groupBy = groups, n = 1, fun = sum))
  expect_equal(r1, c("r3AA", "r4BB"))
  r2 <- rownames(topN(M, groupBy = groups, n = 1, fun = sum, na.rm = TRUE))
  expect_equal(r2, c("r1AA", "r4BB"))
})

test_that("topN for MSnSet", {
    xx <- quantify(itraqdata, reporters = iTRAQ4,
                   method = "max", verbose=FALSE,
                   BPPARAM = SerialParam())
    for (.n in c(1:4)) {
        xx2 <- topN(xx, groupBy = fData(xx)$ProteinAccession, n = .n)
        tmp <- table(fData(xx)$ProteinAccession)
        tmp[tmp > .n] <- .n
        expect_equal(sum(tmp), nrow(xx2))
        expect_equal(ncol(xx), ncol(xx2))
    }
    xx2 <- topN(xx, groupBy = fData(xx)$ProteinAccession, n = 10)
    xo <- featureNames(xx)
    expect_equal(exprs(xx)[xo, ], exprs(xx2)[xo, ])
    expect_equal(fData(xx)[xo, ], fData(xx2)[xo, ])
    expect_equal(pData(xx)[xo, ], pData(xx2)[xo, ])
})
