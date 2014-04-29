context("fragments")

test_that("calculateFragments", {
  pqr <- data.frame(
    mz = c(70.065, 198.124, 354.225,   # a
           98.060, 226.119, 382.220,   # b
           115.087, 243.145, 399.246,  # c
           201.098, 329.157, 426.210,  # x
           175.119, 303.178, 400.230,  # y
           158.092, 286.151, 383.204), # z
    ion = paste0(rep(c("a", "b", "c", "x", "y", "z"), each=3),
                 rep(1:3, times=6)),
    type = rep(c("a", "b", "c", "x", "y", "z"), each=3),
    pos = rep(1:3, 6),
    z = 1,
    seq = c(rep(c("P", "PQ", "PQR"), 3),
            rep(c("R", "QR", "PQR"), 3)),
    stringsAsFactors=FALSE)

  ace <- data.frame(
    mz = c(22.528, 102.544, 167.065,  # a
           36.526, 116.541, 181.062,  # b
           45.039, 125.054, 189.576,  # c
           87.523, 167.539, 203.057,  # x
           74.534, 154.549, 190.068,  # y
           66.021, 146.036, 181.554), # z
    ion = paste0(rep(c("a", "b", "c", "x", "y", "z"), each=3),
                 rep(1:3, times=6)),
    type = rep(c("a", "b", "c", "x", "y", "z"), each=3),
    pos = rep(1:3, 6),
    z = 2,
    seq = c(rep(c("A", "AC", "ACE"), 3),
            rep(c("E", "CE", "ACE"), 3)),
    stringsAsFactors=FALSE)

  pqr. <- MSnbase:::calculateFragments("PQR", type=c("a", "b", "c", "x", "y", "z"))
  pqr.$mz <- round(pqr.$mz, 3)

  ace. <- MSnbase:::calculateFragments("ACE", type=c("a", "b", "c", "x", "y", "z"), z=2)
  ace.$mz <- round(ace.$mz, 3)

  expect_equal(pqr, pqr.)
  expect_equal(ace, ace.)
})

