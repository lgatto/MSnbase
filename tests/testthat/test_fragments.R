context("fragments")

test_that("calculateFragments", {
  pqr <- data.frame(
    mz = c(70.066, 198.124, 354.225,   # a
           98.061, 226.119, 382.220,   # b
           115.087, 243.146, 399.247,  # c
           201.099, 329.157, 426.210,  # x
           175.119, 303.178, 400.231,  # y
           158.093, 286.152, 383.204), # z
    ion = paste0(rep(c("a", "b", "c", "x", "y", "z"), each=3),
                 rep(1:3, times=6)),
    type = rep(c("a", "b", "c", "x", "y", "z"), each=3),
    pos = rep(1:3, 6),
    z = 1,
    seq = c(rep(c("P", "PQ", "PQR"), 3),
            rep(c("R", "QR", "PQR"), 3)),
    stringsAsFactors=FALSE)

  ace <- data.frame(
    mz = c(22.529, 102.544, 167.066,  # a
           36.526, 116.542, 181.063,  # b
           45.040, 125.055, 189.576,  # c
           87.524, 167.539, 203.058,  # x
           74.534, 154.550, 190.068,  # y
           66.021, 146.036, 181.555), # z
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

