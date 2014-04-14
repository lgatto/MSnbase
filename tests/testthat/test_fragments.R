context("fragments")

test_that("calculateFragments", {
  ymass <- c(175.119, 303.1775, 400.2303)
  bmass <- c(98.06, 226.1186, 382.2197)
  mass <- c(bmass, ymass)
  o <- order(mass)
  mass <- mass[o]
  nm <- paste0(rep(c("b", "y"), each=3), rep(1:3, times=2))[o]
  res <- MSnbase:::calculateFragments("PQR")

  expect_equal(nm, res$fragment.str)
  expect_true(all(abs(res$mass - mass) < 1e-4))
  expect_true(all(res$fragment.seq == c("P", "R", "PQ", "QR", "PQR", "PQR")))
})

