context("fragments")

test_that("calculateFragments", {
  pqr <- data.frame(
    mz = c(70.065, 198.124, 354.225,   # a
           98.060, 226.119, 382.220,   # b
           115.087, 243.145, 399.246,  # c
           201.098, 329.157, 426.210,  # x
           175.119, 303.178, 400.230,  # y
           158.092, 286.151, 383.204,  # z
           142.121, 270.180, 367.233,  # y_
           365.193,                    # b*
           286.151, 383.204),          # y*
    ion = c(paste0(rep(c("a", "b", "c", "x", "y", "z"), each=3),
                   rep(1:3, times=6)),
            paste0("y", 1:3, "_"), "b3*", "y2*", "y3*"),
    type = c(rep(c("a", "b", "c", "x", "y", "z", "y_"), each=3),
             "b*", "y*", "y*"),
    pos = c(rep(1:3, 7), 3, 2, 3),
    z = 1,
    seq = c(rep(c("P", "PQ", "PQR"), 3),
            rep(c("R", "QR", "PQR"), 4),
            "PQR", "QR", "PQR"),
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

  expect_message(calculateFragments("PQR", verbose = TRUE),
                 "Modifications used: C=57.02146")
  expect_message(calculateFragments("PQR", modifications = NULL,
                                    verbose = TRUE),
                 "Modifications used: None")

  expect_equal(pqr[1:18,],
               calculateFragments("PQR",
                                  type = c("a", "b", "c", "x", "y", "z"),
                                  neutralLoss = NULL, verbose = FALSE),
               tolerance=1e-5)
  expect_equal(pqr[1:6,],
               calculateFragments("PQR", type = c("a", "b"),
                                  neutralLoss = NULL, verbose = FALSE),
               tolerance=1e-5)
  ## rownames always differ
  expect_equal(pqr[c(10:12, 16:18),],
               calculateFragments("PQR", type = c("x", "z"),
                                  neutralLoss = NULL, verbose = FALSE),
               check.attributes = FALSE, tolerance = 1e-5)

  ## neutral loss
  ## rownames always differ
  expect_equal(pqr[c(4:6, 13:15, 19:24),],
               calculateFragments("PQR", verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss (water=cterm disabled),
  ## rownames always differ
  expect_equal(pqr[c(4:6, 13:15, 22:24),],
               calculateFragments("PQR",
                 neutralLoss=defaultNeutralLoss(disableWaterLoss="Cterm"),
                 verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss (ammonia=Q disabled),
  ## rownames always differ
  expect_equal(pqr[c(4:6, 13:15, 19:21),],
               calculateFragments("PQR",
                 neutralLoss=defaultNeutralLoss(disableAmmoniaLoss="Q"),
                 verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss + nterm mod, rownames always differ
  tpqr <- pqr[c(4:6, 13:15, 19:24),]
  tpqr$mz[c(1:3, 10)] <- tpqr$mz[c(1:3, 10)]+229
  expect_equal(tpqr,
               calculateFragments("PQR", modifications=c(C=57.02146, Nterm=229),
                                                         verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss + nterm + cterm mod, rownames always differ
  tpqr$mz[c(4:9, 11:12)] <- tpqr$mz[c(4:9, 11:12)]-100
  expect_equal(tpqr,
               calculateFragments("PQR", modifications=c(C=57.02146,
                                                         Nterm=229,
                                                         Cterm=-100),
                                                         verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  expect_equal(ace,
               calculateFragments("ACE", type=c("a", "b", "c", "x", "y", "z"),
                                  z=2, neutralLoss=NULL, verbose=FALSE),
               tolerance=1e-5)
  expect_equal(ace[1:9,],
               calculateFragments("ACE", type=letters[1:3], z=2, verbose=FALSE),
               tolerance=1e-5)

})

test_that("defaultNeutralLoss", {
  expect_equal(defaultNeutralLoss(),
               list(water=c("Cterm", "D", "E", "S", "T"),
                    ammonia=c("K", "N", "Q", "R")))
  expect_equal(defaultNeutralLoss(disableWaterLoss=c("T", "E", "S", "D")),
               list(water=c("Cterm"), ammonia=c("K", "N", "Q", "R")))
  expect_equal(defaultNeutralLoss(disableWaterLoss=c("T", "E", "S", "D"),
                                  disableAmmoniaLoss=c("K", "Q")),
               list(water=c("Cterm"), ammonia=c("N", "R")))
  expect_equal(defaultNeutralLoss(disableWaterLoss=c("Cterm",
                                                     "T", "E", "S", "D"),
                                  disableAmmoniaLoss=c("K", "N", "Q", "R")),
               list(water=character(), ammonia=character()))
})
