context("fragments")

test_that("calculateFragments", {
  pqr <- data.frame(
    mz = c(70.065, 198.124,   # a
           98.060, 226.119,   # b
           115.087, 243.145,  # c
           201.098, 329.157,  # x
           175.119, 303.178,  # y
           158.092, 286.151,  # z
           157.108, 285.167,  # y_
           286.151),          # y*
    ion = c(paste0(rep(c("a", "b", "c", "x", "y", "z"), each=2),
                   rep(1:2, times=6)),
            paste0("y", 1:2, "_"), "y2*"),
    type = c(rep(c("a", "b", "c", "x", "y", "z", "y_"), each=2),
             "y*"),
    pos = c(rep(1:2, 7), 2),
    z = 1,
    seq = c(rep(c("P", "PQ"), 3),
            rep(c("R", "QR"), 4),
            "QR"),
    stringsAsFactors=FALSE)

  ace <- data.frame(
    mz = c(22.528, 102.544,  # a
           36.526, 116.541,  # b
           45.039, 125.054,  # c
           87.523, 167.539,  # x
           74.534, 154.549,  # y
           66.021, 146.036), # z
    ion = paste0(rep(c("a", "b", "c", "x", "y", "z"), each=2),
                 rep(1:2, times=6)),
    type = rep(c("a", "b", "c", "x", "y", "z"), each=2),
    pos = rep(1:2, 6),
    z = 2,
    seq = c(rep(c("A", "AC"), 3),
            rep(c("E", "CE"), 3)),
    stringsAsFactors=FALSE)

  expect_message(calculateFragments("PQR", verbose = TRUE),
                 "Modifications used: C=57.02146")
  expect_message(calculateFragments("PQR", modifications = NULL,
                                    verbose = TRUE),
                 "Modifications used: None")

  expect_equal(pqr[1:12,],
               calculateFragments("PQR",
                                  type = c("a", "b", "c", "x", "y", "z"),
                                  neutralLoss = NULL, verbose = FALSE),
               tolerance=1e-5)
  expect_equal(pqr[1:4,],
               calculateFragments("PQR", type = c("a", "b"),
                                  neutralLoss = NULL, verbose = FALSE),
               tolerance=1e-5)
  ## rownames always differ
  expect_equal(pqr[c(7:8, 11:12),],
               calculateFragments("PQR", type = c("x", "z"),
                                  neutralLoss = NULL, verbose = FALSE),
               check.attributes = FALSE, tolerance = 1e-5)

  ## neutral loss
  ## rownames always differ
  expect_equal(pqr[c(3:4, 9:10, 13:15),],
               calculateFragments("PQR", verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss (water=cterm disabled),
  ## rownames always differ
  expect_equal(pqr[c(3:4, 9:10, 15),],
               calculateFragments("PQR",
                 neutralLoss=defaultNeutralLoss(disableWaterLoss="Cterm"),
                 verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss (ammonia=Q disabled),
  ## rownames always differ
  expect_equal(pqr[c(3:4, 9:10, 13:14),],
               calculateFragments("PQR",
                 neutralLoss=defaultNeutralLoss(disableAmmoniaLoss="Q"),
                 verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss + nterm mod, rownames always differ
  tpqr <- pqr[c(3:4, 9:10, 13:15),]
  tpqr$mz[1:2] <- tpqr$mz[1:2] + 229
  expect_equal(tpqr,
               calculateFragments("PQR", modifications=c(C=57.02146, Nterm=229),
                                                         verbose=FALSE),
               check.attributes=FALSE, tolerance=1e-5)

  ## neutral loss + nterm + cterm mod, rownames always differ
  tpqr$mz[3:7] <- tpqr$mz[3:7] - 100
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
  expect_equal(ace[1:6,],
               calculateFragments("ACE", type=letters[1:3], z=2, verbose=FALSE),
               tolerance=1e-5)

  expect_error(calculateFragments("A"), "two or more residues")

  ## issue #200 (mz are not calculated correctly for terminal modifications
  ## and z > 1)
  p <- MSnbase:::get.atomic.mass()["p"]
  expect_equal(calculateFragments("AA", z=2,
                                  modifications=c(Nterm=10),
                                  type="b")$mz - p,
               (calculateFragments("AA", z=1,
                                   modifications=c(Nterm=10),
                                   type="b")$mz - p )/ 2)
  expect_equal(calculateFragments("AA", z=2, neutralLoss=NULL,
                                  modifications=c(Cterm=10),
                                  type="y")$mz - p,
               (calculateFragments("AA", z=1, neutralLoss=NULL,
                                   modifications=c(Cterm=10),
                                   type="y")$mz - p) / 2)

    ## issue #573 (charge is ignored in neutral loss calculation)
    expect_equal(
        subset(calculateFragments("PEPTIDEE", z = 3, type = "b"), pos == 7L),
        data.frame(
            mz = c(261.4570693, 255.4535476),
            ion = c("b7", "b7_"),
            type = c("b", "b_"),
            pos = 7L,
            z = 3,
            seq = "PEPTIDE"
        ),
        check.attributes = FALSE # row.names differ
    )
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
