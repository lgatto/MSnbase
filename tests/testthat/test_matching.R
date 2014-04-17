context("matching")

test_that("relaxedMatch", {
  x <- c(3, 4.2, 5, 7, 2, 8, 9, 6, 10, 1)
  table <- c(1.08, 3.95, 5.11, 11.18, 15.88)

  expect_error(MSnbase:::relaxedMatch(table, x), "sorted non-decreasingly")
  expect_error(MSnbase:::relaxedMatch(x, table, tolerance=2), "must be smaller than 1")
  expect_warning(MSnbase:::relaxedMatch(x, table, tolerance=-2), ".*tolerance.* < 0 is meaningless")
  expect_equal(MSnbase:::relaxedMatch(x, table),
               as.integer(rep(NA, 10)))
  expect_equal(MSnbase:::relaxedMatch(x, table, tolerance=0.2),
               c(NA, 2, 3, NA, NA, NA, 4, 3, 4, 1))
  expect_equal(MSnbase:::relaxedMatch(x, table, tolerance=0.2, nomatch=-1),
               c(-1, 2, 3, -1, -1, -1, 4, 3, 4, 1))
  expect_equal(MSnbase:::relaxedMatch(x, table, tolerance=0.2, relative=FALSE),
               c(NA, NA, 3, NA, NA, NA, NA, NA, NA, 1))
  expect_equal(MSnbase:::relaxedMatch(x, table, tolerance=0.3, relative=FALSE),
               c(NA, 2, 3, NA, NA, NA, NA, NA, NA, 1))
  expect_equal(MSnbase:::relaxedMatch(c(4, 5), 4.8, tolerance=1.1,
                                      relative=FALSE), c(1, 1))
  expect_equal(MSnbase:::relaxedMatch(4.8, c(4, 5), tolerance=1.1,
                                      relative=FALSE), 2)
})

test_that("matchPeaks", {
  s1 <- new("Spectrum2",
            mz=c(1:3, 4.2, 5:10),
            intensity=c(1:4, 10, 15, 7:8, 15, 10))
  s2 <- new("Spectrum2",
            mz=c(1.08, 3.95, 5.11, 11.18, 15.88),
            intensity=c(1, 4, 5, 11, 16))
  expect_error(MSnbase:::matchPeaks(s1, s2, method="foobar"),
               "'arg' should be one")
  expect_equal(MSnbase:::matchPeaks(s1, s1), 1:10)
  expect_equal(MSnbase:::matchPeaks(s1, s2, tolerance=0.2),
               c(1, NA, NA, 2, NA, 3, NA, NA, 4, NA))
  expect_equal(MSnbase:::matchPeaks(s1, s2, tolerance=0.2, method="highest"),
               c(1, NA, NA, 2, NA, 3, NA, NA, 4, NA))
  expect_equal(MSnbase:::matchPeaks(s1, s2, tolerance=0.2, method="closest"),
               c(1, NA, NA, 2, 3, NA, NA, NA, NA, 4))
})

test_that("commonPeaks", {
  s1 <- new("Spectrum2",
            mz=c(1:3, 4.2, 5:10),
            intensity=c(1:4, 10, 15, 7:8, 15, 10))
  s2 <- new("Spectrum2",
            mz=c(1.08, 3.95, 5.11, 11.18, 15.88),
            intensity=c(1, 4, 5, 11, 16))
  expect_error(MSnbase:::commonPeaks(s1, s2, method="foobar"),
               "'arg' should be one")
  expect_equal(MSnbase:::commonPeaks(s1, s1), rep(TRUE, 10))
  expect_equal(MSnbase:::commonPeaks(s1, s2, tolerance=0.2),
               c(TRUE, FALSE, FALSE, TRUE, FALSE,
                 TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(MSnbase:::commonPeaks(s1, s2, tolerance=0.2, method="highest"),
               c(TRUE, FALSE, FALSE, TRUE, FALSE,
                 TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(MSnbase:::commonPeaks(s1, s2, tolerance=0.2, method="closest"),
               c(TRUE, FALSE, FALSE, TRUE, TRUE,
                 FALSE, FALSE, FALSE, FALSE, TRUE))
})

test_that("numberOfCommonPeaks", {
  s1 <- new("Spectrum2",
            mz=c(1:3, 4.2, 5:10),
            intensity=c(1:4, 10, 15, 7:8, 15, 10))
  s2 <- new("Spectrum2",
            mz=c(1.08, 3.95, 5.11, 11.18, 15.88),
            intensity=c(1, 4, 5, 11, 16))
  expect_true(MSnbase:::numberOfCommonPeaks(s1, s1) == 10)
  expect_true(MSnbase:::numberOfCommonPeaks(s1, s2, tolerance=0.2) == 4)
  expect_true(MSnbase:::numberOfCommonPeaks(s1, s2, tolerance=0.2) == 4)
})

