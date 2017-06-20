test_that("[,Chromatograms works", {

    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    colnames(chrs) <- c("a", "b")
    
    ## Subset using indices
    ## extract single element
    expect_true(is(chrs[1, 1], "Chromatogram"))
    expect_true(is(chrs[1, 1, drop = FALSE], "Chromatograms"))
    expect_equal(chrs[1, 2], ch2)
    ## extract a row
    expect_equal(chrs[1, ], list(a = ch, b = ch2))
    expect_equal(unname(chrs[1, , drop = FALSE]),
                 Chromatograms(list(ch, ch2), nrow = 1))    
    ## extract a column
    expect_equal(chrs[, 2], list(ch2, ch3))
    expect_equal(unname(chrs[, 2, drop = FALSE]),
                 Chromatograms(list(ch2, ch3), ncol = 1))

    ## Subset using logical
    expect_true(is(chrs[c(TRUE, FALSE), c(TRUE, FALSE)], "Chromatogram"))
    expect_true(is(chrs[c(TRUE, FALSE), c(TRUE, FALSE), drop = FALSE],
                   "Chromatograms"))
    expect_equal(chrs[c(TRUE, FALSE), c(FALSE, TRUE)], ch2)
    ## extract a row
    expect_equal(chrs[c(TRUE, FALSE), ], list(a = ch, b = ch2))
    expect_equal(unname(chrs[c(TRUE, FALSE), , drop = FALSE]),
                 Chromatograms(list(ch, ch2), nrow = 1))
    ## extract a column
    expect_equal(chrs[, c(FALSE, TRUE)], list(ch2, ch3))
    expect_equal(unname(chrs[, c(FALSE, TRUE), drop = FALSE]),
                 Chromatograms(list(ch2, ch3), ncol = 1))
        
    ## Subset using names
    expect_equal(chrs[, "a"], list(ch, ch1))
    expect_equal(unname(chrs[, "a", drop = FALSE]),
                 Chromatograms(list(ch, ch1), ncol = 1))
})

test_that("[<-,Chromatograms works", {

    ints <- abs(rnorm(12, sd = 20))
    ch <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(20, sd = 14))
    ch1 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(14, sd = 24))
    ch2 <- Chromatogram(rtime = 1:length(ints), ints)
    ints <- abs(rnorm(40, sd = 34))
    ch3 <- Chromatogram(rtime = 1:length(ints), ints)
    chrs <- Chromatograms(list(ch, ch1, ch2, ch3), nrow = 2)
    colnames(chrs) <- c("a", "b")

    ints <- abs(rnorm(94, sd = 200))
    ch4 <- Chromatogram(rtime = 1:length(ints), ints)

    ## errors
    expect_error(chrs[1:2, 1:2] <- list(ch4, ch4, ch4, ch4))
    expect_error(chrs["z", ] <- list(ch4, ch4))

    ## Single element.
    chrs[1, 2] <- ch4
    expect_equal(chrs[1, 2], ch4)
    chrs[, 2] <- list(ch2, ch3)
    expect_equal(chrs[, 2], list(ch2, ch3))
    chrs[2, 1] <- list(ch4)
    expect_equal(chrs[2, 1], ch4)
    
    chrs[, "a"] <- list(ch2, ch1)
    expect_equal(chrs[, 1], list(ch2, ch1))
    expect_error(chrs[, 1] <- list(ch, ch2, ch3))
    
    chrs[, c(TRUE, FALSE)] <- list(ch4, ch4)
    expect_equal(chrs[, 1], list(ch4, ch4))
    
})

