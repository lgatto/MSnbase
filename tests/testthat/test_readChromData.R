test_that("readChromatograms works", {
    library(msdata)
    fl <- proteomics(full.names = TRUE, pattern = "MRM")
    files <- c(fl, fl, fl)
})

test_that(".combine_data.frame works", {
    A <- data.frame(a = c(1, 1, 2, 3, 4, 5, 6), b = c(2, 2, 3, 4, 5, 5, 6))
    B <- data.frame(a = c(1, 3, 3, 3, 3, 4, 5), b = c(2, 4, 4, 4, 5, 5, 5))
    C <- data.frame(a = c(3, 3, 4, 4), b = c(4, 5, 5, 5))
    
    .combine_data.frame(list(A, B, C))
    exp <- data.frame(a = c(1, 1, 2, 3, 4, 5, 6))
    
    expect_error(.combine_data.frame())
    expect_error(.combine_data.frame(list(A, B, C), cols = c("z")))

    ## When I do have factors?
    D <- data.frame(a = c("z", "b", "a", "g", "g"),
                    b = c(1, 2, 3, 4, 4), stringsAsFactors = FALSE)
    E <- data.frame(a = c("g", "a", "d", "g"),
                    b = c(4, 3, 1, 4), stringsAsFactors = FALSE)
    exp <- data.frame(a = c("a", "b", "d", "g", "g", "z"),
                      b = c(3, 2, 1, 4, 4, 1), stringsAsFactors = FALSE)
    expect_equal(.combine_data.frame(list(D, E)), exp)
    expect_equal(.combine_data.frame(list(E, D)), exp)
})
