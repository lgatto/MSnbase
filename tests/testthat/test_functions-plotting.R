test_that(".plotXIC works", {
    x <- data.frame(rt = rep(abs(rnorm(30)), 100),
                    mz = rep(abs(rnorm(30)), 100),
                    i = abs(rnorm(300, mean = 200)))
    od1 <- filterFile(microtofq_on_disk, 1)

    x <- as(filterMz(filterRt(od1, c(270, 290)), c(610, 615)), "data.frame")
    .plotXIC(x)

    ## Test passing additional arguments
    .plotXIC(x, cex = 2)
    .plotXIC(x, pch = 1)
    .plotXIC(x, pch = 23)
})

test_that(".vertical_sub_layout works", {
    exp <- cbind(c(1, 2, 7, 8),
                 c(3, 4, 9, 10),
                 c(5, 6, 11, 12))
    expect_equal(.vertical_sub_layout(5, 2), exp)
    ## 8 plots (will be 3 columns), 3 plots each.
    exp <- cbind(c(1, 2, 3, 10, 11, 12, 19, 20, 21),
                 c(4, 5, 6, 13, 14, 15, 22, 23, 24),
                 c(7, 8, 9, 16, 17, 18, 25, 26, 27))
    expect_equal(.vertical_sub_layout(8, 3), exp)

    ## single plot, 5 each
    exp <- matrix(1:5, ncol = 1)
    expect_equal(.vertical_sub_layout(1, 5), exp)

    ## 19 plots 1 each
    exp <- matrix(1:20, ncol = 5, byrow = TRUE)
    expect_equal(.vertical_sub_layout(19, 1), exp)

    ## 19 plots 2 each
    exp <- cbind(c(1, 2, 11, 12, 21, 22, 31, 32),
                 c(3, 4, 13, 14, 23, 24, 33, 34),
                 c(5, 6, 15, 16, 25, 26, 35, 36),
                 c(7, 8, 17, 18, 27, 28, 37, 38),
                 c(9, 10, 19, 20, 29, 30, 39, 40))
    expect_equal(.vertical_sub_layout(19, 2), exp)
})
