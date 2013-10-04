context("Manual test: plotMzDelta")

test_that("plotMzDelta MSnExp and mzR", {
    f <- "~/Data2/Thermo_HELA_PRT/Thermo_Hela_PRTC_1.mzML"
    ## exp <- readMSData(f)
    ## save(exp, file = "Thermo_Hela_PRTC_1_MSnExp.rda")
    load("Thermo_Hela_PRTC_1_MSnExp.rda")
    
    ms <- openMSfile(f)
    hd <- header(ms)
    ms2 <- which(hd$msLevel == 2)
    hd2 <- hd[ms2, ]
    pl <- MSnbase:::peaksAsLists(ms, ms2)
    expect_equal(length(exp), length(pl))

    i <- which(hd2$acquisitionNum == acquisitionNum(exp[[1]]))
    print(system.time(d1 <- MSnbase:::utils.getMzDelta(exp[[1]], 0.1)))
    print(system.time(d2 <- MSnbase:::utils.getMzDelta_list(pl[[i]], 0.1)))
    expect_equal(d1, d2)

    d1 <- MSnbase:::utils.removePrecMz(exp[[1]])
    d2 <- MSnbase:::utils.removePrecMz_list(pl[[i]], hd2$precursorMZ[i])
    expect_equal(intensity(d1), d2$int)
    expect_equal(mz(d1), d2$mz)

    i <- match(acquisitionNum(exp[1:100]), hd2$acquisitionNum)
    print(system.time(p1 <- plotMzDelta(exp[1:100], plot = FALSE)))
    print(system.time(
        p2 <- MSnbase:::plotMzDelta_list(pl[i],
                                         precMz = hd2$precursorMZ[i],
                                         plot = FALSE)))
    expect_equal(p1, p2)
})
