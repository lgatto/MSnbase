context("Purity correction")

test_that("Purity correction identify", {
    file <- dir(system.file(package = "MSnbase", dir = "extdata"),
                full.name = TRUE, pattern = "mzXML$")
    aa <- readMSData(file, centroided = FALSE)
    bp <- SerialParam()
    msnset <- quantify(aa, method = "trap", reporters = iTRAQ4,
                       BPPARAM = bp)
    impurity0 <- diag(4)
    pc <- purityCorrect(msnset, impurity0)
    expect_true(all(exprs(pc) == exprs(msnset)))
})


test_that("Valid answer", {
              M <- matrix(c(1010, 1990), nrow = 1)
              colnames(M) <- LETTERS[1:2]
              rownames(M) <- "X"
              fd <- data.frame(XX = 1, row.names = "X")
              pd <- data.frame(YY = c(1:2), row.names = colnames(M))
              xx <- MSnSet(M, fd, pd)
              ## exprs(xx)
              ##      A    B
              ## X 1010 1990
              imp <- matrix(c(.97, 0.03,  0.02, 0.98),
                            ncol = 2, byrow = TRUE)
              ##      [,1] [,2]
              ## [1,] 0.97 0.03
              ## [2,] 0.02 0.98
              res <- exprs(purityCorrect(xx, imp))
              expect_equal(as.numeric(res), c(1000, 2000))
              ##      A    B
              ## X 1000 2000
          })
