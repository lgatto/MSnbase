context("MSnSet class")

test_that("MSnSet validity", {
  expect_true(validObject(new("MSnSet")))
})

test_that("MSnSet coersion", {
    library("pRolocdata")
    data(dunkley2006)
    expect_true(validObject(dunkley2006))
    es <- as(dunkley2006, "ExpressionSet")
    expect_true(validObject(es))
})


test_that("Combine MSnSet features", {
  aa <- new("MSnSet",
            exprs=matrix(
              c(rnorm(10, 4, 0.0001),
                rnorm(10, 10, 0.0001)),
              nrow = 10,byrow = TRUE),
            featureData = new("AnnotatedDataFrame",
              data = data.frame(
                A = rep(c("A","B"), each = 5),
                B = paste(
                  rep(c("A","B"), each = 5),
                  1:10,
                  sep = "."))))
  gb <- factor(rep(letters[1:2], each = 5))
  bb <- combineFeatures(aa, gb, "mean")
  cc <- combineFeatures(aa, gb, "sum")
  dd <- combineFeatures(aa, gb, "median")
  ee <- combineFeatures(aa, gb, "weighted.mean", w=rep(1,10))
  ff <- combineFeatures(aa, gb, "medpolish", verbose=FALSE)
  expect_equal(exprs(bb),
               matrix(c(4,10,4,10), ncol = 2),
               tolerance = .001,
               check.attributes = FALSE)
  expect_equal(exprs(dd),
               matrix(c(4,10,4,10), ncol = 2),
               tolerance = .001,
               check.attributes = FALSE)
  expect_equal(exprs(ee),
               matrix(c(4,10,4,10), ncol = 2),
               tolerance = .001,
               check.attributes = FALSE)
  expect_equal(exprs(ff),
               matrix(c(4,10,4,10), ncol = 2),
               tolerance = .001,
               check.attributes = FALSE)
  expect_equal(exprs(cc),
               matrix(c(5*4,5*10,5*4,5*10), ncol = 2),
               tolerance = .001,
               check.attributes = FALSE)
  expect_true(all(fData(bb)[,1] == c("A", "B")))
  expect_true(all(fData(bb)[,2] == c("A.1", "B.6")))
  gb2 <- factor(c("a", "c", "z", "a", "z", "b", "b", "a", "c", "a"))
  gb3 <- factor(rev(c("a", "c", "z", "a", "z", "b", "b", "a", "c", "a")))
  fData(aa)$Z <- gb2  
  zz <- combineFeatures(aa, gb2, fun = "sum")  
  zz3 <- combineFeatures(aa[10:1, ], gb3, fun = "sum")
  expect_true(all.equal(exprs(zz), exprs(zz3)))
  expect_true(all.equal(fData(zz)[, 3:5],
                        fData(zz3)[, 3:5]))
  expect_true(all.equal(as.numeric(exprs(zz["a", ])),
                        colSums(exprs(aa)[gb2 == "a", ]),
                        check.attributes = FALSE))
  expect_true(all.equal(as.numeric(exprs(zz["b", ])),
                        colSums(exprs(aa)[gb2 == "b", ]),
                        check.attributes = FALSE))
  expect_true(all.equal(as.numeric(exprs(zz["c", ])),
                        colSums(exprs(aa)[gb2 == "c", ]),
                        check.attributes = FALSE))
  expect_true(all.equal(as.numeric(exprs(zz["z", ])),
                        colSums(exprs(aa)[gb2 == "z", ]),
                        check.attributes = FALSE))
  expect_true(all.equal(featureNames(zz),
                        as.character(fData(zz)$Z),
                        check.attributes = FALSE))  
})

test_that("Purity correction", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMSData(file,verbose=FALSE)
  msnset <- quantify(aa, method="trap", reporters = iTRAQ4, verbose = FALSE)
  impurity0 <- diag(4)
  pc <- purityCorrect(msnset, impurity0)
  expect_true(all(exprs(pc) == exprs(msnset)))  
})

test_that("makeImpuritiesMatrix", {
  i4 <- dir(system.file("extdata", package = "MSnbase"),
            pattern = "iTRAQ4plexPurityCorrection",
            full.names = TRUE)
  m4 <- makeImpuritiesMatrix(filename = i4, edit = FALSE)
  a4 <- matrix(c(0.929,0.059,0.002,0.000,
                 0.020,0.923,0.056,0.001,
                 0.000,0.030,0.924,0.045,
                0.000,0.001,0.040,0.923),
               nrow=4, byrow = TRUE)
  expect_equal(a4, m4, check.attributes = FALSE)
  t6 <- dir(system.file("extdata", package = "MSnbase"),
            pattern = "TMT6plexPurityCorrection",
            full.names = TRUE)
  m6 <- makeImpuritiesMatrix(filename = t6, edit = FALSE)
  a6 <- matrix(c(0.939, 0.061, 0.000, 0.000, 0.000, 0.000,
                 0.005, 0.928, 0.067, 0.000, 0.000, 0.000,
                 0.000, 0.011, 0.947, 0.042, 0.000, 0.000,
                 0.000, 0.000, 0.017, 0.942, 0.041, 0.000,
                 0.000, 0.000, 0.000, 0.016, 0.963, 0.021,
                 0.000, 0.000, 0.000, 0.002, 0.032, 0.938),
               nrow = 6, byrow = TRUE)
  expect_equal(a6, m6, check.attributes = FALSE)
})

test_that("Normalisation and transpose", {
  bb <- quantify(itraqdata, method="trap", reporters=iTRAQ4, verbose=FALSE)
  bb1 <- normalise(bb, "sum")
  expect_true(all.equal(rowSums(exprs(bb1), na.rm=TRUE),
                        rep(1,nrow(bb1)), check.attributes=FALSE))
  bb2 <- normalise(bb,"max")
  expect_true(all(apply(exprs(bb2), 1, max, na.rm=TRUE) == 1))
  bb3 <- normalise(bb, "quantiles")
  bb4 <- normalise(bb, "quantiles.robust")
  bb5 <- normalise(bb, "vsn")
})


test_that("Transpose and subset", {
  aa <- quantify(itraqdata, method="trap", reporters=iTRAQ4, verbose=FALSE)
  ## transpose
  ##expect_warning(taa <- t(aa),"Dropping protocolData.") ## replaced by message()
  taa <- t(aa)
  expect_true(nrow(aa) == ncol(taa)) 
  expect_true(ncol(aa) == nrow(taa))
  expect_true(all.equal(pData(aa), fData(taa)))
  expect_true(all.equal(pData(taa), fData(aa)))
  ## subset
  bb <- aa[1:2,c(2,4)]
  ## expect_true(all(dim(qual(bb)) == c(4,7)))
  ## expect_true(all(qual(bb)$reporter %in% bb$mz))
})
