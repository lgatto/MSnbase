context("MSnSet class")

test_that("MSnSet validity", {
  expect_true(validObject(new("MSnSet")))
})


test_that("Combine MSnSet features", {
  aa <- new("MSnSet",
            exprs=matrix(c(rnorm(10,4,0.0001),rnorm(10,10,0.0001)),
              nrow=10,byrow=TRUE),
            featureData=new("AnnotatedDataFrame",
              data=data.frame(
                A=rep(c("A","B"),each=5),
                B=paste(rep(c("A","B"),each=5),1:10,sep="."))))
  bb <- combineFeatures(aa,factor(rep(letters[1:2],each=5)),"mean")
  cc <- combineFeatures(aa,factor(rep(letters[1:2],each=5)),"sum")
  dd <- combineFeatures(aa,factor(rep(letters[1:2],each=5)),"median")
  ee <- combineFeatures(aa,factor(rep(letters[1:2],each=5)), "weighted.mean", w=rep(1,10))
  ff <- combineFeatures(aa,factor(rep(letters[1:2], each = 5)), "medpolish", verbose=FALSE)
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
})

test_that("Purity correction", {
  file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
  aa <- readMSData(file,verbose=FALSE)
  msnset <- quantify(aa, method="trap", reporters = iTRAQ4, verbose = FALSE)
  impurity0 <- diag(4)
  pc <- purityCorrect(msnset, impurity0)
  expect_true(all(exprs(pc) == exprs(msnset)))
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
