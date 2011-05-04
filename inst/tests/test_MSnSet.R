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
  bb <- MSnbase:::combineFeatures(aa,factor(rep(letters[1:2],each=5)),"mean")
  cc <- MSnbase:::combineFeatures(aa,factor(rep(letters[1:2],each=5)),"sum")
  dd <- MSnbase:::combineFeatures(aa,factor(rep(letters[1:2],each=5)),"median")
  ee <- MSnbase:::combineFeatures(aa,factor(rep(letters[1:2],each=5)),
                                  "weighted.mean",w=rep(1,10))
  ff <- MSnbase:::combineFeatures(aa,factor(rep(letters[1:2],each=5)),
                                  "medpolish",verbose=FALSE)
  expect_equal(exprs(bb),matrix(c(4,10,4,10),ncol=2),
               tolerance=.001,check.attributes=FALSE)
  expect_equal(exprs(dd),matrix(c(4,10,4,10),ncol=2),
               tolerance=.001,check.attributes=FALSE)
  expect_equal(exprs(ee),matrix(c(4,10,4,10),ncol=2),
               tolerance=.001,check.attributes=FALSE)
  expect_equal(exprs(ff),matrix(c(4,10,4,10),ncol=2),
               tolerance=.001,check.attributes=FALSE)
  expect_equal(exprs(cc),matrix(c(5*4,5*10,5*4,5*10),ncol=2),
               tolerance=.001,check.attributes=FALSE)
  expect_true(all(fData(bb)[,1]==c("A","B")))
  expect_true(all(fData(bb)[,2]==c("A.1","B.6")))
})
