library("pRolocdata")

context("Features of interest")

test_that("FeaturesOfInterest", {
    data(tan2009r1)

    fn <- featureNames(tan2009r1)[1:3]
    foi0 <- FeaturesOfInterest(fn, "small foi")
    expect_null(show(foi0))
    expect_equal(length(foi0@objpar), 0L) ## not traceable
    
    fn <- featureNames(tan2009r1)[1:10]
    desc <- "my description"
    foi1 <- FeaturesOfInterest(fn, desc, tan2009r1)
    expect_null(show(foi1))
    expect_equal(length(foi1@objpar), 4L) ## traceable
    expect_identical(description(foi1), desc)
    foi2 <- FeaturesOfInterest(fn, desc)
    
    expect_equal(length(foi1), length(fn))
    expect_equal(length(foi2), length(fn))

    col <- FoICollection()    
    expect_equal(length(col), 0)
    col <- addFeaturesOfInterest(foi1, col)
    expect_equal(length(col), 1)
    expect_null(show(col))

    expect_message(col1 <- addFeaturesOfInterest(foi2, col),
                   "The features of interest are already present.")
    expect_identical(col, col1)
    
    col2 <- FoICollection(list(foi1))
    expect_identical(col, col2)
    expect_identical(lengths(col), 10L)

    col3 <- FoICollection(list(foi1, foi2))
    expect_identical(foi(col3, 1), list(foi1))
    expect_error(foi(col3, 3),
                 "There are only 2 available FeatureOfInterest instances.")
    
    expect_identical(description(col3), rep(desc, 2))
    
    expect_message(col0 <- addFeaturesOfInterest(foi1, col),
                   "The features of interest are already present.")
    expect_true(all.equal(col, col0))

    expect_message(col0 <- addFeaturesOfInterest(foi2, col),
                   "The features of interest are already present.")
    expect_true(all.equal(col, col0))

    expect_identical(rmFeaturesOfInterest(col3, 2), col2)


    foil <- lapply(replicate(10, sample(nrow(tan2009r1), 15), simplify=FALSE),
                   function(i) FeaturesOfInterest(featureNames(tan2009r1)[i], "my fois"))
    col4 <- FoICollection(foil)
    expect_identical(lengths(col4), rep(15L, 10))
    expect_null(show(col4))

    foim <- as(col4, "matrix")
    expect_identical(foim, as.matrix.FoICollection(col4))
    expect_identical(length(col4), ncol(foim))

    ufns <- unique(unlist(lapply(foi(col4), foi)))
    expect_identical(ufns, rownames(foim))                     
})



test_that("FeaturesOfInterest traceable and fnamesIn", {
    data(tan2009r1)
    fn <- featureNames(tan2009r1)[1:10]
    desc <- "my description"
    expect_error(FeaturesOfInterest(c("AA", fn), desc, tan2009r1))
    f <- FeaturesOfInterest(c("AA", fn), desc)
    expect_equal(length(f), 11)
    expect_true(fnamesIn(f, tan2009r1))
    expect_true(fnamesIn(f, exprs(tan2009r1)))
    tandfr <- as(t(tan2009r1), "data.frame")
    expect_true(fnamesIn(f, tandfr))
    
    expect_equal(fnamesIn(f, tan2009r1, count = TRUE), 10)
    expect_error(FeaturesOfInterest(letters[1:5], desc, tan2009r1))    
})

