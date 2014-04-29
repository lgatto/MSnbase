context("Features of interest")

test_that("FeaturesOfInterest", {
    library("pRolocdata")
    data(tan2009r1)
    fn <- featureNames(tan2009r1)[1:10]
    desc <- "my description"
    foi1 <- FeaturesOfInterest(fn, desc, tan2009r1)
    foi2 <- FeaturesOfInterest(fn, desc)
    
    expect_equal(length(foi1), length(fn))
    expect_equal(length(foi2), length(fn))

    col <- FoICollection()
    expect_equal(length(col), 0)
    col <- addFeaturesOfInterest(foi1, col)
    expect_equal(length(col), 1)
    
    expect_message(col0 <- addFeaturesOfInterest(foi1, col),
                   "The features of interest are already present.")
    expect_true(all.equal(col, col0))

    expect_message(col0 <- addFeaturesOfInterest(foi2, col),
                   "The features of interest are already present.")
    expect_true(all.equal(col, col0))
})


test_that("FeaturesOfInterest traceable and fnamesIn", {
    library("pRolocdata")
    data(tan2009r1)
    fn <- featureNames(tan2009r1)[1:10]
    desc <- "my description"
    expect_error(FeaturesOfInterest(c("AA", fn), desc, tan2009r1))
    f <- FeaturesOfInterest(c("AA", fn), desc)
    expect_equal(length(f), 11)
    expect_true(fnamesIn(f, tan2009r1))
    expect_equal(fnamesIn(f, tan2009r1, count = TRUE), 10)
    expect_error(FeaturesOfInterest(letters[1:5], desc, tan2009r1))    
})
