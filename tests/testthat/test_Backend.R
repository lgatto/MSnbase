context("Backend class")

test_that("validity", {
    b <- BackendMemory()
    b@files <- c(F1="foo", F2="bar")
    expect_true(validObject(b))

    b@files[2] <- "foo"
    expect_error(validObject(b), "Duplicated")
})

test_that(".valid.Backend.files", {
    expect_null(.valid.Backend.files(character()))
    expect_null(.valid.Backend.files(c(F1="foo", F2="bar")))
    expect_match(.valid.Backend.files(c(F1="foo", F2="")), " missing")
    expect_match(.valid.Backend.files(c(F1="foo", F2=NA, F3="bar")), " NA")
    expect_match(.valid.Backend.files(c(F1="foo", F2="foo", F3="bar")),
                 "Duplicated")
    expect_match(.valid.Backend.files(c("foo", "bar")), " missing")
    expect_match(.valid.Backend.files(c(F1="foo", F1="bar")),
                 "Duplicated")
    expect_match(.valid.Backend.files(c(S1="foo", S2="bar")),
                 "don't start")
})

test_that("fileNames", {
    b <- BackendMemory()
    b@files <- c(F1="foo", F2="bar")
    expect_identical(fileNames(b), c(F1="foo", F2="bar"))
})

test_that("show", {
    b <- BackendMemory()
    b@files <- c(F1="foo", F2="bar")
    expect_output(show(b), paste("Backend: BackendMemory ", "Source files:",
                                 "  foo", "  bar", sep="\\n"))
})

test_that("backendSubset,Backend works", {
    be <- BackendMzR()
    be@files <- c("a", "b", "c", "d")
    names(be@files) <- paste0("F", 1:4)
    res <- backendSubset(be, i = 3, file = c(3, 1))
    expect_equal(unname(res@files), c("c", "a"))
    expect_equal(be, backendSubset(be, i = 3, file = 1:4))
})
