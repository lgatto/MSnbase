context("Backend class")

test_that("validity", {
    b <- BackendMemory()
    b@files <- c("foo", "bar")
    expect_error(validObject(b), "counters")
    b@modCount <- 1L:2L
    expect_true(validObject(b))

    b@files[2] <- "foo"
    expect_error(validObject(b), "Duplicated")
})

test_that(".valid.Backend.files", {
    expect_null(.valid.Backend.files(character()))
    expect_null(.valid.Backend.files(c("foo", "bar")))
    expect_match(.valid.Backend.files(c("foo", "")), " missing")
    expect_match(.valid.Backend.files(c("foo", NA, "bar")), " NA")
    expect_match(.valid.Backend.files(c("foo", "foo", "bar")), "Duplicated")
    expect_null(.valid.Backend.files(c("foo", "bar")))
    expect_null(.valid.Backend.files(c("foo", "bar")))
    expect_null(.valid.Backend.files(c(S1="foo", S2="bar")))
})

test_that(".valid.Backend.files", {
    expect_null(.valid.Backend.modCount("foo", 5))
    expect_match(.valid.Backend.modCount("foo", 1:2), " counters")
})

test_that("fileNames", {
    b <- BackendMemory()
    b@files <- c("foo", "bar")
    expect_identical(fileNames(b), c("foo", "bar"))
})

test_that("show", {
    b <- BackendMemory()
    b@files <- c("foo", "bar")
    expect_output(show(b), paste("Backend: BackendMemory ", "Source files:",
                                 "  foo", "  bar", sep="\\n"))
})

test_that("backendSubset,Backend works", {
    be <- BackendMzR()
    be@files <- c("a", "b", "c", "d")
    be@modCount <- rep(0L, 4L)
    spd <- DataFrame(fileIdx = c(3, 3, 1, 3, 1, 1))
    res <- backendSubset(be, spd)
    expect_equal(unname(res@files), c("c", "a"))
    spd <- DataFrame(fileIdx = c(1, 2, 3, 4))
    expect_equal(be, backendSubset(be, spd))
})
