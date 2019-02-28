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

test_that("backendSplitByFile,Backend works", {
    b <- BackendMzR()
    b@files <- c("a", "b", "c")
    b@modCount <- rep(0L, 3L)
    spd <- DataFrame(fileIdx = c(3, 3, 1, 1, 2, 1))
    res <- backendSplitByFile(b, spd)
    bl <- BackendMzR()
    bl@files <- "a"
    bl@modCount <- 0L
    l <- list("1"=bl, "2"=bl, "3"=bl)
    l[[2]]@files <- "b"
    l[[3]]@files <- "c"
    expect_equal(backendSplitByFile(b, spd), l)
    r <- b
    r@files[1] <- "d"
    r@modCount[1L] <- 1L
    l[[1]]@files <- "d"
    l[[1]]@modCount <- 1L
    backendSplitByFile(b, spd) <- l
    expect_equal(b, r)
})
