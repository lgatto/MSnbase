context("Backend class")

test_that("validity", {
    b <- BackendMemory()
    b@files <- c("foo", "bar")
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
