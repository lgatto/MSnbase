test_that(".hdf5_group_name works", {
    res <- .hdf5_group_name(c("a/bb b/ccc", "a/bbb/ccc"))
    expect_equal(length(res), 2)
    expect_true(length(unique(res)) == 2)
    res <- .hdf5_group_name(
        c("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa/bb b/ccc",
          "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa/bbb/ccc"))
    expect_equal(length(res), 2)
    expect_true(length(unique(res)) == 2)
})
