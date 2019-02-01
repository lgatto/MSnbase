test_that("setMSnbaseFastLoad works", {
    orig_value <- MSnbaseOptions()$fastLoad
    setMSnbaseFastLoad(!orig_value)
    expect_equal(MSnbaseOptions()$fastLoad, !orig_value)
    setMSnbaseFastLoad(orig_value)
    expect_equal(MSnbaseOptions()$fastLoad, orig_value)
})

test_that("isMSnbaseFastLoad works", {
    orig_value <- MSnbaseOptions()$fastLoad
    expect_equal(isMSnbaseFastLoad(), orig_value)
    setMSnbaseFastLoad(!orig_value)
    expect_equal(isMSnbaseFastLoad(), !orig_value)
    setMSnbaseFastLoad(orig_value)
})

test_that("Hdf5CompressionLevel works", {
    orig_value <- MSnbaseOptions()$HDF5_COMP_LEVEL
    expect_equal(.hdf5_compression_level(), orig_value)
    setHdf5CompressionLevel(9)
    expect_equal(.hdf5_compression_level(), 9L)
    expect_error(setHdf5CompressionLevel("k"))
    expect_error(setHdf5CompressionLevel(12))
    setHdf5CompressionLevel(orig_value)
})
