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

