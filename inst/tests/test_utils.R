context("utils")

test_that("vec2ssv & ssv2vec", {
  numbers <- 1:3
  string <- "1;2;3"

  expect_equal(MSnbase:::utils.vec2ssv(numbers), string)

  expect_equal(as.numeric(MSnbase:::utils.ssv2vec(string)), numbers)
})

test_that("list2ssv & ssv2list", {
  l <- list(a=1:3, b=4:6)
  string <- c("1;2;3;4;5;6")
  strings <- c(a="1;2;3", b="4;5;6")

  expect_equal(MSnbase:::utils.list2ssv(l), strings)

  expect_equal(lapply(MSnbase:::utils.ssv2list(strings), as.numeric), l)

  expect_equal(MSnbase:::utils.vec2ssv(unlist(l)), string)
})
