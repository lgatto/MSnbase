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

test_that("mergeSpectraAndIdentificationData", {
  ## pseudo fData(MSnSet) output
  fd <- data.frame(spectrum=1:4,
                   acquisitionNum=5:8, 
                   filename="foobar.mzML",
                   uselesscolumn=1,
                   stringsAsFactors=FALSE)
  ## pseudo mzID output
  id <- data.frame(acquisitionnum=c(5, 5, 5, 8), 
                   spectrumFile="foobar.mzML",
                   rank=c(2, 3, 1, 1),
                   accession=paste0("P", 1:4),
                   description=paste0("D", 1:4),
                   spectrumid=paste0("id", 1:4),
                   uselesscolumn=2,
                   stringsAsFactors=FALSE)
  rfd <- data.frame(filename="foobar.mzML",
                    acquisitionNum=5:8,
                    spectrum=1:4,
                    uselesscolumn.spectrum=1,
                    rank=c(1, NA, NA, 1),
                    accession=c("P3;P1;P2", NA, NA, "P4"),
                    description=c("D3;D1;D2", NA, NA, "D4"),
                    uselesscolumn.id=c(2, NA, NA, 2),
                    npsm=c(3, NA, NA, 1),
                    row.names=as.character(1:4),
                    stringsAsFactors=FALSE)
  expect_equal(MSnbase:::utils.mergeSpectraAndIdentificationData(fd, id), rfd)
})
