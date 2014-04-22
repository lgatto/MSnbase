context("utils")

test_that("vec2ssv & ssv2vec", {
  numbers <- 1:3
  string1 <- "1;2;3"
  string2 <- "1,2,3"

  expect_equal(MSnbase:::utils.vec2ssv(numbers), string1)
  expect_equal(MSnbase:::utils.vec2ssv(numbers, sep=","), string2)

  expect_equal(as.numeric(MSnbase:::utils.ssv2vec(string1)), numbers)
  expect_equal(as.numeric(MSnbase:::utils.ssv2vec(string2, sep=",")), numbers)
})

test_that("list2ssv & ssv2list", {
  l <- list(a=1:3, b=4:6)
  string1 <- c("1;2;3;4;5;6")
  string2 <- c("1,2,3,4,5,6")
  strings1 <- c(a="1;2;3", b="4;5;6")
  strings2 <- c(a="1,2,3", b="4,5,6")

  expect_equal(MSnbase:::utils.list2ssv(l), strings1)
  expect_equal(MSnbase:::utils.list2ssv(l, sep=","), strings2)

  expect_equal(lapply(MSnbase:::utils.ssv2list(strings1), as.numeric), l)
  expect_equal(lapply(MSnbase:::utils.ssv2list(strings2, sep=","), as.numeric), l)

  expect_equal(MSnbase:::utils.vec2ssv(unlist(l)), string1)
  expect_equal(MSnbase:::utils.vec2ssv(unlist(l), sep=","), string2)
})

test_that("leftJoin", {
  x <- data.frame(id=1:7,
                  fn=LETTERS[1:7],
                  add=1:7,
                  useless1=1,
                  stringsAsFactors=FALSE)
  rownames(x) <- paste0("R", seq(nrow(x)))
  y1 <- data.frame(id=4:2,
                   fn=LETTERS[4:2],
                   foobar=letters[4:2],
                   useless2=2,
                   stringsAsFactors=FALSE)
  y2 <- data.frame(id=6:7,
                   fn=LETTERS[6:7],
                   foobar=letters[6:7],
                   useless3=3,
                   stringsAsFactors=FALSE)
  z1 <- data.frame(id=1:7,
                   fn=LETTERS[1:7],
                   add=1:7,
                   foobar=c(NA, letters[2:4], rep(NA, 3)),
                   stringsAsFactors=FALSE)
  z2 <- data.frame(id=1:7,
                   fn=LETTERS[1:7],
                   add=1:7,
                   foobar=c(NA, letters[2:4], NA, letters[6:7]),
                   stringsAsFactors=FALSE)
  rownames(z1) <- rownames(z2) <- paste0("R", seq(nrow(x)))
  ## first run
  expect_equal(MSnbase:::utils.leftJoin(x, y1, by=c("id", "fn"),
                                        exclude=c("useless1", "useless2")), z1)
  ## second run (on the results of the first run)
  expect_equal(MSnbase:::utils.leftJoin(z1, y2, by=c("id", "fn"),
                                        exclude=c("useless3")), z2)
})

test_that("mergeSpectraAndIdentificationData", {
  ## pseudo fData(MSnSet) output
  fd <- data.frame(spectrum=1:4,
                   file=c(1, 2, 1, 1),
                   acquisition.number=5:8,
                   row.names=paste0("R", 1:4),
                   stringsAsFactors=FALSE)
  ## pseudo mzID output
  id1 <- data.frame(acquisitionnum=c(5, 5, 5, 8),
                    file=1,
                    spectrumFile="foobar1.mzML",
                    rank=c(2, 3, 1, 1),
                    accession=paste0("P", 1:4),
                    description=paste0("D", 1:4),
                    spectrumid=paste0("id", 1:4),
                    stringsAsFactors=FALSE)
  id2 <- data.frame(acquisitionnum=6,
                    file=2,
                    spectrumFile="foobar2.mzML",
                    rank=1,
                    accession="P9",
                    description="D9",
                    spectrumid="id9",
                    stringsAsFactors=FALSE)
  ## results
  rfd1 <- data.frame(spectrum=1:4,
                     file=c(1, 2, 1, 1),
                     acquisition.number=5:8,
                     rank=c(1, NA, NA, 1),
                     accession=c("P3;P1;P2", NA, NA, "P4"),
                     description=c("D3;D1;D2", NA, NA, "D4"),
                     row.names=paste0("R", 1:4),
                     stringsAsFactors=FALSE)
  rfd2 <- data.frame(spectrum=1:4,
                     file=c(1, 2, 1, 1),
                     acquisition.number=5:8,
                     rank=c(1, 1, NA, 1),
                     accession=c("P3;P1;P2", "P9", NA, "P4"),
                     description=c("D3;D1;D2", "D9", NA, "D4"),
                     row.names=paste0("R", 1:4),
                     stringsAsFactors=FALSE)
  ## first run
  expect_equal(MSnbase:::utils.mergeSpectraAndIdentificationData(fd, id1), rfd1)
  ## second run
  expect_equal(MSnbase:::utils.mergeSpectraAndIdentificationData(rfd1, id2),
               rfd2)
})

test_that("utils.idSummary", {
  ## pseudo fData(MSnSet) output
  fd <- data.frame(file=c(1, 2, 1, 1, 1, 3),
                   identFile=c(4, 5, 4, 4, NA, NA))
  ## results
  rdf <- data.frame(quantFile=c(1, 2, 3),
                    identFile=c(4, 5, NA),
                    coverage=c(0.75, 1, 0))
  expect_error(MSnbase:::utils.idSummary(data.frame(file=1:3, foobar=1:3)),
               "No quantification/identification data found")
  expect_equal(MSnbase:::utils.idSummary(fd), rdf)
})


test_that("formatRt", {
    tc <- c("1:1", "25:24")
    tn <- c(61, 25 * 60 + 24)
    expect_equal(tc, formatRt(tn)) 
    expect_equal(tn, formatRt(tc)) 
})
