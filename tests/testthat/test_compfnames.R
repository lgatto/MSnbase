context("compfnames")

## Class FeatComp, function compfnames and helper functions .compStrings and 
## .calcCompNumbers

test_that("FeatComp", {
    expect_true(validObject(MSnbase:::.FeatComp()))
})

test_that("FeatComp/.compStrings", {
    ## create vectors with simple strings
    string1 <- c("a", "b", "c", "d")
    string2 <- c("c", "d", "e", "f")
    string3 <- c("f", "e")
    x12 <- MSnbase:::.compStrings(string1, string2)
    x13 <- MSnbase:::.compStrings(string1, string3)
    x23 <- MSnbase:::.compStrings(string2, string3)
    
    expect_equal(length(x12), 1)
    expect_true(x12@all, TRUE)
    expect_equal(x12@name, "all")
    
    expect_false(MSnbase:::.compStrings(string1, string2, all = FALSE)@all)
    expect_equal(MSnbase:::.compStrings(string1, string2, 
                                all = FALSE, name = "ut")@name, "ut")
    ## test common features
    expect_equal(x12@common, c("c", "d"))
    expect_equal(x13@common, character())
    expect_equal(x23@common, c("e", "f"))
    ## test unique features for argument string1
    expect_equal(x12@unique1, c("a", "b"))
    expect_equal(x23@unique1, c("c", "d"))
    expect_equal(x13@unique1, c("a", "b", "c", "d"))
    ## test unique features for argument string2
    expect_equal(x12@unique2, c("e", "f"))
    expect_equal(x23@unique2, character())
    expect_equal(x13@unique2, c("f", "e"))    
})

test_that("FeatComp/.calcCompNumbers", {
    string1 <- c("a", "b", "c", "d")
    string2 <- c("c", "d", "e", "f")
    string3 <- c("f", "e")
    ## test .calcCompNumbers which creates matrix 
    ## create objects for input for function .calcCompNumbers
    flist11 <- MSnbase:::.compStrings(string1, string1, FALSE, "ut11")
    flist12 <- MSnbase:::.compStrings(string1, string2, FALSE, "ut12")
    flist13 <- MSnbase:::.compStrings(string1, string3, FALSE, "ut13")
    flist23 <- MSnbase:::.compStrings(string2, string3, FALSE, "ut23")
    flist <- list(flist11, flist12, flist13, flist23)
    res <- MSnbase:::.calcCompNumbers(flist)
    ## start testing
    expect_true(is.matrix(res))
    expect_equal(dim(res), c(4, 3))    
    expect_equal(rownames(res), c("ut11", "ut12", "ut13", "ut23"))
    expect_equal(colnames(res), c("common", "unique1", "unique2"))    
    expect_equal(as.numeric(res[1, ]), c(4, 0, 0))    
    expect_equal(as.numeric(res[2, ]), c(2, 2, 2))
    expect_equal(as.numeric(res[3, ]), c(0, 4, 2))
    expect_equal(as.numeric(res[4, ]), c(2, 2, 0))
})

test_that("FeatComp/compfnames", {
    ## test function compfnames
    ## create three test MSnSets
    obj1 <- new("MSnSet",
                exprs = matrix(1, ncol = 4, nrow = 4, 
                    dimnames = list(letters[1:4])),
                featureData = new("AnnotatedDataFrame",
                    data = data.frame(markers = LETTERS[1:4],
                        row.names = letters[1:4])))
    obj2 <- new("MSnSet", 
                exprs = matrix(1, ncol = 4, nrow = 4, 
                    dimnames = list(letters[3:6])),
                featureData = new("AnnotatedDataFrame",
                    data = data.frame(markers = LETTERS[3:6],
                        row.names = letters[3:6])))
    obj3 <- new("MSnSet", 
                exprs = matrix(1, ncol = 4, nrow = 4, 
                    dimnames = list(letters[5:8])),
                featureData = new("AnnotatedDataFrame",
                    data = data.frame(markers = rep(LETTERS[6:7],2),
                        row.names = letters[5:8])))    
    expect_true(validObject(obj1))
    expect_true(validObject(obj2))
    expect_true(validObject(obj3))
    
    ## start testing function compfnames
    expect_error(compfnames(obj1, obj2, "pd.markers", "markers"))
    expect_error(compfnames(obj1, obj2, "markers", "pd.markers"))
    ## create object which compares featureNames of obj1 and obj2
    comp12 <- compfnames(obj1, obj2, NULL, "markers", FALSE)
    expect_equal(length(comp12), 1)
    expect_equal(comp12[[1]]@common, c("c", "d"))
    expect_equal(comp12[[1]]@unique1, c("a", "b"))
    expect_equal(comp12[[1]]@unique2, c("e", "f"))
    expect_true(comp12[[1]]@all)
    comp12 <- compfnames(obj1, obj2, "markers", verbose = FALSE)
    expect_equal(comp12[[4]]@common, "c")
    expect_equal(comp12[[4]]@unique1, character())
    expect_equal(comp12[[4]]@unique2, character())
    expect_equal(comp12[[5]]@common, "d")
    ## create object which compares featureNames of obj1 and obj3
    comp13 <- compfnames(obj1, obj3, "markers", "markers", FALSE)
    expect_equal(length(comp13), 7)
    expect_equal(comp13[[1]]@common, character())
    expect_equal(comp13[[1]]@unique1, c("a", "b", "c", "d"))
    expect_equal(comp13[[1]]@unique2, c("e", "f", "g", "h"))
    expect_equal(comp13[[2]]@common, character())
    expect_equal(comp13[[2]]@unique1, "a")
    expect_equal(comp13[[2]]@unique2, character())
    expect_equal(comp13[[7]]@common, character())
    expect_equal(comp13[[7]]@unique1, character())
    expect_equal(comp13[[7]]@unique2, c("f", "h"))
    ## create object whcih compares featureNames of obj2 and obj3
    comp23 <- compfnames(obj2, obj3, "markers", "markers", FALSE)
    expect_equal(length(comp23), 6)
    expect_equal(comp23[[1]]@common, c("e", "f"))
    expect_equal(comp23[[1]]@unique1, c("c", "d"))
    expect_equal(comp23[[1]]@unique2, c("g", "h"))
    expect_equal(comp23[[2]]@common, character())
    expect_equal(comp23[[2]]@unique1, "c")
    expect_equal(comp23[[2]]@unique2, character())
    expect_equal(comp23[[6]]@common, character())
    expect_equal(comp23[[6]]@unique1, character())
    expect_equal(comp23[[6]]@unique2, c("f", "h"))                        
})

