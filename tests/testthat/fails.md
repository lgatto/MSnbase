# `filterAcquisitionNum(ondisk, n = 1:100)` (`@test_MSnExpFilters.R#78`) 2016-08-01

Notes: The `initialize` calls in 11 and 12 are from
`Biobase::"[",AnnotatedDataFrame` (called via `NAnnotatedDataFrame`, which extends it).

```r
> library(testthat)
> test_package("MSnbase", filter = "MSnExp")[OnDisk]MSnExp equality: .......
MSnExp class: ....................................................
MSnExp processing:
MSnExp data: ..............
MSnExp identification data: .........................
MSnExp filter functions: ......................1Error in signalCondition(e) :
  no function to return from, jumping to top level
  ..................2.....
  OnDiskMSnExp class, other methods: ^C

Failed -------------------------------------------------------------------------
1. Error: filterAcquisitionNum (@test_MSnExpFilters.R#78) ----------------------
attempt to apply non-function
1: expect_warning(ll <- length(filterAcquisitionNum(ondisk, n = 1:100))) at /home/lg390/R/x86_64-pc-linux-gnu-library/3.3/MSnbase/tests/testthat/test_MSnExpFilters.R:78
2: capture_warnings(object)
3: withCallingHandlers(code, warning = function(condition) {
    out$push(condition)
    invokeRestart("muffleWarning")
})
4: filterAcquisitionNum(ondisk, n = 1:100)
5: filterAcquisitionNum(ondisk, n = 1:100)
6: .local(object, ...)
7: object[selAcqN | !selFile]
8: object[selAcqN | !selFile]
9: .local(x, i, j, ..., drop)
10: phenoData(x)[file, , drop = FALSE]
11: phenoData(x)[file, , drop = FALSE]
12: initialize(x, data = pD, varMetadata = mD)
13: initialize(x, data = pD, varMetadata = mD)
14: .local(.Object, ...)
15: callNextMethod(.Object, data = data, varMetadata = varMetadata, ...)
16: eval(call, callEnv)
17: eval(expr, envir, enclos)
18: .nextMethod(.Object, data = data, varMetadata = varMetadata, ...)
19: .local(.Object, ...)
20: callNextMethod(.Object, ...)
21: addNextMethod(method, f, envir = methodEnv)
22: addNextMethod(method, f, envir = methodEnv)
23: .findNextFromTable(method, f, optional, envir)
24: .findInheritedMethods(defined, fdef, mtable = NULL, excluded = excluded)
25: length(classes)
26: (function (x)
    x$.self$finalize())(<environment>)

2. Failure: filterMz (@test_MSnExpFilters.R#125) -------------------------------
all.equal(onDisk2F, filterMsLevel(onDiskF, msLevel. = 2)) isn't true.
```
