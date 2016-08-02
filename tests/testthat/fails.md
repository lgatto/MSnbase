## `filterAcquisitionNum(ondisk, n = 1:100)` (`@test_MSnExpFilters.R#78`) 

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
## `readMSData with pdata` (`@test_MSnExp.R#36`)

Note: again something related to `AnnotatedDataFrame`.

```r
Failed -------------------------------------------------------------------------
1. Error: readMSData with pdata (@test_MSnExp.R#36) ----------------------------
attempt to apply non-function
1: readMSData(file, pdata = pd, verbose = FALSE) at /home/lg390/R/x86_64-pc-linux-gnu-library/3.3/MSnbase/tests/testthat/test_MSnExp.R:36
2: new("AnnotatedDataFrame", data = data.frame(spectrum = 1:length(nms), row.names = nms))
3: initialize(value, ...)
4: initialize(value, ...)
5: .local(.Object, ...)
6: callNextMethod(.Object, data = data, varMetadata = varMetadata, ...)
7: eval(call, callEnv)
8: eval(expr, envir, enclos)
9: .nextMethod(.Object, data = data, varMetadata = varMetadata, ...)
10: .local(.Object, ...)
11: callNextMethod(.Object, ...)
12: addNextMethod(method, f, envir = methodEnv)
13: addNextMethod(method, f, envir = methodEnv)
14: .findNextFromTable(method, f, optional, envir)
15: .findInheritedMethods(defined, fdef, mtable = NULL, excluded = excluded)
16: is(m, "MethodDefinition")
17: getClassDef(cl)
18: .getClassFromCache(Class, where, package = package, resolve.confl = "none")
19: (function (x)
   x$.self$finalize())(<environment>)
```

## `readMSData` and dummy MSnExp msLevel 2 instance (`@test_MSnExp.R#90`)

Note: `[` subsetting and initialise `AnnotatedDataFrame`.

```r
2. Error: readMSData and dummy MSnExp msLevel 2 instance (@test_MSnExp.R#90) ---
attempt to apply non-function
1: aa[1:2] at /home/lg390/R/x86_64-pc-linux-gnu-library/3.3/MSnbase/tests/testthat/test_MSnExp.R:90
2: aa[1:2]
3: .local(x, i, j, ..., drop)
4: featureData(x)[i, ]
5: featureData(x)[i, ]
6: initialize(x, data = pD, varMetadata = mD)
7: initialize(x, data = pD, varMetadata = mD)
8: .local(.Object, ...)
9: callNextMethod(.Object, data = data, varMetadata = varMetadata, ...)
10: eval(call, callEnv)
11: eval(expr, envir, enclos)
12: .nextMethod(.Object, data = data, varMetadata = varMetadata, ...)
13: .local(.Object, ...)
14: callNextMethod(.Object, ...)
15: addNextMethod(method, f, envir = methodEnv)
16: addNextMethod(method, f, envir = methodEnv)
17: .findNextFromTable(method, f, optional, envir)
18: .findInheritedMethods(defined, fdef, mtable = NULL, excluded = excluded)
19: is(m, "MethodDefinition")
20: possibleExtends(cl, class2, class1Def, class2Def)
21: .identC(class1[[1L]], class2)
22: (function (x)
   x$.self$finalize())(<environment>)
```   