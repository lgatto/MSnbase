## Some code lines to evaluate performance of accessing SQLite MS data.

library(MSnbase)
library(msdata)
library(mzR)
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))
inMem <- readMSData(files = mzf, msLevel. = 1, centroided. = TRUE)
onDisk <- readMSData2(files = mzf, msLevel. = 1, centroided. = TRUE)


file_con <- mzR::openMSfile(mzf[1])
pks <- peaks(file_con)
hdr <- header(file_con)

data <- do.call(rbind, pks)
data <- cbind(rep(1:length(pks), lengths(pks)/2), data)
colnames(data) <- c("acquisition_num", "mz", "intensity")

## Create the SQL database
library(RSQLite)
db_con <- dbConnect(SQLite(), dbname = "test.sqlite")
dbWriteTable(db_con, name = "mzint", data.frame(data), row.names = FALSE)
dbGetQuery(db_con, "create index mzint_acqnum_idx on mzint (acquisition_num);")
dbDisconnect(db_con)

db_con <- dbConnect(SQLite(), dbname = "test.sqlite")

library(microbenchmark)
library(testthat)
a <- peaks(file_con)
a <- data.frame(do.call(rbind, a))
b <- dbGetQuery(db_con, "select mz, intensity from mzint;")

expect_equal(unname(a), unname(b))

## The full data file
microbenchmark(do.call(rbind, peaks(file_con)),
               dbGetQuery(db_con, "select mz, intensity from mzint;"),
               times = 5)
## Unit: milliseconds
##                                                    expr       min        lq
##                         do.call(rbind, peaks(file_con)) 119.24763 120.53678
##  dbGetQuery(db_con, "select mz, intensity from mzint;")  87.99631  89.53832
##      mean    median        uq      max neval cld
##  122.4439 122.64996 124.73163 125.0535     5   b
##   91.4046  90.48934  92.60442  96.3946     5  a


## random samples:
set.seed(18011977)
idxs <- sort(sample(1:length(pks), 20))
microbenchmark(do.call(rbind, peaks(file_con, idxs)),
               dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",
                                         "acquisition_num in (",
                                        paste0(idxs, collapse = ", "),");")),
               times = 5)
## Unit: milliseconds
##                                                                                                                                    expr
##                                                                                                   do.call(rbind, peaks(file_con, idxs))
##  dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num in (", paste0(idxs, collapse = ", "), ");"))
##       min       lq     mean   median       uq      max neval cld
##  20.30434 20.76485 20.70102 20.78625 20.82145 20.82823     5  a
##  24.28117 24.48217 24.90876 24.73960 24.96314 26.07773     5   b

## A single sample:
microbenchmark(a <- peaks(file_con, 13),
               b <- dbGetQuery(db_con,
                               paste0("select mz, intensity from mzint where ",
                                      "acquisition_num = 13;")),
               times = 5)
## Unit: microseconds
##                                                                                                     expr
##                                                                                 a <- peaks(file_con, 13)
##  b <- dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num = 13;"))
##       min       lq     mean   median       uq      max neval cld
##   141.080  152.425  388.716  165.853  200.209 1284.013     5  a
##  1748.661 1822.920 1911.104 1882.123 1923.760 2178.055     5   b

## All spectra: RSQLite faster.
## limited spectra: mzR faster.
## single spectrum: mzR faster.


dbDisconnect(db_con)
mzR::close(file_con)
## Seems getting the full data is better in SQLite, but extracting subsets is
## considerable faster with mzR.

## Make that an issue in MSnbase.


############################################################
#### Large data file
mz_f <- "/Volumes/Ext64/data/2016/2016-11/NoSN/190516_BLANK_POS_1.mzML"

file_con <- mzR::openMSfile(mz_f)
pks <- peaks(file_con)
hdr <- header(file_con)

data <- do.call(rbind, pks)
data <- cbind(rep(1:length(pks), lengths(pks)/2), data)
colnames(data) <- c("acquisition_num", "mz", "intensity")

## Store the data into the SQLite database.
library(RSQLite)
db_con <- dbConnect(SQLite(), dbname = "large_file.sqlite")
dbWriteTable(db_con, name = "mzint", data.frame(data), row.names = FALSE)
dbGetQuery(db_con, "create index mzint_acqnum_idx on mzint (acquisition_num);")
dbDisconnect(db_con)


## Just connect to the pre-created file
db_con <- dbConnect(SQLite(), dbname = "large_file.sqlite")
## Note: the original mzML is ~ 45M
##       the sqlite database: ~ 129M

## The full data set:
microbenchmark(do.call(rbind, peaks(file_con)),
               dbGetQuery(db_con, "select mz, intensity from mzint;"),
               times = 5)
## Unit: seconds
##                                                    expr      min       lq
##                         do.call(rbind, peaks(file_con)) 3.849484 3.855993
##  dbGetQuery(db_con, "select mz, intensity from mzint;") 2.173048 2.186558
##      mean   median       uq      max neval cld
##  3.874310 3.867167 3.883957 3.914946     5   b
##  2.206608 2.204896 2.221304 2.247237     5  a

## random samples:
set.seed(18011977)
idxs <- sort(sample(1:length(pks), 20))
microbenchmark(do.call(rbind, peaks(file_con, idxs)),
               dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",
                                         "acquisition_num in (",
                                         paste0(idxs, collapse = ", "),");")),
               times = 5)
## Unit: milliseconds
##                                                                                                                                    expr
##                                                                                                   do.call(rbind, peaks(file_con, idxs))
##  dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num in (", paste0(idxs, collapse = ", "), ");"))
##       min       lq     mean   median       uq      max neval cld
##  46.80579 47.60341 48.26999 47.98958 48.20668 50.74451     5   a
##  44.76770 47.62399 49.21587 50.45550 50.93661 52.29557     5   a

## A single Spectrum:
microbenchmark(a <- peaks(file_con, 13),
               b <- dbGetQuery(db_con,
                               paste0("select mz, intensity from mzint where ",
                                      "acquisition_num = 13;")),
               times = 5)
## Unit: microseconds
##                                                                                                     expr
##                                                                                 a <- peaks(file_con, 13)
##  b <- dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num = 13;"))
##       min       lq     mean   median       uq      max neval cld
##   117.159  130.465  398.175  142.672  172.915 1427.664     5  a
##  1591.109 1599.142 1706.570 1655.031 1668.107 2019.460     5   b


## So:
## 1) Full data: RSQlite faster.
## 2) Limited number of samples: same performance.
## 3) Single spectrum: mzR faster.


dbDisconnect(db_con)
mzR::close(file_con)

############################################################
#### gzipped file with uncompressed binary data
mz_f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")

file_con <- mzR::openMSfile(mz_f)
pks <- peaks(file_con)
hdr <- header(file_con)

data <- do.call(rbind, pks)
data <- cbind(rep(1:length(pks), lengths(pks)/2), data)
colnames(data) <- c("acquisition_num", "mz", "intensity")

## Create the SQL database
db_con <- dbConnect(SQLite(), dbname = "prote_1.sqlite")
dbWriteTable(db_con, name = "mzint", data.frame(data), row.names = FALSE)
dbGetQuery(db_con, "create index mzint_acqnum_idx on mzint (acquisition_num);")
dbDisconnect(db_con)

db_con <- dbConnect(SQLite(), dbname = "prote_1.sqlite")

## Size: mzML.gz: 16MB
##       RSQLite: 63MB

## The full data set:
microbenchmark(do.call(rbind, peaks(file_con)),
               dbGetQuery(db_con, "select mz, intensity from mzint;"),
               times = 5)
## Unit: seconds
##                                                    expr     min       lq
##                         do.call(rbind, peaks(file_con)) 1.26702 1.269837
##  dbGetQuery(db_con, "select mz, intensity from mzint;") 1.03514 1.035995
##      mean   median       uq      max neval cld
##  1.282256 1.277391 1.281634 1.315399     5   b
##  1.070749 1.091303 1.094308 1.096998     5  a


## random samples:
set.seed(18011977)
idxs <- sort(sample(1:length(pks), 20))
microbenchmark(do.call(rbind, peaks(file_con, idxs)),
               dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",
                                         "acquisition_num in (",
                                         paste0(idxs, collapse = ", "),");")),
               times = 5)
## Unit: milliseconds
##                                                                                                                                    expr
##                                                                                                   do.call(rbind, peaks(file_con, idxs))
##  dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num in (", paste0(idxs, collapse = ", "), ");"))
##        min        lq      mean    median        uq       max neval cld
##  74.252536 75.527560 76.482204 76.565536 77.697335 78.368055     5   b
##   5.155087  5.281799  5.402167  5.346787  5.528892  5.698268     5  a


## A single Spectrum:
microbenchmark(a <- peaks(file_con, 13),
               b <- dbGetQuery(db_con,
                               paste0("select mz, intensity from mzint where ",
                                      "acquisition_num = 13;")),
               times = 5)
## Unit: microseconds
##                                                                                                     expr
##                                                                                 a <- peaks(file_con, 13)
##  b <- dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num = 13;"))
##       min       lq     mean   median       uq      max neval cld
##   168.211  175.492 1335.600  203.130  232.770 5898.396     5   a
##  2213.052 2229.595 2310.885 2243.433 2248.521 2619.822     5   a

## Same as before.
dbDisconnect(db_con)
mzR::close(file_con)


############################################################
#### unzipped file with uncompressed binary data
mz_f <- "/Users/jo/tmp/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"

file_con <- mzR::openMSfile(mz_f)
pks <- peaks(file_con)
hdr <- header(file_con)

data <- do.call(rbind, pks)
data <- cbind(rep(1:length(pks), lengths(pks)/2), data)
colnames(data) <- c("acquisition_num", "mz", "intensity")

## Create the SQL database
db_con <- dbConnect(SQLite(), dbname = "prote_2.sqlite")
dbWriteTable(db_con, name = "mzint", data.frame(data), row.names = FALSE)
dbGetQuery(db_con, "create index mzint_acqnum_idx on mzint (acquisition_num);")
dbDisconnect(db_con)

db_con <- dbConnect(SQLite(), dbname = "prote_2.sqlite")

## Size: mzML.gz: 31MB
##       RSQLite: 63MB

## The full data set:
microbenchmark(do.call(rbind, peaks(file_con)),
               dbGetQuery(db_con, "select mz, intensity from mzint;"),
               times = 5)
## Unit: seconds
##                                                    expr      min       lq
##                         do.call(rbind, peaks(file_con)) 1.066969 1.067947
##  dbGetQuery(db_con, "select mz, intensity from mzint;") 1.043946 1.062582
##      mean   median       uq      max neval cld
##  1.080404 1.071158 1.072499 1.123449     5   a
##  1.078063 1.082835 1.092297 1.108655     5   a

## random samples:
set.seed(18011977)
idxs <- sort(sample(1:length(pks), 20))
microbenchmark(do.call(rbind, peaks(file_con, idxs)),
               dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",
                                         "acquisition_num in (",
                                         paste0(idxs, collapse = ", "),");")),
               times = 5)
## Unit: milliseconds
##                                                                                                                                    expr
##                                                                                                   do.call(rbind, peaks(file_con, idxs))
##  dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num in (", paste0(idxs, collapse = ", "), ");"))
##       min       lq      mean   median        uq       max neval cld
##  9.417424 9.507414 10.068734 9.835241 10.325901 11.257690     5   b
##  5.114088 5.125851  5.323403 5.313672  5.408287  5.655118     5  a

## A single Spectrum:
microbenchmark(a <- peaks(file_con, 13),
               b <- dbGetQuery(db_con,
                               paste0("select mz, intensity from mzint where ",
                                      "acquisition_num = 13;")),
               times = 5)
## Unit: microseconds
##                                                                                                     expr
##                                                                                 a <- peaks(file_con, 13)
##  b <- dbGetQuery(db_con, paste0("select mz, intensity from mzint where ",      "acquisition_num = 13;"))
##       min       lq      mean   median       uq      max neval cld
##   169.158  177.188  432.3446  205.424  243.988 1365.965     5  a
## 2262.097 2295.035 2372.3092 2295.061 2301.724 2707.629     5   b

## Interestingly, same performance for full data access.

dbDisconnect(db_con)
mzR::close(file_con)


## MAKE THAT AN ISSUE!!!
