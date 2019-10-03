## Test combineSpectra on real data sets.
library(MSnbase)
library(doParallel)
library(testthat)
MZML_PROFILE <- "~/data/2017/2017_10/profile/"
MZML_CENTR <- "~/data/2017/2017_10/centroid/"
registerDoParallel(3)
register(DoparParam(), default = TRUE)
fls <- c("20171020_POOL_POS_1.mzML",    # 2017_10
         "20171020_POOL_POS_3.mzML",
         "20171020_POOL_POS_7.mzML")

## Read the data in profile mode.
data_profile_none <- readMSData(paste0(MZML_PROFILE, fls),
                                centroided = FALSE, mode = "onDisk")
data_profile <- data_profile_none

file1 <- filterFile(data_profile, 1)
## Estimate mz resolution on the spectrum with the most peaks
file1 <- as(file1, "MSnExp")
which.max(unlist(spectrapply(file1, peaksCount)))
mzres <- estimateMzResolution(file1[[701]])
mzres_all <- estimateMzResolution(file1)
plot(density(unlist(mzres_all)))
abline(v = mzres)
abline(v = 0.0012)
## strange - have multiple mz resolutions being predicted.

mzsc_all <- estimateMzScattering(file1)
mean(unlist(mzsc_all))
plot(density(unlist(mzsc_all)))
## most have below 0.0008, some go up to 0.0012 - that's about half of mzres.
abline(v = 0.0012)

## Compare the results from combineSpectra and combineSpectra_for.
res1 <- combineSpectra(spectra(file1[700:702]))
res2 <- MSnbase:::combineSpectra_for(spectra(file1[700:702]))  ## Much slower.
all.equal(res1@intensity, res2@intensity)
all.equal(res1@mz, res2@mz)
plot(res1@intensity, file1[[701]]@intensity)
abline(0, 1)

## Now run the code that uses the estimate_mz_scattering
system.time(
    tmp <- combineSpectraMovingWindow(file1)
)
##     user   system  elapsed 
## 2447.395   91.265 2544.897 
##     user   system  elapsed 
## 1520.200   67.344 1587.714 


## Using the pre-computed values - how much can we improve here?
system.time(
    tmp2 <- combineSpectraMovingWindow(file1, mzd = min(unlist(mzres_all)) / 2)
)
##     user   system  elapsed 
## 2260.998   63.988 2329.366 
##     user   system  elapsed 
## 1510.594   59.905 1574.037 


## Check the number of peaks before and after!!!
pcount_1 <- unlist(spectrapply(file1, peaksCount))
pcount_2 <- unlist(spectrapply(tmp, peaksCount))


sps <- spectra(file1[2:4])


## DONE HERE


sp_orig <- spectra(data_profile)

sps <- split(sp_orig, f = fromFile(data_profile))

z <- sps[[1]]


sps_unsp <- unsplit(sps, f = fromFile(x))
names(sps_unsp) <- featureNames(x)

expect_equal(sps_unsp, sp_orig)

z <- split(spectra(x), fromFile(x))[[1]]
intF <- sum
mzF <- mean
hws <- 1L

## combineSpectra
x <- z[1:5]
mzs <- unlist(lapply(x, function(z) z@mz))
mz_order <- order(mzs)

dens <- density(diff(mzs[mz_order]), n = length(mzs))
idxs <- which(diff(sign(diff(dens$y))) == 2)
plot(dens)
abline(v = dens$x[idxs])

plot(dens, xlim = c(0, 0.1))
abline(v = dens$x[idxs])


mz1 <- mz(x[[1]])
mz2 <- mz(x[[2]])
mz3 <- mz(x[[3]])
plot(y = c(mz1, mz2, mz3), x = c(rep(1, length(mz1)), rep(2, length(mz2)),
                                 rep(3, length(mz3))),
     ylim = c(51, 53))

d1 <- density(diff(mz1), n = length(mz1))
plot(d1, xlim = c(0, 0.1))
d1$x[which.max(d1$y)]

d2 <- density(diff(mz2), n = length(mz2))
plot(d2, xlim = c(0, 0.5))
d2$x[which.max(d2$y)]

d3 <- density(diff(mz3), n = length(mz3))
plot(d3, xlim = c(0, 0.5))
d3$x[which.max(d3$y)]


all_mz <- sort(c(mz1, mz2, mz3))
all_d <- density(diff(all_mz), n = length(all_mz))
plot(all_d, xlim = c(0, 0.05))
abline(v = c(d1$x[which.max(d1$y)], d2$x[which.max(d2$y)], d3$x[which.max(d3$y)]))

## So, either estimate the default distance between m/z values within each sample
## and use the half of the max peak, or use the first minimum.


sps <- spectra(data_profile)
sps_f1 <- split(sps, fromFile(data_profile))[[1]]
system.time(
    mz_res <- lapply(sps_f1, function(z) {
        .estimate_mz_resolution(z@mz)
    })
)  # Takes 70 seconds!

which.max(unlist(lapply(sps_f1, peaksCount)))
.estimate_mz_resolution(sps_f1[[701]]@mz)

median(unlist(mz_res))  # not identical but comparable


## Now compare if the .estimate_mz_scattering function is somewhat similar to
## the .estimate_mz_resolution / 2.
mzs <- list(sps_f1[[700]]@mz, sps_f1[[701]]@mz, sps_f1[[702]]@mz)

d <- density(diff(sort(unlist(mzs))), n = length(unlist(mzs))/2)
plot(d, xlim = c(0, 0.01))

## Density per file
d1 <- density(diff(mzs[[1]]), n = length(mzs[[1]]))
d2 <- density(diff(mzs[[2]]), n = length(mzs[[2]]))
d3 <- density(diff(mzs[[3]]), n = length(mzs[[3]]))

plot(d1, xlim = c(0, 0.01))
points(d2, col = "red", type = "l")
points(d3, col = "green", type = "l")

abline(v = .estimate_mz_scattering(sort(unlist(mzs))))
abline(v = .estimate_mz_resolution(mzs[[1]]), col = "blue")
abline(v = .estimate_mz_resolution(mzs[[1]])/2, col = "blue")

d <- density(diff(sort(unlist(mzs))), n = max(c(512, length(sort(unlist(mzs)))/2)))


## With the test file.
f <- msdata::proteomics(full.names = TRUE,
                        pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")


prof <- readMSData(f, mode = "onDisk")
sps <- spectra(prof)

which.max(unlist(lapply(sps, peaksCount)))
d <- density(diff(sps[[135]]@mz), n = length(sps[[135]]@mz)/2)
plot(d, xlim = c(0, 0.1))

## Some testing on estimate m/z scattering and m/z resolution.
od1 <- filterFile(data_profile, 1)

mzs <- mz(od1[[3]])
mzdiffs <- diff(mzs)
plot(density(mzdiffs, n = length(mzdiffs) / 4), xlim = c(0, 0.01))

mzsqrdiffs <- diff(sqrt(mzs))
plot(density(mzsqrdiffs, n = length(mzdiffs) / 4), xlim = c(0, 0.001))

mzs1 <- mz(od1[[3]])
mzs2 <- mz(od1[[4]])
mzs3 <- mz(od1[[5]])


plot(density(diff(mzs1), n = length(mzs1)/2), xlim = c(0, 0.01))
points(density(diff(mzs2), n = length(mzs1)/2), col = "green", type = "l")
points(density(diff(mzs3), n = length(mzs1)/2), col = "blue", type = "l")

mzdiffs <- diff(sort(c(mzs1, mzs2, mzs3)))
## How about sqrt mz.
mzdiffs_sqrt <- diff(sqrt(sort(c(mzs1, mzs2, mzs3))))

par(mfrow = c(1, 2))
plot(density(mzdiffs, n = length(mzs1)/2), xlim = c(0, 0.01))
mzd <- MSnbase:::.estimate_mz_scattering(mzdiffs, is_diff = TRUE)
abline(v = mzd, col = "grey")
plot(density(mzdiffs_sqrt, n = length(mzs1) / 2), xlim = c(0, 0.001))
mzd_sqrt <- MSnbase:::.estimate_mz_scattering(mzdiffs_sqrt, is_diff = TRUE)
abline(v = mzd_sqrt, col = "grey")

mzd
mzd_sqrt
sqrt(mzd)

## Define whether we should use sqrt or just mz.
qnts <- quantile(mzs1, prob = c(0, 0.1, 0.9, 1))
mzs_low <- mzs1[mzs1 >= qnts[1] & mzs1 <= qnts[2]]
mzs_high <- mzs1[mzs1 >= qnts[3] & mzs1 <= qnts[4]]

par(mfrow = c(1, 2))
plot(density(diff(mzs_low), n = length(mzs_low)), xlim = c(0, 0.01))
points(density(diff(mzs_high), n = length(mzs_high)), col = "red", type = "l")
## The same for the sqrt data.
plot(density(diff(sqrt(mzs_low)), n = length(mzs_low)), xlim = c(0, 0.0001))
points(density(diff(sqrt(mzs_high)), n = length(mzs_high)), col = "red",
       type = "l")

a <- MSnbase:::.estimate_mz_resolution(mzs_low)
b <- MSnbase:::.estimate_mz_resolution(mzs_high)

diff(c(a, b)) < 1e-6


a <- MSnbase:::.estimate_mz_resolution(sqrt(mzs_low))
b <- MSnbase:::.estimate_mz_resolution(sqrt(mzs_high))
diff(c(a, b)) < 1e-6

