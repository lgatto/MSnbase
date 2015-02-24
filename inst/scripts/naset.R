library("MSnbase")
library("pRolocdata")
data(dunkley2006, package = "pRolocdata")

x <- dunkley2006
set.seed(1)
i <- sample(prod(dim(x)), 500)
exprs(x)[i] <- NA
j <- sample(nrow(x), 15)
k <- sample(nrow(x), 20)
exprs(x)[j, 1:8] <- NA
exprs(x)[k, 9:16] <- NA

randna <- rep(TRUE, nrow(x))
randna[c(j,k)] <- FALSE

fd <- data.frame(nNA = apply(exprs(x), 1,
                     function(xx) sum(is.na(xx))),
                 randna = randna)
rownames(fd) <- featureNames(x)

pd <- data.frame(nNA = apply(exprs(x), 2,
                     function(xx) sum(is.na(xx))))
rownames(pd) <- sampleNames(x)


naset <- MSnSet(exprs = exprs(x),
                fd, pd)

save(naset, file = "../../data/naset.rda")
