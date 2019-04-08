##' Feature-based weighting of peptides for protein ratio estimation
##' (called by iPQF main function))
##'
##' @param pos list of proteins with corresponding peptide spectra
##' position assignment (in object)
##' @param ma ratio matrix/intensity matrix
##' @param features data.frame of peptide features
##' @return A matrix with estimated protein ratios (rows are proteins,
##' columns are samples).
##' @author Martina Fisher
##' @noRd
iPQF.method  <- function(pos, mat, features, feature.weight) {
    ## calculated protein ratios based on feature ranking
    feat.trend <- matrix(data = NA_real_, nrow = length(pos),
                         ncol = dim(mat)[2], byrow = FALSE)
    rownames(feat.trend) <- names(pos)
    #weight.list <- vector("list", length(pos))   ## peptide weights
    #avrank.list <- vector("list", length(pos))   ## peptide feature average ranks
    #rankmat.list <- vector("list", length(pos))  ## feature rank matrix

    for (i in 1:length(pos)) {
        Mat <- mat[pos[[i]], ]
        ## Ranking: smaller rank is better (~more reliable peptide)
        ## combined redundant-unique-dist vector
        ru.pep <- features$ru.dist[pos[[i]]]
        ru.rank <- rank(ru.pep)
        charge.rank <- rank(features$charge[pos[[i]]])
        ll.rank <- rank(features$seq.l[pos[[i]]])
        sc.rank <- rank(-features$score[pos[[i]]])
        mo.rank <- rank(features$mod.stat[pos[[i]]])
        mass.rank <- rank(features$prec.mass[pos[[i]]])
        int.rank <- rank(-features$mean.ionInt[pos[[i]]])

        ## ranking matrix
        rank.mat <- cbind(ru.rank, charge.rank,
                          ll.rank, sc.rank,
                          mo.rank, mass.rank,
                          int.rank)
        #rankmat.list[[i]] <- rank.mat

        ## feature weighting  (based on correlation study in manuscript)
        rank.mat2 <- (t(apply(rank.mat, 1,
                              function(x) feature.weight * x)))

        ## normalizing: divide by peptide number (rank zw 0-1)
        rmat.n <- apply(rank.mat2, 2, function(x) x/length(which(!is.na(x))))

        ## sum ranks for each peptide:
        pep.sumrank <- rowSums(rmat.n, na.rm = TRUE)

        ## divide by number of features: "average rank"
        av.rank <- pep.sumrank/sum(c(1:7)^2)     ## worst av.rank = 1
        av.rank <- 1 - av.rank  ## reverse: best = 1 "weight" = high = reliable peptide
        av.rank <- av.rank^2
        #avrank.list[[i]] <- av.rank

        ## Approach: Feature-Weighting
        weight <- av.rank
        weight <- weight / sum(weight)
        trend <- apply(Mat, 2, function(x) weighted.mean(x, w=weight))

        feat.trend[i,] <- trend
        #weight.list[[i]] <- weight

    }
    return(feat.trend)
}


#####  Internal - FUNCTIONS    < called within iPQF main function !>
## to build internal objects required for iPQF.method

## Ratio Matrix Construction Function
# Define Ratio Calculation: all Channels to individual channel, or sum of all channels
ratio.mat <- function(mat, method) {
    if (method == "sum") {
        ## equivalent to normalise(., "sum")
        ## ration.mat <- mat/rowSums(mat)
        ratio.mat <- t(apply(mat, 1, function(x) x/sum(x)))
    } else {
        base.channel <- match(method, colnames(mat))
        ## ratio.mat <- mat/mat[, base.channel]
        ratio.mat <- t(apply(mat, 1,
                             function(x) x /x[base.channel]))
  }
  return(ratio.mat)
}


## example call:
## exprs(object) <- ratio.mat(exprs(object), method="sum")
## ratio.mat(head(mat), method="X114_ions")

## Function: 'Build uniques.all'
## list elements= proteins
## entries per list.element: discrete numbers presenting spectra
## assigned to the specific protein with same numbers for identical
## sequences
uniques.list <- function(pos, sequence) {
    pep.all <- lapply(pos, function(i) sequence[i])
    un.pep.all <- lapply(pos, function(i) unique(sequence[i]))
    uniques.all <- lapply(1:length(pos),
                         function(x) match(pep.all[[x]], un.pep.all[[x]]))
    names(uniques.all) <- names(pos)
    return(uniques.all)
}



## Function: 'Redundant peptides' (multiply measured sequence):
## pos.r: position of peptides with redundantly measured status
## pos.rd: Mean Distance of each peptide to other peptides within a redundant group  (follows pos.r structure)
redundant.dist <- function(pos, uniques.all, mat) {
    uni.tab <- lapply(uniques.all, table)    ## list: name=individual sequence, value= ## appearance of same sequence
    anz.uni <- lapply(uni.tab,length)        ## number of different sequences in a protein profile

    red.name <- lapply(uni.tab, function(x) which(x>1))        ## which protein has redundant sequences?
    red.prot <- which(lapply(red.name, length)>0)              ## list index (protein profiles) with redundant sequences

    ## redundant group of peptides:
    pos.r <- vector("list", length(pos))                ## pos.r: position of redundant peptide spectra
    pos.rd <- vector("list", length(pos))           ## distance of redundant peptides
    for (i in red.prot) {
        nn <- names(red.name[[i]])
        a <- dist <- c()
        for (k in nn) {
            posr <- which(uniques.all[[i]] == k)
            a <- c(a, pos[[i]][posr])
            Mat <- mat[pos[[i]][posr], ]
            red.vec <- vector("numeric", length(posr))
            for (j in 1:length(posr)){
                if (length(posr) == 2)
                    red.vec[j] <-  sqrt(sum((Mat[-j,]-Mat[j,])^2))
                if (length(posr) > 2)
                    red.vec[j] <- mean(apply(Mat[-j,], 1,
                                             function(x) sqrt(sum((x-Mat[j,])^2))))
            }
            dist <- c(dist, red.vec)
        }
        pos.r[[i]] <- a
        pos.rd[[i]] <- dist
    }
    return(list(pos.r, pos.rd))
}



## Function: Unique peptide spectrum match (sequence only measured once)
## pos.u: position of peptides with uniquely measured status
## pos.ud: mean distance of a 'unique' peptide to other 'unique' sequences assigned to the protein (follows pos.u structure)

uni.measured.dist <- function(pos, uniques.all, mat) {
    uni.tab <- lapply(uniques.all, table)
    uni.name <- lapply(uni.tab, function(x) which(x==1)) ## which protein has single measured peptides?
    uni.prot <- which(lapply(uni.name, length) > 0)      ## index protein profiles with single measured Peptides

    num.unis <- unlist(lapply(uni.name, length))  ## number single uniques for each protein
    table(num.unis)

    pos.u <- vector("list", length(pos))
    pos.ud <- vector("list", length(pos))

    for (i in uni.prot) {
        uni.p <- match(names(uni.name[[i]]), uniques.all[[i]])
        pos.u[[i]] <- pos[[i]][uni.p]
        Mat <- mat[pos.u[[i]], ]
        uni.vec <- vector("numeric", length(pos.u[[i]]))
        for (j in 1:length(uni.vec)){
            if (length(uni.vec) == 1) {
                Mat <- mat[pos[[i]], ]
                k <- match(pos.u[[i]], pos[[i]]) ## for singles: mean distance to all other peptides of the protein
                uni.vec[j] <- mean(apply(Mat[-k,], 1,
                                         function(x) sqrt(sum((x-Mat[k,])^2))))
            }
            if (length(uni.vec) == 2)
                uni.vec[j] <-  sqrt(sum((Mat[-j,]-Mat[j,])^2))
            if (length(uni.vec) > 2)
                uni.vec[j] <- mean(apply(Mat[-j,], 1, function(x) sqrt(sum((x-Mat[j,])^2))))
        }
        pos.ud[[i]] <- uni.vec
    }
    return(list(pos.u, pos.ud))
}


##' The iPQF spectra-to-protein summarisation method integrates
##' peptide spectra characteristics and quantitative values for protein
##' quantitation estimation. Spectra features, such as charge state,
##' sequence length, identification score and others, contain valuable
##' information concerning quantification accuracy. The iPQF algorithm
##' assigns weights to spectra according to their overall feature reliability
##' and computes a weighted mean to estimate protein quantities.
##' See also \code{\link{combineFeatures}} for a more
##' general overview of feature aggregation and examples.
##'
##' @title iPQF: iTRAQ (and TMT) Protein Quantification based on
##'     Features
##' @param object An instance of class \code{MSnSet} containing
##'     absolute ion intensities.
##' @param groupBy Vector defining spectra to protein
##'     matching. Generally, this is a feature variable such as
##'     \code{fData(object)$accession}.
##' @param low.support.filter A \code{logical} specifying if proteins
##'     being supported by only 1-2 peptide spectra should be filtered
##'     out. Default is \code{FALSE}.
##' @param ratio.calc Either \code{"none"} (don't calculate any
##'     ratios), \code{"sum"} (default), or a specific channel (one of
##'     \code{sampleNames(object)}) defining how to calculate relative
##'     peptides intensities.
##' @param feature.weight Vector \code{"numeric"} giving weight to the
##'     different features. Default is the squared order of the
##'     features redundant -unique-distance metric, charge state, ion
##'     intensity, sequence length, identification score, modification
##'     state, and mass based on a robustness analysis.
##' @param method.combine A \code{logical} defining whether to further
##'     use median polish to combine features.
##' @return A \code{matrix} with estimated protein ratios.
##' @author Martina Fischer
##' @references iPQF: a new peptide-to-protein summarization method
##'     using peptide spectra characteristics to improve protein
##'     quantification. Fischer M, Renard BY.  Bioinformatics. 2016
##'     Apr 1;32(7):1040-7. doi:10.1093/bioinformatics/btv675. Epub
##'     2015 Nov 20. PubMed PMID:26589272.
##' @examples
##' data(msnset2)
##' head(exprs(msnset2))
##' prot <- combineFeatures(msnset2,
##'                         groupBy = fData(msnset2)$accession,
##'                         method = "iPQF")
##' head(exprs(prot))
iPQF <- function(object, groupBy,
                 low.support.filter = FALSE,
                 ratio.calc = "sum",
                 method.combine = FALSE,
                 feature.weight = c(7,6,4,3,2,1,5)^2
                 ) {

    if (!inherits(object,"MSnSet"))
        stop("'object' is required to be of class MSnSet")
    ## Check NA/Zero values still in data set?
    rm.pos <- apply(exprs(object), 2,
                    function(x) which(is.na(x) | x==0))
    rm.rows <- unique(unlist(rm.pos))
    if (length(rm.rows) > 0)
        stop("Remove NA/Zero Intensities in you data",
             "before peptide summarization. ",
             length(rm.rows),
             " spectr[a/um] should be removed.")

    ## Check mzTab standard names are provdied?
    mzTab.names <- c("sequence", "accession",
                     "charge", "modifications",
                     "mass_to_charge",
                     "search_engine_score")
    wrong.colnames <- which( mzTab.names %in% colnames(fData(object)) == FALSE)
    if (length(wrong.colnames) > 0) {
        stop(" In FeatureData the following column names, according to the mzTab standard, are required: ",
             "\n", paste(mzTab.names[wrong.colnames], collapse = ", ")) }

    ## Extract individual features
    sequence <- as.character(fData(object)$sequence)
    accession <- as.character(fData(object)$accession)

    charge <- as.integer(as.vector(fData(object)$charge))
    seq.l  <- sapply(sequence, nchar)
    prec.mass   <- as.numeric(as.vector(fData(object)$mass_to_charge)) * charge

    ses.na <- which(is.na(fData(object)$search_engine_score))
    score  <- as.numeric(as.vector(fData(object)$search_engine_score))
    if (length(which(is.na(score))) > length(ses.na))
        stop("Search engine scores are expected as column of numeric values only.")

    mod.stat.org <- as.character(fData(object)$modifications)
    mod.stat <- ifelse(mod.stat.org == "null" | mod.stat.org == "0" , 0, 1) ## FIXME

    ## Should the above not be NULL and 0? These are probably badly
    ## imported values from mzTab.

    ion.mat <- exprs(object)               ## absolute ion intensity
    mean.ionInt <- apply(ion.mat, 1, mean) ## mean absolute ion intensity

    ## Calculate Ratio Matrix
    if (ratio.calc != "none") {
        mat <- ratio.mat(exprs(object), method=ratio.calc)  ## relative intensities
    }
    else mat <- ion.mat  ## absolute intensities

    ## Protein-Peptide Spectra position assignment
    ## Build list 'pos.pep': Protein ID (accession name)- Peptide spectra assignment (position in object)

    ## group.by character vector! groupBy = feature pep reihenfolge
    uni.ids <- levels(as.factor(groupBy))
    pos.all <- sapply(uni.ids, function(y) which(groupBy==y))
    names(pos.all) <- uni.ids

    ## remove low-supported proteins
    singles <- which(unlist(lapply(pos.all, length)) < 3) ## proteins supported by only 1-2 peptide spectra

    if (length(singles) == length(pos.all)) {
        warning("No proteins with > 2 peptide spectra - keeping everything.")
        pos.pep <- pos.all
    } else {
        pos.pep <- pos.all[-singles] ## proteins supported by >2 peptide spectra
    }

    ## Redundantly measured peptide spectrum status in a protein profile (uniques.all):
    uniques.all <- uniques.list(pos.pep, sequence)

    ## Redundant peptides
    red.result <- redundant.dist(pos.pep, uniques.all, mat)
    pos.r <- red.result[[1]]
    pos.rd <- red.result[[2]]
    ## Single measured peptides
    uni.result <- uni.measured.dist(pos.pep, uniques.all, mat)
    pos.u <- uni.result[[1]]
    pos.ud <- uni.result[[2]]

    ## Combined distance vector:
    ru.dist <- rep(NA, nrow(object))
    ru.dist[unlist(pos.r)] <- unlist(pos.rd)
    ru.dist[unlist(pos.u)] <- unlist(pos.ud)


    ## Data.frame of selected peptide features:
    features <- data.frame(charge, seq.l, prec.mass,
                           score, mod.stat, mean.ionInt, ru.dist)


    ## Protein Quantification - Peptide Summarization
    iPQF.result  <- iPQF.method(pos.pep, mat, features, feature.weight=feature.weight)

    if (!low.support.filter) {
        single.prots <- unique(accession[unlist(pos.all[singles])])
        msg <- paste0("The following ", length(single.prots), " proteins are only supported by 1 or 2 peptide spectra,\n",
                      "hence, protein quantification is not reliable and can only be calculated\n",
                      "by the 'mean' in these cases, corresponding protein accessions are:\n  ",
                      paste(single.prots, collapse = ", "))
        message(msg)
        single.quant <- lapply(pos.all[singles],
                               function(i) {
                                   if (length(i) == 2) apply(mat[i,],2,mean) else mat[i,]
                               })
        single.quant <- do.call(rbind, single.quant)
        quant.result <- rbind(iPQF.result, single.quant)
    } else quant.result <- iPQF.result

    quant.result <- quant.result[match(uni.ids, rownames(quant.result)), ]
    ## Method combination: iPQF with MedianPolish
    if (method.combine) {
        MP.quant <- lapply(1:length(pos.all),
                          function(i) {
                              medpol <- medpolish(mat[pos.all[[i]],],
                                                  trace.iter = FALSE)
                              if (length(pos.all[[i]]) == 1)
                                  result <- medpol$overall + medpol$row
                              else
                                  result <- medpol$overall + medpol$col
                              return(result)
                          } )
        MP.quant <- do.call(rbind, MP.quant)
        quant.result <- t(sapply(1:length(pos.all),
                                 function(k) apply(rbind(MP.quant[k,], quant.result[k,]), 2, mean)))
        rownames(quant.result) <- names(pos.all)
    }
    return(quant.result)
}
