# A class for either numeric or logical
setClassUnion("numericORlogical", c("numeric", "logical"))

#' @title quantro
#'
#' @exportClass quantro 
#' 
setClass(Class = "quantro", 
         representation = representation(
             summary = "list", B = "numeric",
             anova = "data.frame", MSbetween = "numeric",
             MSwithin = "numeric", quantroStat = "numeric",
             quantroStatPerm = "numericORlogical", 
             quantroPvalPerm = "numericORlogical")
)

setMethod("show", "quantro", 
          function(object){
            cat("quantro: Test for global differences in distributions\n")
            cat("   nGroups: ", paste(object@summary$nGroups), "\n")
            cat("   nTotSamples: ", paste(object@summary$nTotSamples), "\n")
            cat("   nSamplesinGroups: ", 
                paste(object@summary$nSamplesinGroups), "\n")
            cat("   anovaPval: ", paste(round(object@anova$`Pr(>F)`[1],5)), 
                "\n") 
            cat("   quantroStat: ", paste(round(object@quantroStat,5)), "\n")
            if(is.na(object@quantroPvalPerm)){
                noPermutation <- "Use B > 0 for permutation testing."
                cat("   quantroPvalPerm: ", paste(noPermutation), "\n")
            } else {
                showPval <- ifelse(object@quantroPvalPerm < (1/object@B), 
                       paste0("< ", 1 / object@B), object@quantroPvalPerm)
                cat("   quantroPvalPerm: ", paste(showPval), "\n")  
            }
          }
)

#' @title quantro
#'
#' @description This is a function that tests for global differences between 
#' groups of distributions which asses whether global normalization 
#' methods such as quantile normalization should be applied. 
#' This function defines the quantro class and constructor. 
#'
#' @param object an object which is inherited from an 
#' \code{eSet} such as an \code{ExpressionSet} or 
#' \code{MethylSet} object. The \code{object} can also be a 
#' data frame or matrix with observations
#' (e.g. probes or genes) on the rows and samples as the 
#' columns.  
#' @param groupFactor a group level factor associated with 
#' each sample or column in the \code{object}. The order of the
#' \code{groupFactor} must match the order of the columns in
#' \code{object}. 
#' @param B number of permutations to assess statistical significance 
#' in a permutation test. Default \code{B}=0. 
#' @param qRange the range of quantiles to consider. Default is 
#' \code{seq(0, 1, length.out = nrow(object))}. 
#' @param useMedianNormalized TRUE/FALSE argument specifying if the 
#' median normalized data should be used or not as input to test for 
#' global differences between distributions. Default is TRUE. 
#' @param verbose TRUE/FALSE argument specifying if verbose messages 
#' should be returned or not. Default is TRUE. 
#' 
#' @return A \code{quantro} S4 class object
#' \item{summary}{Returns a list of three elements 
#' related to a summary of the experiment: 
#' (1) the number of groups (nGroups), 
#' (2) total number of samples (nTotSamples), 
#' (3) number of samples in each group (nSamplesinGroups).}
#' \item{B}{Number of permutations for permutation testing.}
#' \item{anova}{ANOVA to test if the medians of the 
#' distributions (averaged across groups) are different across groups.}
#' \item{quantroStat}{A test statistic which is a ratio of the mean squared 
#' error between groups of distributions to the mean squared error within 
#' groups of distributions (psuedo F-statistic).}
#' \item{quantroStatPerm}{If \code{B} is not equal to 0, then a permutation 
#' test was performed to assess the statistical significance of \code{quantroStat}. 
#' These are the test statistics resulting from the permuted samples.}
#' \item{quantroPvalPerm}{If \code{B} is not equal to 0, then this is the 
#' p-value associated with the proportion of times the test statistics 
#' resulting from the permuted samples were larger than \code{quantroStat}.}
#' 
#' @details 
#' Quantile normalization is one of the most widely used normalization tools 
#' for data analysis in genomics. Although it was originally developed for 
#' gene expression microarrays it is now used across many different 
#' high-throughput applications including RNAseq and ChIPseq. The methodology 
#' relies on the assumption that observed changes in the empirical 
#' distribution of samples are due to unwanted variability. Because the data is 
#' transformed to remove these differences it has the potential to remove 
#' interesting biologically driven global variation. Therefore, applying 
#' quantile normalization, or other global normalization methods
#' that rely on similar assumptions, may not be an appropriate depending 
#' on the type and source of variation. 
#' 
#' This function can be used to test a priori to the data analysis whether 
#' global normalization methods such as quantile normalization should be 
#' applied. The \code{quantro} function uses the raw unprocessed high-throughput 
#' data to test for global differences in the distributions across a set of groups. 
#' 
#' The \code{quantro} function will perform two tests:
#' 
#' 1. An ANOVA to test if the medians of the distributions are different across 
#' groups. Differences across groups could be attributed to unwanted technical 
#' variation (such as batch effects) or real global biological variation. 
#' This is a helpful step for the user to verify if there is some unaccounted
#' technical variation. 
#' 
#' 2. A test for global differences between the distributions across groups.
#' The main output is a test statistic called \code{quantroStat}. This test 
#' statistic is a ratio of two variances and is similar to the idea of ANOVA. 
#' The main idea of the test is to compare the variability of distributions 
#' within the groups to the variability of distributions between the groups. 
#' If the variance between the groups is sufficiently larger than the variance 
#' within the groups, quantile normalization may not be an appropriate 
#' normalization technique depending on the source of variation 
#' (technical or biological variation). As a default, we perform this test on 
#' after a median normalization, but this option may be changed.
#' 
#' To assess the statistical significance of \code{quantroStat}, we use 
#' permutation testing. To perform a permutation test, set \code{B} to the 
#' number of permutations which will create a null distribution.  If the number
#' of samples is large, this number can be a large number such as 1000. This 
#' step can be very slow, but a parallelization has been implemented 
#' throught the \code{foreach} package. Register the number of cores using 
#' the \code{doParallel} package. 
#' 
#' See the vignette for more details. 
#' 
#' @aliases quantro
#' 
#' @docType methods
#' @name show
#' @import Biobase minfi doParallel foreach iterators
#' @importFrom methods show
#' 
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' p <- getBeta(flowSorted, offset = 100)
#' pd <- pData(flowSorted)
#' 
#' qtest <- quantro(object = p, groupFactor = pd$CellType)
#' 
#' @rdname quantro
#' @export
quantro <- function(object, groupFactor = NULL, B = 0, qRange = NULL,
                    useMedianNormalized = TRUE, verbose = TRUE)
{ 
    if(inherits(object, "eSet")){
        if(is(object, "ExpressionSet")){ object <- exprs(object) }
        if(is(object, "MethylSet")){ object <- getBeta(object, offset = 100) }
    }

    if(is.null(groupFactor)){  
        stop("Must provide groupFactor to specify the group 
                level information associated with each sample or 
                or column in object.")
    } else {
        groupFactor <- as.factor(groupFactor)
    }

    if(ncol(object) != length(groupFactor)){
        stop("Number of columns in object does not match 
             length of groupFactor.")
    }

    nT <- ncol(object)
    groupLevels <- levels(groupFactor)
    K <- length(groupLevels)
    nk <- c(table(groupFactor))
    objectMedians <- apply(object, 2, median)
    objectMedians <- round(objectMedians, 7)

    if(length(unique(objectMedians)) == 1){
        anovaFit <- NA
        if(verbose){
            message("[quantro] All median values equal. No ANOVA performed. 
                    No median normalization.")
        }
    } else {
        anovaFit <- anova(lm(objectMedians ~ groupFactor))
        anovaPval <- (anovaFit$`Pr(>F)`[1] < 0.05)
        if(verbose){
            if(anovaPval){ 
                message("[quantro] Average medians of the distributions are 
                        not equal across groups.")
            } else { 
                message("[quantro] Average medians of the distributions are 
                        equal across groups.") 
            }
        }
    }

    if(useMedianNormalized){
        objectMedianNormalized <- sweep(object, 2, objectMedians, FUN = "-")
        objectNorm <- objectMedianNormalized
    } else {
        objectNorm <- object
    }

    if(verbose){ message("[quantro] Calculating the quantro test statistic.") }

    if(is.null(qRange)){
        Fnik = apply(objectNorm, 2, sort) 
    } else {
        Fnik = apply(objectNorm, 2, quantile, probs = qRange, na.rm = TRUE)
    }

    Fndotk = sapply( groupLevels, function(x){ 
        rowMeans(Fnik[, which(groupFactor %in% x)]) 
    } )
    
    Fndotdot = rowMeans(Fnik)
    
    betweenDiff = colMeans( (Fndotk - replicate(K, Fndotdot))^2 )
    MSb = sum(betweenDiff*nk) / (K - 1)

    withinDiff <- do.call("cbind", lapply(groupLevels, function(x){
        Fnik[, which(groupFactor %in% x) ] - 
            rowMeans(Fnik[, which(groupFactor %in% x) ]) }))

    withinDiff <- colMeans((withinDiff)^2)
    MSe <- sum(withinDiff) / (nT - K)	
    quantroStat <- MSb / MSe

    results <- new("quantro")
    results@summary <- list("nGroups" = K, 
                            "nTotSamples" = nT,
                            "nSamplesinGroups" = nk)
    results@B <- B
    results@anova <- anovaFit
    results@MSbetween <- MSb
    results@MSwithin <- MSe
    results@quantroStat <- quantroStat
    results@quantroStatPerm <- NA
    results@quantroPvalPerm <- NA

    if (B == 0){ message("[quantro] No permutation testing performed. 
                         Use B > 0 for permutation testing.") }
    if (B < 0){ stop("Must pick B greater than or equal to 0.") }
    if (B > 0){
        if(verbose){ message("[quantro] Starting permutation testing.") }
        if(!getDoParRegistered()){ registerDoSEQ() }
        
        workers <- getDoParWorkers()
        backend <- getDoParName() 
        version <- getDoParVersion()

        if (verbose) {
            if (workers == 1) {
                mes <- "[quantro] Using a single core (backend: %s, 
                        version: %s)."
                message(sprintf(mes, backend, version))
            } else {
                mes <- "[quantro] Parallelizing using %s workers/cores 
                        (backend: %s, version: %s)."
                message(sprintf(mes, workers, backend, version))
            }
        }

    one.test <- function(xstar){ 
        FndotkPerm <- sapply(groupLevels, function(x){ 
            rowMeans(Fnik[, which(xstar %in% x) ]) })
        
        betweenDiff <- colMeans( (FndotkPerm - replicate(K, Fndotdot))^2 )
        MSb <- sum(betweenDiff*nk) / (K - 1)
        
        withinDiff <- do.call("cbind", lapply(groupLevels, function(x){ 
            Fnik[, which(xstar %in% x) ] - 
                rowMeans(Fnik[, which(xstar %in% x) ])}))
        withinDiff <- colMeans( (withinDiff)^2 )
        MSe <- sum(withinDiff) / (nT - K)	

        FstatsinglePerm <- MSb / MSe
        return(data.frame(MSb, MSe, FstatsinglePerm))
    }

    xperm <- replicate(B, sample(groupFactor))
    chunksize <- ceiling(B/workers)
    subMat <- NULL
    permResults <- foreach(subMat = iter(xperm, by = "col", 
                                         chunksize = chunksize), 
                           .combine = "c") %dopar% { 
                               apply(subMat, 2, function(x){ one.test(x) }) 
                           }
    permResultDataFrame <- do.call(rbind, permResults)
    quantroPvalPerm <- mean(permResultDataFrame$FstatsinglePerm > quantroStat)

    results@quantroStatPerm <- permResultDataFrame[,3]
    results@quantroPvalPerm <- quantroPvalPerm
    }

    return(results)
}
