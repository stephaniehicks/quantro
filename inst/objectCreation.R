#' @name flowSorted
#' @title A subset of FlowSorted.DLPFC.450k data set
#'
#' @description 
#' This is the script used to create a subset of the 
#' FlowSorted.DLPFC.450k data set from Bioconductor. 
#' The purpose is to create an example object for the man
#' pages and vignette in this package. 
#' 
#' @docType data
#' @keywords datasets
#' @format A MethylSet object with 1e4 rows (probes) 
#' and 58 columns (samples).
#' 
#' @rdname flowSorted
#' @export
#' 
library(FlowSorted.DLPFC.450k)
library(minfi)
data(FlowSorted.DLPFC.450k)

Mset <- preprocessRaw(FlowSorted.DLPFC.450k)
pd <- pData(Mset)
meth <- getMeth(Mset)
unmeth <- getUnmeth(Mset)

# Subset probes to speed up computation time for vignette and man pages
set.seed(123)
subRows <- sample(seq_len(nrow(meth)), size = 1e4)
meth <- meth[subRows,]
unmeth <- unmeth[subRows,]

# Create MethylSet object
flowSorted <- MethylSet(Meth = meth, Unmeth = unmeth, 
                        phenoData = AnnotatedDataFrame(pd))

# Save MethylSet object to use in vignette and man pages
save(flowSorted, file = "data/flowSorted.RData")

