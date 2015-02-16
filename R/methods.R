#' @title Accessors for the 'summary' slot of a quantro object.
#' 
#' @description Accessors for the 'summary' slot of a quantro object.
#' 
#' @usage
#' \S4method{summary}{quantro}(object)
#'
#' @docType methods
#' @name summary
#' @rdname summary
#' @aliases summary summary,quantro-method
#' @param object a \code{quantro} object
#' @param ... other 
#' 
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' summary(qtest)
summary.quantro <- function(object) object@summary

#' @rdname summary
#' @export
setMethod("summary", signature(object="quantro"), summary.quantro)



#' Accessors for the 'anova' slot of a quantro object.
#' 
#' Accessors for the 'anova' slot of a quantro object.
#' 
#' @usage
#' \S4method{anova}{quantro}(object)
#'
#' @docType methods
#' @name anova
#' @rdname anova
#' @aliases anova anova,quantro-method
#' @param object a \code{quantro} object
#' @param ... other
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' anova(qtest)
anova.quantro <- function(object) object@anova

#' @rdname anova
#' @export
setMethod("anova", signature(object="quantro"), anova.quantro)


#' Accessors for the 'MSbetween' slot of a quantro object.
#' 
#' Accessors for the 'MSbetween' slot of a quantro object.
#' 
#' @usage
#' \S4method{MSbetween}{quantro}(object)
#'
#' @docType methods
#' @name MSbetween
#' @rdname MSbetween
#' @aliases MSbetween MSbetween,quantro-method
#' @param object a \code{quantro} object
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' MSbetween(qtest)
MSbetween.quantro <- function(object) object@MSbetween

#' @rdname MSbetween
#' @export
setMethod("MSbetween", signature(object="quantro"), MSbetween.quantro)


#' Accessors for the 'MSwithin' slot of a quantro object.
#' 
#' Accessors for the 'MSwithin' slot of a quantro object.
#' 
#' @usage
#' \S4method{MSwithin}{quantro}(object)
#'
#' @docType methods
#' @name MSwithin
#' @rdname MSwithin
#' @aliases MSwithin MSwithin,quantro-method
#' @param object a \code{quantro} object
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' MSwithin(qtest)
MSwithin.quantro <- function(object) object@MSwithin

#' @rdname MSwithin
#' @export
setMethod("MSwithin", signature(object="quantro"), MSwithin.quantro)




#' Accessors for the 'quantroStat' slot of a quantro object.
#' 
#' Accessors for the 'quantroStat' slot of a quantro object.
#' 
#' @usage
#' \S4method{quantroStat}{quantro}(object)
#'
#' @docType methods
#' @name quantroStat
#' @rdname quantroStat
#' @aliases quantroStat quantroStat,quantro-method
#' @param object a \code{quantro} object
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' quantroStat(qtest)
quantroStat.quantro <- function(object) object@quantroStat

#' @rdname quantroStat
#' @export
setMethod("quantroStat", signature(object="quantro"), quantroStat.quantro)


#' Accessors for the 'quantroStatPerm' slot of a quantro object.
#' 
#' Accessors for the 'quantroStatPerm' slot of a quantro object.
#' 
#' @usage
#' \S4method{quantroStatPerm}{quantro}(object)
#'
#' @docType methods
#' @name quantroStatPerm
#' @rdname quantroStatPerm
#' @aliases quantroStatPerm quantroStatPerm,quantro-method
#' @param object a \code{quantro} object
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' quantroStatPerm(qtest)
quantroStatPerm.quantro <- function(object) object@quantroStatPerm

#' @rdname quantroStatPerm
#' @export
setMethod("quantroStatPerm", signature(object="quantro"), 
          quantroStatPerm.quantro)


#' Accessors for the 'quantroPvalPerm' slot of a quantro object.
#' 
#' Accessors for the 'quantroPvalPerm' slot of a quantro object.
#' 
#' @usage
#' \S4method{quantroPvalPerm}{quantro}(object)
#'
#' @docType methods
#' @name quantroPvalPerm
#' @rdname quantroPvalPerm
#' @aliases quantroPvalPerm quantroPvalPerm,quantro-method
#' @param object a \code{quantro} object
#'
#' @examples
#' library(minfi)
#' data(flowSorted)
#' pd <- pData(flowSorted)
#' qtest <- quantro(flowSorted, groupFactor = pd$CellType)
#' quantroPvalPerm(qtest)
quantroPvalPerm.quantro <- function(object) object@quantroPvalPerm

#' @rdname quantroPvalPerm
#' @export
setMethod("quantroPvalPerm", signature(object="quantro"), 
          quantroPvalPerm.quantro)


