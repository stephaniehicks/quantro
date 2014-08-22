#' @rdname summary
#' @export
setGeneric("summary", function(object) standardGeneric("summary"))

#' @rdname anova
#' @export
setGeneric("anova", function(object) standardGeneric("anova"))

#' @rdname MSbetween
#' @export
setGeneric("MSbetween", function(object) standardGeneric("MSbetween"))

#' @rdname MSwithin
#' @export
setGeneric("MSwithin", function(object) standardGeneric("MSwithin"))

#' @rdname quantroStat
#' @export
setGeneric("quantroStat", function(object) standardGeneric("quantroStat"))

#' @rdname quantroStatPerm
#' @export
setGeneric("quantroStatPerm", function(object) 
    standardGeneric("quantroStatPerm"))

#' @rdname quantroPvalPerm
#' @export
setGeneric("quantroPvalPerm", function(object) 
    standardGeneric("quantroPvalPerm"))
