#' @title Density plots of columns in a matrix
#' 
#' @description Plots the density of the columns of a matrix
#' 
#' @param object object an object which is inherited from an 
#' \code{eSet} such as an \code{ExpressionSet} or 
#' \code{MethylSet} object. The \code{object} can also be a 
#' data frame or matrix with observations
#' (e.g. probes or genes) on the rows and samples as the 
#' columns.  
#' @param groupFactor an optional factor variable representing which 
#' group each column in \code{object} belongs to. It is important 
#' that values in \code{groupFactor} be in the same 
#' order as the columns in \code{object}. 
#' @param type the type of lines to plot. Default type is line ("l"). 
#' @param col a numeric of colors to plot each group.
#' @param lty the line type. Default is the solid line. 
#' @param ... other arguments that can be passed to the 
#' code{matplot} function. 
#' 
#' @return A density plot for each column in \code{object}
#' 
#' @export
#' @examples
#' library(minfi)
#' data(flowSorted)
#' p <- getBeta(flowSorted, offset = 100)
#' pd <- pData(flowSorted)
#' matdensity(object = p, groupFactor = pd$CellType, col = c(2,3), 
#' xlab = "beta values", ylab = "density")
matdensity <- function(object, groupFactor = NULL, col = NULL, type = "l", 
                       lty = 1, ...){
    
    if(inherits(object, "eSet")){
        if(is(object, "ExpressionSet")){ object <- exprs(object) }
        if(is(object, "MethylSet")){ object <- getBeta(object, offset = 100) }
    }
    
    getdensity  <- function(object){
        min.object <- min(object[is.finite(object)], na.rm = TRUE)
        max.object <- max(object[is.finite(object)], na.rm = TRUE)
        densityMat <- apply(object, 2, function(z){
            density(z, from = min.object, to = max.object, na.rm = TRUE)$y
        })
        x = seq(from=min.object, to=max.object, length.out = nrow(densityMat))
        list(densityMat = densityMat, x = x)
    }
    output <- getdensity(object)

    if(!is.null(groupFactor)){
        groupFactor <- as.factor(groupFactor)
        if(!is.null(col)){
            col = col[as.integer(groupFactor)]
        } else {
            col = seq_len(length(levels(groupFactor)))[as.integer(groupFactor)]
        }
    } else { 
        if(!is.null(col)){
            col = rep(col, ncol(object))
        } else {
            col = rep(1, ncol(object))
        }
    }

    matplot(x = output$x, output$densityMat, type = type, col = col, 
            lty = lty, ...)
}

#' @title Box plots of columns in a matrix
#' 
#' @description Box plots of the columns of a matrix, but 
#' the columns are ordered and colored by a group-level 
#' variable 
#' 
#' @param object object an object which is inherited from an 
#' \code{eSet} such as an \code{ExpressionSet} or 
#' \code{MethylSet} object. The \code{object} can also be a 
#' data frame or matrix with observations
#' (e.g. probes or genes) on the rows and samples as the 
#' columns.  
#' @param groupFactor a factor variable representing which 
#' group each column in \code{object} belongs to. It is important 
#' that values in \code{groupFactor} be in the same 
#' order as the columns in \code{object}.
#' @param col a numeric of colors to plot each group.
#' @param las a numeric in (0, 1, 2, 3) to orient the axis labels. 
#' Default is 3 (always vertical). 
#' @param ... other arguments that can be passed to the 
#' code{boxplot} function. 
#' 
#' @return A box plot for each column in \code{object}
#' 
#' @export
#' @examples
#' library(minfi)
#' data(flowSorted)
#' 
#' p <- getBeta(flowSorted, offset = 100)
#' pd <- pData(flowSorted)
#' matboxplot(object = p, groupFactor = pd$CellType, col = c(2,3))
matboxplot <- function (object, groupFactor, col = NULL, las = 3, ...){
    
    if(inherits(object, "eSet")){
        if(is(object, "ExpressionSet")){ object <- exprs(object) }
        if(is(object, "MethylSet")){ object <- getBeta(object, offset = 100) }
    }

    groupFactor <- as.factor(as.character(groupFactor))
    objectordered <- object[, order(groupFactor)]
    groupFactorOrdered <- groupFactor[order(groupFactor)]

    if(is.null(col)){
        col = seq_len(length(levels(groupFactor)))[as.integer(
            groupFactorOrdered)]
    } else {
        col = col[as.integer(groupFactorOrdered)]
    }
    
    boxplot(objectordered, col = col, las = las, ...)
}


#' @title Plot results from \code{quantro} function.
#' 
#' @description This function plots the histogram of the
#' null test statistics from permutation test in the \code{quantro}
#' function. 
#' 
#' @param object a quantro object from \code{quantro} 
#' @param savePlot a TRUE/FALSE object argument
#' determining if the plot will be saved for further use or
#' immediately displayed on the screen.
#' @param xLab label for x-axis
#' @param yLab label for y-axis
#' @param mainLab title of plot
#' @param binWidth binwidth or histogram. Default is the stat_bin default. 
#' 
#' @return A histogram will be plotted containing the null 
#' test statistics when using boostrapped samples. 
#' The red line is the observed test statistic
#' \code{quantroStat} from \code{quantro()}.
#' 
#' @import ggplot2 
#' 
#' @export
#' @examples
#' \donttest{
#' library(minfi)
#' data(flowSorted)
#' p <- getBeta(flowSorted, offset = 100)
#' pd <- pData(flowSorted)
#' 
#' library(doParallel)
#' registerDoParallel(cores=4)
#' qtest <- quantro(p, pd$CellType, B = 100)
#' quantroPlot(qtest)
#' }
quantroPlot <- function(object, savePlot = FALSE, xLab = NULL, 
                        yLab = NULL, mainLab = NULL, binWidth = NULL) 
{
    B <- object@B
    quantroStat <- quantroStat(object)
    quantroStatPerm <- quantroStatPerm(object)
    quantroPvalPerm <- quantroPvalPerm(object)

    if(B == 0){
        stop("quantro was run without any permutation testing. Must define 
             B > 0 in quantro() to plot distribution of test 
             statistics from permutation testing. ")
    }  

    if(is.null(xLab)){ 
        xLab <- paste0("Histogram of test statistics \n from permutation 
                       test with B = ", B)
    }

    if(is.null(yLab)){ yLab <- "density" }

    if(is.null(mainLab)){ mainLab <- paste0(
        "quantroStat = ", round(quantroStat, 3), 
        ", quantroPvalPerm = ", round(quantroPvalPerm, 3))
    }

    value = NULL
    dat = data.frame("value" = quantroStatPerm)
    gmat <- ggplot(dat, aes(x = value)) + 
                geom_histogram(binwidth = binWidth) + 
                geom_vline(xintercept = quantroStat, colour = "red") +
                xlab(xLab) + ylab(yLab) + labs(title = mainLab)

    if (savePlot) { 
        suppressMessages(gmat)
    } else { 
        suppressMessages(print(gmat)) 
    }
}