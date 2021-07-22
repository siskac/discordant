#' @import Biobase
#' @import biwt
#' @import gtools
#' @import MASS
#' @import stats
#' @import tools
#' 
#' @export
createVectors <- function(x, y = NULL, groups, cor.method = c("spearman")) {
    print(x)
    if(checkInputs(x,y,groups)) {
        stop("Please fix inputs.")
    }
    
    index1 <- which(groups == 1)
    index2 <- which(groups == 2)
    
    check <- match(cor.method, c("spearman", "pearson", "bwmc", "sparcc"))
    if(is.na(check)) {
        stop("Please enter spearman, pearson, bwmc or sparcc for correlation 
            metric.")
    }
    
    x <- exprs(x)
    
    if(is.null(y) == FALSE) {
        y <- exprs(y)
        data <- rbind(x, y)
        data1 <- data[,index1]
        data2 <- data[,index2]
        featureSize = dim(x)[1]
        if(cor.method == c("spearman") || cor.method == c("pearson")) {
            statMatrix1 <- cor(t(data1), method = cor.method)
            statMatrix2 <- cor(t(data2), method = cor.method)
        }
        if(cor.method == c("bwmc")) {
            statMatrix1 <- biwt.cor(data1)
            statMatrix2 <- biwt.cor(data2)
        }
        if(cor.method == c("sparcc")) {
            if(min(data) < 0) {
                stop("SparCC can only be applied to sequencing data.")
            }
            statMatrix1 <- SparCC.count(t(data1))
            statMatrix2 <- SparCC.count(t(data2))
        }
        statMatrix1 <- statMatrix1[1:featureSize,
                                   (featureSize + 1):dim(data1)[1]]
        statMatrix2 <- statMatrix2[1:featureSize,
                                   (featureSize + 1):dim(data1)[1]]
        statVector1 <- as.vector(statMatrix1)
        statVector2 <- as.vector(statMatrix2)
        vector_names <- getNames(x,y)
    }
    
    if(is.null(y)) {
        data <- x # NOTE: Line added by MM to fix min(data) error below
        data1 <- x[,index1]
        data2 <- x[,index2]
        if(cor.method == c("spearman") || cor.method == c("pearson")) {
            statMatrix1 <- cor(t(data1), method = cor.method)
            statMatrix2 <- cor(t(data2), method = cor.method)
        }   
        if(cor.method == c("bwmc")) {
            statMatrix1 <- biwt.cor(data1)
            statMatrix2 <- biwt.cor(data2)
        }   
        if(cor.method == c("sparcc")) {
            if(min(data) < 0) {
                stop("SparCC can only be applied to sequencing data.")
            }   
            # NOTE: $cor.w added to each statement below by MM as SparCC.count
            #    returns two vectors, one for correlation and one for covariance
            statMatrix1 <- SparCC.count(t(data1))$cor.w
            statMatrix2 <- SparCC.count(t(data2))$cor.w
        }   
        statVector1 <- as.vector(statMatrix1)
        statVector2 <- as.vector(statMatrix2)
        diag <- lower.tri(statMatrix1, diag = FALSE)
        statVector1 <- statMatrix1[diag]
        statVector2 <- statMatrix2[diag]
        vector_names <- getNames(x)
    }
    
    names(statVector1) <- vector_names
    names(statVector2) <- vector_names
    
    return(list(v1 = statVector1, v2 = statVector2))
}