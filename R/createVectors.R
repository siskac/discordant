#' @import Biobase
#' @import biwt
#' @import gtools
#' @import MASS
#' @import stats
#' @import tools
#' 
#' @export
createVectors <- function(x, y = NULL, groups, 
                          cor.method = c("spearman", "pearson", "bwmc", 
                                         "sparcc")) {
    cor.method <- match.arg(cor.method)
    .checkCreateVectorsInputs(x, y, groups, cor.method)
    #print(x)
    
    index1 <- which(groups == 1)
    index2 <- which(groups == 2)
    x <- exprs(x)
    
    if(is.null(y)) {
        data <- x
    } else {
        y <- exprs(y)
        data <- rbind(x, y)
        featureSize <- dim(x)[1]
    }
    
    data1 <- data[,index1]
    data2 <- data[,index2]
    
    if (cor.method == "spearman" || cor.method == "pearson") {
        statMatrix1 <- cor(t(data1), method = cor.method)
        statMatrix2 <- cor(t(data2), method = cor.method)
    } else if (cor.method == "bwmc") {
        statMatrix1 <- biwt.cor(data1)
        statMatrix2 <- biwt.cor(data2)
    } else if (cor.method == "sparcc") {
        statMatrix1 <- SparCC.count(t(data1))$cor.w
        statMatrix2 <- SparCC.count(t(data2))$cor.w
    }
    
    if (is.null(y)) {
        statVector1 <- as.vector(statMatrix1)
        statVector2 <- as.vector(statMatrix2)
        diag <- lower.tri(statMatrix1, diag = FALSE)
        statVector1 <- statMatrix1[diag]
        statVector2 <- statMatrix2[diag]
    } else {
        statMatrix1 <- statMatrix1[1:featureSize,
                                   (featureSize + 1):dim(data1)[1]]
        statMatrix2 <- statMatrix2[1:featureSize,
                                   (featureSize + 1):dim(data1)[1]]
        statVector1 <- as.vector(statMatrix1)
        statVector2 <- as.vector(statMatrix2)
    }
    
    vector_names <- getNames(x, y)
    names(statVector1) <- vector_names
    names(statVector2) <- vector_names
    return(list(v1 = statVector1, v2 = statVector2))
}

.checkCreateVectorsInputs <- function(x, y = NULL, groups, cor.method) {
    
    if (unique(unique(groups) != c(1,2))) {
        stop("groups vector must consist of 1s and 2s corresponding to first
             and second group")
    }
    
    if (!is(x, "ExpressionSet") || (!is.null(y) && !is(y, "ExpressionSet"))) {
        stop("x and y (if present) must be type ExpressionSet")
    }
    
    x <- exprs(x)
    if(is.null(y)) {
        data <- x
    } else {
        y <- exprs(y)
        data <- rbind(x, y)
    }
    if (cor.method == "sparcc" && min(data) < 0) {
        stop("Negative values found. SparCC can only be applied to ",
             "sequencing data.")
    }
}