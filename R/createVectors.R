#' Create correlation coefficient vectors based on bivariate data
#' 
#' Calculates correlation coefficients based on two groups of -omics bivariate 
#' data. Currently, only two groups of samples can be specified. Used to make 
#' input for discordantRun().
#' 
#' @param x ExpressionSet of -omics data
#' @param y Optional second ExpressionSet of -omics data, induces dual -omics 
#' analysis
#' @param groups n-length vector of 1s and 2s matching samples belonging to 
#' groups 1 and 2
#' @param cor.method Correlation method to measure association. Options are 
#' "spearman", "pearson", "bwmc" and "sparcc"
#' @details Creates vectors of correlation coefficents based on feature pairs 
#' within x or between x and y. The names of the vectors are the feature pairs 
#' taken from x and y.
#' @return List of two named numeric vectors. Vectors give the correlation
#' coefficients for groups 1 and 2 respectively, and vector names give
#' the each feature for the resptive feature pair seperated by an underscore.
#' @author Charlotte Siska \email{siska.charlotte@@gmail.com}
#' @author Max McGrath \email{max.mcgrath@@ucdenver.edu}
#' @examples 
#' 
#' ## load data
#' data("TCGA_GBM_miRNA_microarray")
#' data("TCGA_GBM_transcript_microarray")
#' print(colnames(TCGA_GBM_transcript_microarray)) # look at groups
#' groups <- c(rep(1,10), rep(2,20))
#' # transcript-transcript pairs
#' vectors <- createVectors(TCGA_GBM_transcript_microarray, 
#'                          groups = groups, cor.method = c("pearson"))
#' # miRNA-transcript pairs
#' vectors <- createVectors(TCGA_GBM_transcript_microarray, 
#'                          TCGA_GBM_miRNA_microarray, groups = groups)
#'                          
#' @references 
#' Siska C, Bowler R and Kechris K. The Discordant Method: A Novel Approach for 
#' Differential Correlation. (2015) Bioinformatics. 32(5): 690-696.
#' 
#' Friedman J and Alm EJ. Inferring Correlation Networks from Genomic Survey 
#' Data. (2012) PLoS Computational Biology. 8:9, e1002687.
#' 
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
    
    if(is.null(y)) {
        data <- exprs(x)
        vector_names <- .getNames(exprs(x), y)
    } else {
        data <- rbind(exprs(x), exprs(y))
        vector_names <- .getNames(exprs(x), exprs(y))
        featureSize <- dim(exprs(x))[1]
    }
    
    data1 <- data[, which(groups == 1)]
    data2 <- data[, which(groups == 2)]
    
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