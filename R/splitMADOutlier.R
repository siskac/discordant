#' Outliers using left and right MAD
#' 
#' Identify features with outliers using left and right median absolute 
#' deviation (MAD).
#' 
#' @param mat m by n matrix of -omics data, where rows are features and columns 
#' samples.
#' @param filter0 Option to filter out features if they have at least one 0 
#' value. Default is TRUE.
#' @param threshold Threshold of how many MADs outside the left or right median 
#' is used to determine features with outliers.
#' @return
#' \describe{
#' \item{mat.filtered}{Input matrix where features with outliers filtered out.}
#' \item{index}{Index of features that have no outliers.}
#' }
#' @details The purpose of this function is to determine outliers in 
#' non-symmetric distributions. The distribution is split by the median. 
#' Outliers are identifed by being however many median absolute deviations (MAD)
#' from either split distribution.
#' 
#' @references 
#' Leys C, Klein O, Bernard P and Licata L. "Detecting Outliers: Do Not Use 
#' Standard Deviation Around the Mean, Use Absolute Deivation Around the 
#' Median." Journal of Experimental Social Psychology, 2013. 49(4), 764-766.
#' 
#' Magwene, PM, Willis JH, Kelly JK and Siepel A. "The Statistics of Bulk 
#' Segregant Analysis Using Next Generation Sequencing." PLoS Computational 
#' Biology, 2011. 7(11), e1002255.
#' 
#' 
#' @examples 
#' ## Simulate matrix of continuous -omics data.
#' data(TCGA_Breast_miRNASeq)
#' 
#' ## Filter matrix based on outliers.
#' mat.filtered <- splitMADOutlier(TCGA_Breast_miRNASeq)$mat.filtered
#' 
#' @export
splitMADOutlier <- function(mat, filter0 = TRUE, threshold = 2) {
    if(mode(mat) != "S4") {
        stop("data matrix mat must be type ExpressionSet")
    }
    mat <- exprs(mat)
    maxMAD <- c()
    mat.filtered <- NULL
    index <- c()
    for(i in 1:nrow(mat)) {
        y <- mat[i,]
        x <- y
        if(length(which(is.infinite(x))) >0 ){
            x <- x[-which(is.infinite(x))]
        }
        if(length(which(is.na(x)))> 0) {
            x <- x[-which(is.na(x))]
        }
        median.x <- median(x)
        left <- x[x <= median.x]
        right <- x[x > median.x]
        left.mad <- mad(left)
        right.mad <- mad(right)
        leftThresh <- median(left) - left.mad*threshold
        rightThresh <- median(right) + right.mad*threshold 
        left.check <- sum(left < leftThresh)
        right.check <- sum(right > rightThresh)
        if(left.check == 0 & right.check == 0) {
            mat.filtered <- rbind(mat.filtered, y)
            rownames(mat.filtered)[nrow(mat.filtered)] <- rownames(mat)[i]
            index <- c(index,i)
        }
    }
    colnames(mat.filtered) = colnames(mat)
    return(list(mat.filtered = mat.filtered, index = index))
}