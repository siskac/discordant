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