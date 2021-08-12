#' Fisher Transformation of Pearson Correlation Coefficients to Z Scores
#' 
#' Transforms Pearsons correlation coefficients into z scores using Fishers 
#' method.
#' 
#' @param rho Integer or numeric vector of Pearson's correlation coefficients
#' @return Returns Fisher-transformed correlation coefficients
#' @references 
#' Fisher, R.A. (1915). "Frequency distribution of the values of the correlation
#' coefficient in samples of an indefinitely large population". Biometrika
#' (Biometrika Trust) 10 (4).
#' @details Fisher's transformation is when correlation coefficients are 
#' transformed into a z score. These z scores have an approximately normal 
#' distribution.
#' @examples 
#' ## Create integer or list of Pearson's correlation coefficients.
#' 
#' library(MASS)
#' rhoV <- as.vector(cor(t(mvrnorm(10,rep(3,100),diag(100)))))
#' 
#' ## Determine Fisher-Transformed z scores of rho
#' zV <- fishersTrans(rhoV)
#' 
#' @export
fishersTrans <- function(rho) {
    r = (1 + rho) / (1 - rho)
    z = 0.5 * log(r, base = exp(1))
    return(z)
}

.subSampleData <- function(pdata, class, mu, sigma, nu, tau, pi, components) {
    n <- as.integer(dim(pdata)[1])
    g <- as.integer(nlevels(as.factor(class)))
    
    yl.outer <- function(k, zx, zy){
        return( c(zx[k,] %o% zy[k,]) )
    }
    
    zx <- unmap(class[,1], components = components)
    zy <- unmap(class[,2], components = components)
    zxy <- sapply(1:dim(zx)[1], yl.outer, zx, zy)
    
    results <- subsampling_cpp(as.double(pdata[,1]), as.double(pdata[,2]), 
                               as.double(t(zxy)), n, as.double(pi), 
                               as.double(mu), as.double(sigma), as.double(nu), 
                               as.double(tau), g)
    
    return(list(pi = t(array(results[[5]], dim=c(g,g))), 
                mu_sigma = rbind(results[[6]], results[[7]]), 
                nu_tau = rbind(results[[8]], results[[9]]), 
                class = apply(array(results[[3]], dim = c(n, g*g)), 1, order,
                              decreasing=TRUE)[1,], 
                z = array(results[[3]], dim=c(n,g*g))))
} 

# modified from package mclust
.unmap <- function(classification, components){
    n <- length(classification)
    # u <- sort(unique(classification)) # OG Code
    u <- 0:(components - 1) # Max's potential fix
    labs <- as.character(u)
    k <- length(u)
    z <- matrix(0, n, k)
    for (j in 1:k) z[classification == u[j], j] <- 1
    dimnames(z) <- list(NULL, labs)
    return(z)
}

.getNames <- function(x, y = NULL) {
    if(is.null(y) == FALSE) {
        y <- exprs(y)
        namesMatrix <- NULL
        for(i in 1:nrow(x)) {
            tempMatrix <- cbind(rep(rownames(x)[i], nrow(y)), rownames(y))
            namesMatrix <- rbind(namesMatrix, tempMatrix)
        }
    } else {
        temp <- matrix(NA,nrow = nrow(x), ncol = nrow(x))
        diag <- lower.tri(temp, diag = FALSE)
        temp[diag] <- rep(1, sum(diag == TRUE))
        
        namesMatrix <- NULL
        
        for(i in 1:dim(temp)[1]) {
            outputCol <- temp[,i]
            index <- which(is.na(outputCol) == FALSE)
            if(length(index) > 0) {
                tempMatrix <- cbind(rep(rownames(x)[i],length(index)), 
                                    rownames(x)[index])
                namesMatrix <- rbind(namesMatrix, tempMatrix)
            }
        }
    }
    
    vector_names <- apply(namesMatrix, 1, function(k) paste(k[1],"_",k[2],
                                                            sep = ""))
    return(vector_names)
}