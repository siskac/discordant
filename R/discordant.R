# Applied from Lai et al 2007 Bioinformatics.

# Code written by Charlotte Siska
# R functions fisherTrans, createVectors, discordantRun and makeTable
# Email: charlotte.siska@ucdenver.edu

# Code written by Yinglei Lai
# C code and R functions unmap and em.normal.partial.concordant
# E-mail: ylai@gwu.edu

# Code written by Fang Huaying
# R function for SparCC, adapted from https://bitbucket.org/yonatanf/sparcc
# Email: hyfang@pku.edu.cn

################################################################################
# File: SparCC.R
# Aim : SparCC 
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 11/12/2014
#-------------------------------------------------------------------------------
# SparCC for counts known
#   function: SparCC.count
#   input:
#          x ------ nxp count data matrix, row is sample, col is variable
#       imax ------ resampling times from posterior distribution. default 20
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
#

SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # dimension for w (latent variables)
  p <- ncol(x);
  n <- nrow(x);
  # posterior distribution (alpha)
  x <- x + 1;
  # store generate data
  y <- matrix(0, n, p);
  # store covariance/correlation matrix
  cov.w <- cor.w <- matrix(0, p, p);
  indLow <- lower.tri(cov.w, diag = TRUE);
  # store covariance/correlation for several posterior samples
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
  for(i in 1:imax) {
    # generate fractions from posterior distribution
    y <- t(apply(x, 1, function(x) rdirichlet(n = 1, alpha = x)));
    # estimate covariance/correlation
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
    # store variance/correlation only low triangle 
    covs[, i] <- cov_cor$cov.w[indLow];
    cors[, i] <- cov_cor$cor.w[indLow];
  }
  # calculate median for several posterior samples
  cov.w[indLow] <- apply(covs, 1, median); 
  cor.w[indLow] <- apply(cors, 1, median);
  #
  cov.w <- cov.w + t(cov.w);
  diag(cov.w) <- diag(cov.w) / 2;
  cor.w <- cor.w + t(cor.w);
  diag(cor.w) <- 1;
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}

#-------------------------------------------------------------------------------
# SparCC for fractions known
#   function: SparCC.frac
#   input:
#          x ------ nxp fraction data matrix, row is sample, col is variable
#       kmax ------ max iteration steps for SparCC. default is 10
#      alpha ------ the threshold for strong correlation. default is 0.1
#       Vmin ------ minimal variance if negative variance appears. default is 1e-4
#   output: a list structure
#      cov.w ------ covariance estimation
#      cor.w ------ correlation estimation
SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # Log transformation
  x <- log(x);
  p <- ncol(x);
  # T0 = var(log(xi/xj)) variation matrix
  TT <- stats::var(x);
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
  # Variance and correlation coefficients for Basic SparCC  
  rowT0 <- rowSums(T0);
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
  var.w[var.w < Vmin] <- Vmin;
  #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
  #  sqrt(outer(var.w, var.w, "*")) / 2;
  Is <- sqrt(1/var.w);
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
  # Truncated correlation in [-1, 1]
  cor.w[cor.w <= - 1] <- - 1; 
  cor.w[cor.w >= 1] <- 1;
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1; 
  # Remove pairs
  rp <- NULL;
  # Left components
  cp <- rep(TRUE, p);
  # Do loops until max iteration or only 3 components left
  k <- 0;  
  while(k < kmax && sum(cp) > 3) {
    # Left T0 = var(log(xi/xj)) after removing pairs
    T02 <- T0;
    # Store current correlation to find the strongest pair
    curr_cor.w <- cor.w;
    # Remove diagonal
    diag(curr_cor.w) <- 0;
    # Remove removed pairs
    if(!is.null(rp)) {
      curr_cor.w[rp] <- 0;
    }
    # Find the strongest pair in vector form
    n_rp <- which.max(abs(curr_cor.w));
    # Remove the pair if geater than alpha
    if(abs(curr_cor.w[n_rp]) >= alpha) {
      # Which pair in matrix form
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
      # Update remove pairs
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
      rp <- c(rp, n_rp);
      # Update T02
      T02[rp] <- 0;
      # Which component left
      cp <- (diag(Lmat) > 0);
      # Update variance and truncated lower by Vmin
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
      var.w[var.w <= Vmin] <- Vmin;
      # Update correlation matrix and truncated by [-1, 1]
      #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
      #  sqrt(outer(var.w, var.w, "*")) / 2;    
      Is <- sqrt(1/var.w);
      cor.w <- (var.w + rep(var.w, each = p) - T0) * 
        Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor.w[cor.w <= - 1] <- - 1;
      cor.w[cor.w >= 1] <- 1;
    }
    else {
      break;
    }
    # 
    k <- k + 1;
  }
  # Covariance
  Is <- sqrt(var.w);
  cov.w <- cor.w * Is * rep(Is, each = p);
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}



#modified from package mclust
unmap <- function(classification){
    n <- length(classification)
    u <- sort(unique(classification))
    labs <- as.character(u)
    k <- length(u)
    z <- matrix(0, n, k)
    for (j in 1:k) z[classification == u[j], j] <- 1
    dimnames(z) <- list(NULL, labs)
    return(z)
}


em.normal.partial.concordant <- function(data, class, tol=0.001, restriction=0, constrain=0, iteration=1000){
    n <- as.integer(dim(data)[1])
    g <- as.integer(nlevels(as.factor(class)))

    yl.outer <- function(k, zx, zy){
        return( c(zx[k,] %o% zy[k,]) )
    }
    yl.diag <- function(k, z){
        return( c(diag(z[k,])) )
    }

    zx <- unmap(class[,1])
    zy <- unmap(class[,2])
    zxy <- sapply(1:dim(zx)[1], yl.outer, zx, zy)

    pi <- double(g*g)
    mu <- double(g)
    sigma <- double(g)
    nu <- double(g)
    tau <- double(g)
    loglik <- double(1)
    convergence <- integer(1)
    results <- .C("em_normal_partial_concordant", as.double(data[,1]), as.double(data[,2]), as.double(t(zxy)), n, pi, mu, sigma, nu, tau, g, loglik, as.double(tol), as.integer(restriction), as.integer(constrain), as.integer(iteration), convergence)
    return(list(model="PCD", convergence=results[[16]], pi=t(array(results[[5]],dim=c(g,g))), mu_sigma=rbind(results[[6]], results[[7]]), nu_tau=rbind(results[[8]], results[[9]]), loglik=results[[11]], class=apply(array(results[[3]], dim=c(n,g*g)),1,order,decreasing=TRUE)[1,], z=array(results[[3]], dim=c(n,g*g))))
}

subSampleData <- function(pdata, class, mu, sigma, nu, tau, pi) {
    n <- as.integer(dim(pdata)[1])
    g <- as.integer(nlevels(as.factor(class)))

    yl.outer <- function(k, zx, zy){
        return( c(zx[k,] %o% zy[k,]) )
    }
    yl.diag <- function(k, z){
        return( c(diag(z[k,])) )
    }

    zx <- unmap(class[,1])
    zy <- unmap(class[,2])
    zxy <- sapply(1:dim(zx)[1], yl.outer, zx, zy)

    results <- .C("subsampling", as.double(pdata[,1]), as.double(pdata[,2]), as.double(t(zxy)), n, as.double(pi), as.double(mu), as.double(sigma), as.double(nu), as.double(tau), g)
    return(list(pi=t(array(results[[5]],dim=c(g,g))), mu_sigma=rbind(results[[6]], results[[7]]), nu_tau=rbind(results[[8]], results[[9]]), class=apply(array(results[[3]], dim=c(n,g*g)),1,order,decreasing=TRUE)[1,], z=array(results[[3]], dim=c(n,g*g))))

} 

fishersTrans <- function(rho) {
    r = (1+rho)/(1-rho)
    z = 0.5*log(r,base = exp(1))
    return(z)
}

getNames <- function(x, y = NULL) {
    if(is.null(y) == FALSE) {
        namesMatrix <- NULL
        for(i in 1:nrow(x)) {
            tempMatrix <- cbind(rep(rownames(x)[i],nrow(y)), rownames(y))
            namesMatrix <- rbind(namesMatrix, tempMatrix)
        }
    }  

    if(is.null(y)) {
        temp <- matrix(NA,nrow = nrow(x), ncol = nrow(x))
        diag <- lower.tri(temp, diag = FALSE)
        temp[diag] <- rep(1, sum(diag == TRUE))

        namesMatrix <- NULL

        for(i in 1:dim(temp)[1]) {
            outputRow <- temp[i,]
            index <- which(is.na(outputRow) == FALSE)
            if(length(index) > 0) {
                tempMatrix <- cbind(rep(rownames(x)[i],length(index)), rownames(x)[index])
                namesMatrix <- rbind(namesMatrix, tempMatrix)
            }
        }
    }

    vector_names <- apply(namesMatrix, 1, function(k) paste(k[1],"_",k[2],sep = ""))
    return(vector_names)
}

createVectors <- function(x, y = NULL, groups, cor.method = c("spearman")) {
    index1 <- which(groups == 1)
    index2 <- which(groups == 2)

    check <- match(cor.method, c("spearman", "pearson", "bwmc", "sparcc"))
    if(is.na(check)) {
        stop("Please enter spearman, pearson, bwmc or sparcc for correlation metric.")
    }

    if(is.null(y) == FALSE) {
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
        statMatrix1 <- statMatrix1[1:featureSize,(featureSize + 1):dim(data1)[1]]
        statMatrix2 <- statMatrix2[1:featureSize,(featureSize + 1):dim(data1)[1]]
        statVector1 <- as.vector(statMatrix1)
        statVector2 <- as.vector(statMatrix2)
        vector_names <- getNames(x,y)
    }

    if(is.null(y)) {
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
        statMatrix1 <- SparCC.count(t(data1))
        statMatrix2 <- SparCC.count(t(data2))
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

discordantRun <- function(v1, v2, x, y = NULL, transform = TRUE, subsampling = FALSE, subSize = dim(x)[1], iter = 100, components = 3) {

    if(transform == TRUE) {
        if(range(v1)[1] < -1 || range(v1)[2] > 1 || range(v2)[1] < -1 || range(v2)[2] > 1) {
            stop("correlation vectors have values less than -1 and/or greater than 1.")
        }
        v1 <- fishersTrans(v1)
        v2 <- fishersTrans(v2)
    }

    check <- match(components, c(3,5))
    if(is.na(check)) {
        stop("components must be 3 or 5")
    }

    featureSize = dim(x)[1]

    pdata <- cbind(v1, v2)
	
    param1 <- sd(v1)
    param2 <- sd(v2)
	
    if(components == 3) {
        class <- cbind(rep(0, dim(pdata)[1]), rep(0, dim(pdata)[1]))
        class[pdata[,1]>0+param1,1] <- 2
        class[pdata[,1]<0-param1,1] <- 1
        class[pdata[,2]>0+param2,2] <- 2
        class[pdata[,2]<0-param2,2] <- 1
        discordClass <- c(2,3,4,6,7,8)
    }

    if(components == 5) {
        class <- cbind(rep(0, dim(pdata)[1]), rep(0, dim(pdata)[1]))
        class[pdata[,1]>0+param1,1] <- 2
        class[pdata[,1]<0-param1,1] <- 1
        class[pdata[,1]>0+(2*param1),1] <- 4
        class[pdata[,1]<0-(2*param1),1] <- 3
        class[pdata[,2]>0+param2,2] <- 2
        class[pdata[,2]<0-param2,2] <- 1
        class[pdata[,2]>0+(2*param2),2] <- 4
        class[pdata[,2]<0-(2*param2),2] <- 3
        concordClass <- c(1,7,13,19,25)
        discordClass <- setdiff(1:25,concordClass)
    }

    if(subsampling == TRUE) {
        subSize = nrow(x)
        if(is.null(y) == FALSE & nrow(y) < subSize) {
            subSize = nrow(y)
        }
        total_mu <- rep(0,3)
        total_sigma <- rep(0,3)
        total_nu <- rep(0,3)
        total_tau <- rep(0,3)
        total_pi <- rep(0,3)
        for(i in 1:iter) {
        # make sure pairs are independent
            rowIndex <- sample(nrow(x), subSize)
            colIndex <- sample(nrow(y), subSize)
            mat1 <- matrix(v1, nrow = nrow(x), byrow = FALSE)
            mat2 <- matrix(v2, nrow = nrow(x), byrow = FALSE)

            subSampV1 <- sapply(1:subSize, function(x) mat1[rowIndex[x], colIndex[x]])
            subSampV2 <- sapply(1:subSize, function(x) mat2[rowIndex[x], colIndex[x]])

            sub.pdata <- cbind(subSampV1, subSampV2)
            sub.class <- cbind(rep(0, subSize), rep(0, subSize))
            if(components == 3) {
                sub.class[sub.pdata[,1]>0+param1,1] <- 2
                sub.class[sub.pdata[,1]<0-param1,1] <- 1
                sub.class[sub.pdata[,2]>0+param2,2] <- 2
                sub.class[sub.pdata[,2]<0-param2,2] <- 1
            }
            if(components == 5) {
                sub.class[sub.pdata[,1]>0+param1,1] <- 2
                sub.class[sub.pdata[,1]<0-param1,1] <- 1
                sub.class[sub.pdata[,1]>0+(2*param1),1] <- 4
                sub.class[sub.pdata[,1]<0-(2*param1),1] <- 3
                sub.class[pdata[,2]>0+param2,2] <- 2
                sub.class[pdata[,2]<0-param2,2] <- 1
                sub.class[pdata[,2]>0+(2*param2),2] <- 4
                sub.class[pdata[,2]<0-(2*param2),2] <- 3
            }
            pd <- em.normal.partial.concordant(sub.pdata, sub.class, tol=0.001, restriction=0, constrain=c(0,-sd(pdata),sd(pdata)), iteration=1000)
            total_mu <- total_mu + pd$mu_sigma[1,]
            total_sigma <- total_sigma + pd$mu_sigma[2,]
            total_nu <- total_nu + pd$nu_tau[1,]
            total_tau <- total_tau + pd$nu_tau[2,]
            total_pi <- total_pi + pd$pi
        }

        mu <- total_mu/iter
        sigma <- total_sigma/iter
        nu <- total_nu/iter
        tau <- total_tau/iter
        pi <- total_pi/iter

        finalResult <- subSampleData(pdata, class, mu, sigma, nu, tau, pi)
        zTable <- finalResult$z
        classVector <- finalResult$class
    } else {
        pd <- em.normal.partial.concordant(pdata, class, tol=0.001, restriction=0, constrain=c(0,-sd(pdata),sd(pdata)), iteration=1000)
        zTable <- pd$z
        classVector <- pd$class
    }

    discordPPV <- apply(zTable, 1, function(x) sum(x[discordClass])/sum(x))

    if(is.null(y) == FALSE) {
        discordPPMatrix <- matrix(discordPPV, nrow = featureSize, byrow = FALSE)
        classMatrix <- matrix(classVector, nrow = featureSize, byrow = FALSE)
        rownames(discordPPMatrix) <- rownames(x)
        colnames(discordPPMatrix) <- rownames(y)
        rownames(classMatrix) <- rownames(x)
        colnames(classMatrix) <- rownames(y)

        vector_names <- getNames(x,y)
        names(discordPPV) <- vector_names
        names(classVector) <- vector_names
    }   

    if(is.null(y)) {
        discordPPMatrix <- matrix(NA,nrow = featureSize, ncol = featureSize)
        classMatrix <- discordPPMatrix
        diag <- lower.tri(discordPPMatrix, diag = FALSE)
        discordPPMatrix[diag] <- discordPPV
        classMatrix[diag] <- classVector
        rownames(discordPPMatrix) <- rownames(x)
        colnames(discordPPMatrix) <- rownames(x)
        rownames(classMatrix) <- rownames(x)
        colnames(classMatrix) <- rownames(x)
        vector_names <- getNames(x)
        names(discordPPV) <- vector_names
        names(classVector) <- vector_names
    }

    zTable <- t(apply(zTable, 1, function(x) x/sum(x)))
    rownames(zTable) <- vector_names

    return(list(discordPPMatrix = discordPPMatrix, discordPPVector = discordPPV, classMatrix = classMatrix, classVector = classVector, probMatrix = zTable, loglik = pd$loglik))
}

splitMADOutlier <- function(mat, filter0 = TRUE, threshold = 2) {
    if(is.matrix(mat) == FALSE) {
        mat <- as.matrix(mat)
    }
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
