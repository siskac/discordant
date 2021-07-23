# Applied from Lai et al 2007 Bioinformatics.
#
# Code written by Charlotte Siska and Max McGrath
# R functions fishersTrans, createVectors, discordantRun and makeTable
# Email: charlotte.siska@ucdenver.edu
#
# Code written by Yinglei Lai
# C code and R functions unmap and em.normal.partial.concordant
# E-mail: ylai@gwu.edu
#
# Code written by Fang Huaying
# R function for SparCC, adapted from https://bitbucket.org/yonatanf/sparcc
# Email: hyfang@pku.edu.cn

#' @export
discordantRun <- function(v1, v2, x, y = NULL, transform = TRUE, 
                          subsampling = FALSE, subSize = dim(x)[1], iter = 100, 
                          components = 3) {
  
    .checkDiscordantInputs(v1, v2, x, y, transform, subsampling, subSize, iter, 
                           components)
    
    if (transform) {
        v1 <- fishersTrans(v1)
        v2 <- fishersTrans(v2)
    }
  
    x <- exprs(x)
    if (!is.null(y)) { y <- exprs(y) }
    featureSize = dim(x)[1]
    
    pdata <- cbind(v1, v2)
    param1 <- sd(v1)
    param2 <- sd(v2)
    class <- cbind(.assignClass(v1, param1, components),
                   .assignClass(v2, param2, components))
    
    if (components == 3) {
        discordClass <- c(2,3,4,6,7,8)
    } else {
        discordClass <- setdiff(1:25, c(1, 7, 13, 19, 25))
    }
    
    if (subsampling) {
        subSize <- min(nrow(x), nrow(y))
        total_mu <- total_sigma <- total_nu <- 
          total_tau <- total_pi <- rep(0, components) 
        # NOTE: changed this from default 3 to components... need to figure out 
        #   if that's right
        
        for(i in 1:iter) {
            # make sure pairs are independent
            rowIndex <- sample(nrow(x), subSize)
            colIndex <- sample(nrow(y), subSize)
            mat1 <- matrix(v1, nrow = nrow(x), byrow = FALSE)
            mat2 <- matrix(v2, nrow = nrow(x), byrow = FALSE)
            
            subSampV1 <- sapply(1:subSize, function(x) mat1[rowIndex[x], 
                                                            colIndex[x]])
            subSampV2 <- sapply(1:subSize, function(x) mat2[rowIndex[x], 
                                                            colIndex[x]])
            
            sub.pdata <- cbind(subSampV1, subSampV2)
            sub.class <- cbind(.assignClass(subSampV1, param1, components),
                               .assignClass(subSampV2, param2, components))
            
            pd <- em.normal.partial.concordant(sub.pdata, sub.class, components)
            total_mu <- total_mu + pd$mu_sigma[1,]
            total_sigma <- total_sigma + pd$mu_sigma[2,]
            total_nu <- total_nu + pd$nu_tau[1,]
            total_tau <- total_tau + pd$nu_tau[2,]
            total_pi <- total_pi + pd$pi
      }
      
      mu <- total_mu / iter
      sigma <- total_sigma / iter
      nu <- total_nu / iter
      tau <- total_tau / iter
      pi <- total_pi / iter
      
      finalResult <- subSampleData(pdata, class, mu, sigma, nu, tau, pi, 
                                   components)
      zTable <- finalResult$z
      classVector <- finalResult$class
    } else {
        pd <- em.normal.partial.concordant(pdata, class, components)
        zTable <- pd$z
        classVector <- pd$class
    }
    
    discordPPV <- apply(zTable, 1, function(x) sum(x[discordClass])/sum(x))
    
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
    } else {
        discordPPMatrix <- matrix(discordPPV, nrow = featureSize, 
                                  byrow = FALSE)
        classMatrix <- matrix(classVector, nrow = featureSize, byrow = FALSE)
        rownames(discordPPMatrix) <- rownames(x)
        colnames(discordPPMatrix) <- rownames(y)
        rownames(classMatrix) <- rownames(x)
        colnames(classMatrix) <- rownames(y)
        
        vector_names <- getNames(x,y)
        names(discordPPV) <- vector_names
        names(classVector) <- vector_names
    }
    
    zTable <- t(apply(zTable, 1, function(x) x/sum(x)))
    rownames(zTable) <- vector_names
    
    return(list(discordPPMatrix = discordPPMatrix, discordPPVector = discordPPV,
                classMatrix = classMatrix, classVector = classVector, 
                probMatrix = zTable, loglik = pd$loglik))
}

em.normal.partial.concordant <- function(data, class, components) {
    tol <- 0.001
    restriction <- 0
    constrain <- 0
    iteration <- 1000
    n <- as.integer(dim(data)[1])
    g <- as.integer(nlevels(as.factor(class)))

    yl.outer <- function(k, zx, zy){
        return( c(zx[k,] %o% zy[k,]) )
    }

    zx <- unmap(class[,1], components = components)
    zy <- unmap(class[,2], components = components)
    zxy <- sapply(1:dim(zx)[1], yl.outer, zx, zy)

    pi <- double(g*g)
    mu <- double(g)
    sigma <- double(g)
    nu <- double(g)
    tau <- double(g)
    loglik <- double(1)
    convergence <- integer(1)
    
    results <- em_normal_partial_concordant_cpp(as.double(data[,1]), 
                                                as.double(data[,2]), 
                                                as.double(t(zxy)), n, pi, mu, 
                                                sigma, nu, tau, g, loglik, 
                                                as.double(tol), 
                                                as.integer(restriction), 
                                                as.integer(constrain), 
                                                as.integer(iteration), 
                                                convergence)
    
    return(list(model = "PCD", 
                convergence = results[[16]],
                pi = t(array(results[[5]], dim=c(g,g))), 
                mu_sigma = rbind(results[[6]], results[[7]]), 
                nu_tau = rbind(results[[8]], results[[9]]),
                loglik = results[[11]], 
                class = apply(array(results[[3]], dim = c(n,g*g)), 
                              1, order, decreasing = TRUE)[1,], 
                z = array(results[[3]], dim = c(n, g*g))))
}


# Internal function to assign class to vector based on number of components and 
#   each elements distance from zero
#' @importFrom dplyr case_when
.assignClass <- function(x, param, components) {
  if (components == 3) {
    rtn <- case_when(x < -param ~ 1, x > param ~ 2, TRUE ~ 0)
  } else {
    rtn <- case_when(x < -2*param ~ 3, x > 2*param ~ 4, x < -param ~ 1,
                     x > param ~ 2, TRUE ~ 0)
  }
  return(rtn)
}

# Internal function to validate user inputs for discordantRun()
.checkDiscordantInputs <- function(v1, v2, x, y, transform, 
                                   subsampling, subSize, iter, 
                                   components) {
  
  if (!is(x, "ExpressionSet") || (!is.null(y) && !is(y, "ExpressionSet"))) {
    stop("x and y (if present) must be type ExpressionSet")
  }
  
  # Need to double check if this is true w/ Katerina
  if (is.null(y) && subsampling) {
    stop("y cannot be NULL if subsampling is TRUE")
  }
  
  if (!(components %in% c(3, 5))) {
    stop ("components must be equal to 3 or 5")
  }
  
  if (transform && (range(v1)[1] < -1 || range(v1)[2] > 1 || 
                    range(v2)[1] < -1 || range(v2)[2] > 1)) {
    stop ("correlation vectors have values less than -1 and/or greater than 1.")
  }
}