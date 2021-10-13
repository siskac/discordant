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

#' Run Discordant Algorithm
#' 
#' Runs discordant algorithm on two vectors of correlation coefficients.
#' 
#' @param v1 Vector of correlation coefficients in group 1
#' @param v2 Vector of correlation coefficients in group 2
#' @param x ExpressionSet of -omics data
#' @param y ExpressionSet of -omics data, induces dual -omics analysis
#' @param transform If TRUE v1 and v2 will be Fisher transformed
#' @param subsampling If TRUE subsampling will be run
#' @param subSize Indicates how many feature pairs to be used for subsampling. 
#' Default is the feature size in x
#' @param iter Number of iterations for subsampling. Default is 100
#' @param components Number of components in mixture model.
#' 
#' @return
#' \describe{
#'   \item{discordPPVector}{Vector of differentially correlated posterior 
#'   probabilities.}
#'   \item{discordPPMatrix}{Matrix of differentially correlated posterior 
#'   probabilities where rows and columns reflect features}
#'   \item{classVector}{Vector of classes that have the highest posterior 
#'   probability}
#'   \item{classMatrix}{Matrix of classes that have hte highest posterior 
#'   probability where rows and columns reflect features}
#'   \item{probMatrix}{Matrix of posterior probabilities where rows are each 
#'   molecular feature pair and columns are nine different classes}
#'   \item{loglik}{Final log likelihood}
#' }
#' @details 
#' The discordant algorithm is based on a Gaussian mixture model. If there are 
#' three components, correlation coefficients are clustered into negative 
#' correlations (-), positive correlations (+) and no correlation (0). If there 
#' are five components, then there are two more classes for very negative 
#' correlation (--) and very positive correlations (++). All possible 
#' combinations for these components are made into classes. If there are three 
#' components, there are 9 classes. If there are five components, there are 25 
#' classes.
#' 
#' The posterior probabilities for each class are generated and outputted into 
#' the value probMatrix. The value probMatrix is a matrix where each column is a
#'  class and each row is a feature pair. The values discordPPVector and 
#'  discordPPMatrix are the summed differential correlation posterior 
#'  probability for each feature pair. The values classVector and classMatrix 
#'  are the class with the highest posterior probability for each feature pair.
#' @references 
#' Siska C, Bowler R and Kechris K. The Discordant Method: A Novel Approach for 
#' Differential Correlation (2015), Bioinformatics. 32 (5): 690-696.
#' 
#' Lai Y, Zhang F, Nayak TK, Modarres R, Lee NH and McCaffrey TA. Concordant 
#' integrative gene set enrichment analysis of multiple large-scale two-sample 
#' expression data sets. (2014) BMC Genomics 15, S6.
#' 
#' Lai Y, Adam B-l, Podolsky R, She J-X. A mixture model approach to the tests 
#' of concordance and discordancd between two large-scale experiments with two 
#' sample groups. (2007) Bioinformatics 23, 1243-1250.
#' 
#' @author Charlotte Siska \email{siska.charlotte@@gmail.com}
#' @author Max McGrath \email{max.mcgrath@@ucdenver.edu}
#' 
#' @examples
#' # Load Data
#' data(TCGA_GBM_miRNA_microarray)
#' data(TCGA_GBM_transcript_microarray)
#' print(colnames(TCGA_GBM_transcript_microarray)) # look at groups
#' groups <- c(rep(1,10), rep(2,20))
#' 
#' ## DC analysis on only transcripts pairs
#' 
#' vectors <- createVectors(TCGA_GBM_transcript_microarray, 
#'                          groups = groups)
#' result <- discordantRun(vectors$v1, vectors$v2, 
#'                         TCGA_GBM_transcript_microarray)
#' 
#' ## DC analysis on miRNA-transcript pairs
#' 
#' vectors <- createVectors(TCGA_GBM_transcript_microarray, 
#'                          TCGA_GBM_miRNA_microarray, groups = groups, 
#'                          cor.method = c("pearson"))
#' result <- discordantRun(vectors$v1, vectors$v2, 
#'                         TCGA_GBM_transcript_microarray, 
#'                        TCGA_GBM_miRNA_microarray)
#' @export
discordantRun <- function(v1, v2, x, y = NULL, transform = TRUE, 
                          subsampling = FALSE, subSize = NULL, iter = 100, 
                          components = 3) {
  
    .checkDiscordantInputs(v1, v2, x, y, transform, subsampling, subSize, iter, 
                           components)
    
    if (transform) {
        v1 <- fishersTrans(v1)
        v2 <- fishersTrans(v2)
    }
    x <- exprs(x)
    if (!is.null(y)) { y <- exprs(y) }
    pdata <- cbind(v1, v2)
    param1 <- sd(v1)
    param2 <- sd(v2)
    class <- cbind(.assignClass(v1, param1, components),
                   .assignClass(v2, param2, components))
    
    if (subsampling) {
        subSize <- .setSubSize(x, y, subSize)
        total_mu <- total_sigma <- total_nu <- 
          total_tau <- total_pi <- rep(0, components)
        
        i <- 0
        repeats <- 0
        
        while(i < iter && repeats < floor(iter * .1)) {
            subSamples <- .createSubsamples(x, y, v1, v2, subSize)
            sub.pdata <- cbind(subSamples$v1, subSamples$v2)
            sub.class <- cbind(.assignClass(subSamples$v1, param1, components),
                               .assignClass(subSamples$v2, param2, components))
            
            pd <- tryCatch({em.normal.partial.concordant(sub.pdata, sub.class, 
                                                   components)},
                error = function(unused) return(NULL))
            
            if(is.null(pd)) { 
              repeats <- repeats + 1 
            } else {
              i <- i + 1
              total_mu <- total_mu + pd$mu_sigma[1,]
              total_sigma <- total_sigma + pd$mu_sigma[2,]
              total_nu <- total_nu + pd$nu_tau[1,]
              total_tau <- total_tau + pd$nu_tau[2,]
              total_pi <- total_pi + pd$pi
            }
        }
        
        if (repeats >= floor(iter * .1)) {
          stop(paste0("\nInsufficient data for subsampling. Increase number of",
               "\nfeatures, reduce numberof components used, or increase", 
               "\nsubSize if not at default value. Alternatively, set",
               "\nsubsampling=FALSE."))
        }
      
      mu <- total_mu / iter
      sigma <- total_sigma / iter
      nu <- total_nu / iter
      tau <- total_tau / iter
      pi <- total_pi / iter
      
      finalResult <- .subSampleData(pdata, class, mu, sigma, nu, tau, pi, 
                                   components)
      zTable <- finalResult$z
      classVector <- finalResult$class
    } else {
        pd <- tryCatch({em.normal.partial.concordant(pdata, class, components)},
                       error = function(unused) return(NULL))
        if (is.null(pd)) {
          stop(
            paste0("\nInsufficient data for component estimation. Increase",
                   "\nnumber of features or reduce number of components used."))
        }
        zTable <- pd$z
        classVector <- pd$class
    }
    
    rtn <- .prepareOutput(x, y, pd, zTable, classVector, components)
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

    zx <- .unmap(class[,1], components = components)
    zy <- .unmap(class[,2], components = components)
    # .checkForMissingComponents(zx, zy)
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

.prepareOutput <- function(x, y, pd, zTable, classVector, components) {
  
    if (components == 3) {
        discordClass <- c(2,3,4,6,7,8)
    } else {
        discordClass <- setdiff(1:25, c(1, 7, 13, 19, 25))
    }
    featureSize = dim(x)[1]
    discordPPV <- apply(zTable, 1, function(x) sum(x[discordClass]) / sum(x))
    
    if(is.null(y)) {
        discordPPMatrix <- matrix(NA, nrow = featureSize, ncol = featureSize)
        classMatrix <- discordPPMatrix
        diag <- lower.tri(discordPPMatrix, diag = FALSE)
        discordPPMatrix[diag] <- discordPPV
        classMatrix[diag] <- classVector
        colnames(discordPPMatrix) <- rownames(x)
        colnames(classMatrix) <- rownames(x)
        vector_names <- .getNames(x)
    } else {
        discordPPMatrix <- matrix(discordPPV, nrow = featureSize, 
                                  byrow = FALSE)
        classMatrix <- matrix(classVector, nrow = featureSize, byrow = FALSE)
        colnames(discordPPMatrix) <- rownames(y)
        colnames(classMatrix) <- rownames(y)
        vector_names <- .getNames(x,y)
    }

    rownames(discordPPMatrix) <- rownames(x)
    rownames(classMatrix) <- rownames(x)
    names(discordPPV) <- vector_names
    names(classVector) <- vector_names
    
    zTable <- t(apply(zTable, 1, function(x) x / sum(x)))
    rownames(zTable) <- vector_names
    
    return(list(discordPPMatrix = discordPPMatrix, discordPPVector = discordPPV,
                classMatrix = classMatrix, classVector = classVector, 
                probMatrix = zTable, loglik = pd$loglik))
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

# Internal function to independently subsample data according to whether the
#   analysis is within -omics or between -omics. 
.createSubsamples <- function(x, y = NULL, v1, v2, subSize) {
    if(is.null(y)) {
        sampleNames <- rownames(x)
        nameSet1 <- sample(sampleNames, subSize)
        nameSet2 <- sample(setdiff(sampleNames, nameSet1), subSize)
        nameSetA <- paste0(nameSet1, "_", nameSet2)
        nameSetB <- paste0(nameSet2, "_", nameSet1)
        subSampV1 <- ifelse(is.na(v1[nameSetA]), v1[nameSetB], v1[nameSetA])
        subSampV2 <- ifelse(is.na(v2[nameSetA]), v2[nameSetB], v2[nameSetA])
    } else {
        # make sure pairs are independent
        rowIndex <- sample(nrow(x), subSize)
        colIndex <- sample(nrow(y), subSize)
        mat1 <- matrix(v1, nrow = nrow(x), byrow = FALSE)
        mat2 <- matrix(v2, nrow = nrow(x), byrow = FALSE)
        subSampV1 <- sapply(1:subSize, function(x) mat1[rowIndex[x], 
                                                        colIndex[x]])
        subSampV2 <- sapply(1:subSize, function(x) mat2[rowIndex[x], 
                                                        colIndex[x]])
    }
    return(list(v1 = subSampV1,
                v2 = subSampV2))
}

# Internal function to validate user inputs for discordantRun()
#' @importFrom methods is
.checkDiscordantInputs <- function(v1, v2, x, y, transform, 
                                   subsampling, subSize, iter, 
                                   components) {
  
  if (!is(x, "ExpressionSet") || (!is.null(y) && !is(y, "ExpressionSet"))) {
    stop("x and y (if present) must be type ExpressionSet")
  }
  
  if (!(components %in% c(3, 5))) {
    stop ("Components must be equal to 3 or 5")
  }
  
  if (transform && (range(v1)[1] < -1 || range(v1)[2] > 1 || 
                    range(v2)[1] < -1 || range(v2)[2] > 1)) {
    stop ("Correlation vectors have values less than -1 and/or greater than 1.")
  }
}

# Internal function that checks whether all types of component are present
#   in given vectors. If a certain component is not present, we run into a
#   divide-by-zero error that crashes R
# .checkForMissingComponents <- function(zx, zy) {
#   sumZx <- colSums(zx)
#   sumZy <- colSums(zy)
#   if(any(c(sumZx, sumZy) == 0)) {
#     stop("Insufficient data for component estimation. Increase number of 
#          features or reduce number of components used. If subsampling=TRUE,
#          you may need set subsampling=FALSE.")
#   }
# }

# Internal function to set sub size based on data
.setSubSize <- function(x, y, subSize) {
  if (is.null(y)) {
    if (is.null(subSize)) {
      subSize <- floor(length(rownames(x)) / 2)
    } else if (subSize > floor(length(rownames(x)) / 2)) {
      subSize <- floor(length(rownames(x)) / 2)
      warning(paste0("subSize argument too large. Using subSize ", subSize,
                     " See vignette for more information."))
    }
  } else {
    if (is.null(subSize)) {
      subSize <- min(nrow(x), nrow(y))
    } else if (subSize > min(nrow(x), nrow(y))) {
      subSize <- min(nrow(x), nrow(y))
      warning(paste0("subSize argument to large. Using subSize ", subSize,
                     " See vignette for more information."))
    }
  }
  return(subSize)
}