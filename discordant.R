dyn.load("discordant.so")

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


em.normal.partial.concordant <- function(data, class, tol=0.000001, restriction=0, constrain=0, iteration=1000){
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
	if(restriction>0){
		if(length(constrain)==g){
			results <- .C("em_normal_partial_concordant", as.double(data[,1]), as.double(data[,2]), as.double(t(zxy)), n, pi, mu, sigma, nu, tau, g, loglik, as.double(tol), as.integer(restriction), as.integer(constrain), as.integer(iteration), convergence)
		}
		else{
			print("Error with constrain!")
			return(0)
		}
	}
	else{
		results <- .C("em_normal_partial_concordant", as.double(data[,1]), as.double(data[,2]), as.double(t(zxy)), n, pi, mu, sigma, nu, tau, g, loglik, as.double(tol), as.integer(restriction), as.integer(constrain), as.integer(iteration), convergence)
	}
	print(paste("convergence within", results[[15]], "run?", results[[16]], sep=" "))
	return(list(model="PCD", convergence=results[[16]], pi=t(array(results[[5]],dim=c(g,g))), mu_sigma=rbind(results[[6]], results[[7]]), nu_tau=rbind(results[[8]], results[[9]]), loglik=results[[11]], class=apply(array(results[[3]], dim=c(n,g*g)),1,order,decreasing=T)[1,], z=array(results[[3]], dim=c(n,g*g))))
}

fishersTrans <- function(rho) {
	r = (1+rho)/(1-rho)
	z = 0.5*log(r,base = exp(1))
	return(z)
}

createVectors <- function(data1, data2, multOmics, featureSize) {
	if(multOmics != TRUE && multOmics != FALSE) {
		stop("multOmics is not a boolean value")
	}
	
	if(dim(data1)[1] != dim(data2)[1]) {
		stop("Datasets do not have same feature size.")
	}
	
	statMatrix1 <- cor(t(data1))
	statMatrix2 <- cor(t(data2))
	if(multOmics == TRUE) {
		if(length(featureSize) == 0) {
			stop("Need an input for feature size.")
		}
		statMatrix1 <- statMatrix1[1:featureSize,(featureSize + 1):dim(data1)[1]]
		statMatrix2 <- statMatrix2[1:featureSize,(featureSize + 1):dim(data1)[1]]
	}

	statVector1 <- as.vector(statMatrix1)
	statVector2 <- as.vector(statMatrix2)
	
	if(multOmics == FALSE) {
		diagMatrix <- lower.tri(statMatrix1, diag = FALSE)
		diagVector <- as.vector(diagMatrix)
		indexVector <- which(diagVector == TRUE)
		statVector1 <- statVector1[indexVector]
		statVector2 <- statVector2[indexVector]
	}

	return(list(v1 = statVector1, v2 = statVector2))
}

discordantRun <- function(v1, v2, multOmics, transform, featureSize) {

	if(multOmics != TRUE && multOmics != FALSE) {
		stop("multOmics is not a boolean value")
	}
	if(transform != TRUE && transform != FALSE) {
		stop("transform is not a boolean value.")
	}

	if(transform == TRUE) {
		if(range(v1)[1] < -1 || range(v1)[2] > 1 || range(v2)[1] < -1 || range(v2)[2] > 1) {
			stop("correlation vectors have values less than -1 and/or greater than 1.")
		}
		v1 <- fishersTrans(v1)
		v2 <- fishersTrans(v2)
	}
	
	if(multOmics == TRUE && length(v1)%%featureSize != 0) {
		stop("featureSize is not a multiple of vector length.")
	}

	pdata <- cbind(v1, v2)
	
	param1 <- sd(v1)
	param2 <- sd(v2)

	class <- cbind(rep(0, dim(pdata)[1]), rep(0, dim(pdata)[1]))
	class[pdata[,1]>0+param1,1] <- 2
	class[pdata[,1]<0-param1,1] <- 1
	class[pdata[,2]>0+param2,2] <- 2
	class[pdata[,2]<0-param2,2] <- 1

	pd <- em.normal.partial.concordant(pdata, class, tol=0.001, restriction=0, constrain=c(0,-sd(pdata),sd(pdata)), iteration=1000)

	discordClass <- c(2,3,4,6,7,8)

	discordSum <- apply(pd$z,1,function(x) sum(x[discordClass]))
	totalSum <- apply(pd$z,1,function(x) sum(x))

	discordPPV <- discordSum/totalSum
	
	if(multOmics == TRUE) {
		discordPPMatrix <- matrix(discordPPV, nrow = featureSize, byrow = FALSE)
	}
	
	if(multOmics == FALSE) {
		tempMatrix <- matrix(NA,nrow = featureSize, ncol = featureSize)
		diagMatrix <- lower.tri(tempMatrix, diag = FALSE)
		diagVector <- as.vector(diagMatrix)
		indexVector <- which(diagVector == TRUE)
		diagVector <- rep(NA, length(diagVector))
		diagVector[indexVector] <- discordPPV
		discordPPMatrix <- matrix(diagVector, nrow = featureSize, byrow = FALSE)
	}
	
	return(list(discordPPMatrix = discordPPMatrix, class = pd$class, probMatrix = pd$z, convergence = pd$convergence, loglik = pd$loglik))
}

makeTable <- function(discordPPMatrix, multOmics, featureNames1, featureNames2 = NA) {

	if(length(featureNames1) != dim(discordPPMatrix)[1] || multOmics == TRUE && length(featureNames2) != dim(discordPPMatrix)[2] || multOmics == FALSE && length(featureNames1) != dim(discordPPMatrix)[2]) {
		stop("length of feature names does not meet dimension size")
	}
	
	if(multOmics == FALSE) {
		featureNames2 = featureNames1
	}

	outMatrix = NULL
	
	for(i in 1:dim(discordPPMatrix)[1]) {
		for(j in 1:dim(discordPPMatrix)[2]) {
			if(is.na(discordPPMatrix[i,j]) == FALSE) {
				row <- c(featureNames1[i], featureNames2[j], discordPPMatrix[i,j])
				outMatrix <- rbind(outMatrix, row)
			}
		}
	}
	
				
	
	return(outMatrix)
}
