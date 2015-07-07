##############################################################	
##############################################################	
##############################################################	

pilotQV <- function(x, y, lambda.grid=NULL, lambda=NULL, weights = NULL){
	ifelse(is.matrix(x), ncolx <- ncol(x), ncolx <- 1)
	pilot <- Tps(x,y, lambda=lambda, weights =weights)	
    
    if(is.null(lambda.grid))lambda.grid <- pilot$gcv.grid[,1]
    
	nh <- length(lambda.grid)
	p_fitted <- pilot$fitted.values
	p_sigmahat <- pilot$shat.GCV
	
	V <- diag(pilot$matrices$np)
	V[(pilot$matrices$nt+1):pilot$matrices$np,(pilot$matrices$nt+1):pilot$matrices$np] <- pilot$matrices$V
	QV <- qr.qy(pilot$matrices$qr.T, V)
	
	xM <- pilot$xM
	yM <- pilot$yM
	weightsM <- pilot$weightsM
	uniquerows <- pilot$uniquerows
	dupinfo <- ifelse(pilot$weightsM < length(y), TRUE, FALSE)
	
	return(list(matrices = pilot$matrices, QV = QV, lambda.grid = lambda.grid, 
				p_fitted = p_fitted, p_sigmahat = p_sigmahat, ncolx = ncolx, 
				eff.df = pilot$eff.df, p_lambda = pilot$lambda, y=y, x=x, 
				xM=xM, yM=yM, uniquerows=uniquerows, weightsM = weightsM, dupinfo=dupinfo))
}

##############################################################	
##############################################################	
##############################################################	

localriskR <- function(lambda, pilot, rmethod=c("estimated", "true", "corrected"), pilotf=NULL, pilotsigma=NULL){
	weights2 <- sqrt(pilot$weightsM)
	D0 <- pilot$matrices$D[(pilot$matrices$nt+1):pilot$matrices$np]
	D <- c(c(rep(1, pilot$matrices$nt), 1/(1+D0*lambda)))
	ID <- c(c(rep(0, pilot$matrices$nt), D0*lambda/(1+D0*lambda)))
	
	QVD2 <- t(apply(pilot$QV, 1, function(x,d) x*d, d=D^2) )
	QVID <- t(apply(pilot$QV, 1, function(x,d) x*d, d=ID) )
	#QVD <- t(apply(pilot$QV, 1, function(x,d) x*d, d=D) )
	trA2 <- c(apply((QVD2) * pilot$QV, 1, sum))
	IA <- tcrossprod(QVID, pilot$QV)
	IA <- apply(IA, 2, function(x, wt2) x/wt2, wt2=weights2)
	IA <- t(apply(IA, 1, function(x, wt2) x*wt2, wt2=weights2))
	#A <- tcrossprod(QVD, pilot$QV)
	if(rmethod[1]=="estimated"){
		
		pilotf <- c(pilot$p_fitted) 
		pilotsigma <- pilot$p_sigmahat
		out <- c((IA%*%pilotf)^2 + (pilotsigma^2)*trA2)
	}
	if(rmethod[1]=="true"){
		if(is.null(pilotf) | is.null(pilotsigma)){
			if(is.null(pilot$tf)) {
				warning("True f and variance must be...")
				return("Error")	
			}else{
				pilotf <- pilot$tf
				pilotsigma <- pilot$tsigma
			}
		}
		out <- c((IA%*%pilotf)^2 + (pilotsigma^2)*trA2)
	}
	
	if(rmethod[1]=="corrected"){
		
		QVID2 <- t(apply(pilot$QV, 1, function(x,d) x*d, d=ID^2) )
		trIA2 <- c(apply((QVID2) * pilot$QV, 1, sum))	
		pilotf <- pilot$y
		pilotsigma <- pilot$p_sigmahat
		out <-c((IA%*%y)^2 - (pilotsigma^2)*trIA2 + (pilotsigma^2)*trA2)
	}
	out
}

#R --arch=x86_64 CMD SHLIB -dynamiclib -O3 -fopenmp localrisk.f90 localfit.f90  -o localmethod64.so 
#R --arch=x86_64 CMD SHLIB -dynamiclib -O3 localrisk.f90 localfit.f90 -o localmethod64.so 
#rm *.o
#R --arch=i386 CMD SHLIB -dynamiclib -O3 localrisk.f90 localfit.f90 -o localmethod32.so 
#R --arch=i386 CMD SHLIB --libs-only -dynamiclib -O3 localrisk.f90 localfit.f90 -o localmethod32.so 

path <- "~/Dropbox/Research/CP/"
if(Sys.info()[1] == "Linux"){
  dyn.load(paste(path, "localmethod64_unix.so", sep=""))
} else{
  #if(.Platform$r_arch == "x86_64"){
    dyn.load(paste(path, "localmethod64.so", sep=""))
  #}else{
  #  dyn.load(paste(path, "localmethod32.so", sep=""))
  #}
}





localrisk <- function(pilot, rmethod=c("estimated", "true", "corrected"), pilotf=NULL, pilotsigma=NULL){
	
	yM <- pilot$y[pilot$uniquerows]
	n <- length(yM)
	
	
	nh <- length(pilot$lambda.grid)
	
	if(rmethod[1]=="estimated"){
		
		if(is.null(pilotf)){
			if(length(pilot$p_fitted) != n){
				pilotf <- c(pilot$p_fitted)[pilot$uniquerows]
			}else{
				pilotf <- c(pilot$p_fitted)
			}
		}
		if(is.null(pilotsigma)){
			pilotsigma <- pilot$p_sigmahat
		}
		
		out <- .Fortran("localrisk1", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
		dim(out) <- c(n, nh)
	}
	
	if(rmethod[1]=="true"){
		if(is.null(pilotf) | is.null(pilotsigma)){
			if(is.null(pilot$tf)) {
				warning("True f and variance must be...")
				return("Error")	
			}else{
				if(length(pilot$tf) != n){
					pilotf <- pilot$tf[pilot$uniquerows]
				}else{
					pilotf <- pilot$tf
				}
				
				pilotsigma <- pilot$tsigma
			}
		}
		
		out <- .Fortran("localrisk1", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
		dim(out) <- c(n, nh)
	}
	
	if(rmethod[1]=="corrected"){
		if(is.null(pilotf)){
			pilotf <- yM
		}
		if(is.null(pilotsigma)){
			pilotsigma <- pilot$p_sigmahat
		}
		
		out <- .Fortran("localrisk2", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
		dim(out) <- c(n, nh)
	}
    if(!is.null(dim(pilot$x))){
        xM = pilot$x[pilot$uniquerows, ]
    }else{
        xM = pilot$x[pilot$uniquerows]
    }
	
	return(list(risk=out, xM=xM, lambda.grid=pilot$lambda.grid))
}

##############################################################	
##############################################################	
##############################################################	

#blockrisk <- function(lambda, pilot, rmethod=c("estimated", "true", "corrected"), bind){
#	out <- localrisk(lambda, pilot, rmethod)
#	return(sum(out[bind])/length(bind))
#}
movingriskR <- function(lrisk, newx, rval){
    ndim <- ncol(newx)
    if(is.null(ndim)) ndim <- 1
    
    x <- lrisk$xM
    lambda.grid <- lrisk$lambda.grid
    nh <- length(lambda.grid)
    
    if(ndim==1){
        nn <- length(newx)
        brisk <- matrix(NA, nn, nh)
        
        for(i in 1:nn){
            ind <- which( abs(newx[i] - x) < rval)
            if(length(ind) > 1){
                brisk[i,] <- apply(lrisk$risk[ind, ], 2, mean)   
            }else{
                brisk[i,] <- lrisk$risk[ind, ]
            }
        }
    } else{
        stop("The univariate case is supported only. \n" )
    }
    return(brisk)
}

movinglambda <- function(lrisk, newx, rval=0.05, boundary=c("symmetric", "none")){
    ndim <- ncol(newx)
    if(is.null(ndim)) ndim <- 1
    x <- lrisk$xM
    lambda.grid <- lrisk$lambda.grid
    lrisk <- lrisk$risk
    nh <- length(lambda.grid)
    
    if(ndim==1){
        
        n <- length(x)
        nn <- length(newx)
        if(boundary[1] == "none"){
            minloc <- .Fortran("movingrisk", as.double(lrisk), as.integer(n), as.integer(nh), 
                    as.double(x), as.double(newx), as.integer(nn), 
                    as.double(rval), movrisk=double(nn), minloc=integer(nn))
        }else{
            minloc <- .Fortran("movingrisk3", as.double(lrisk), as.integer(n), as.integer(nh), 
                    as.double(x), as.double(newx), as.integer(nn), 
                    as.double(rval), movrisk=double(nn), minloc=integer(nn))        
        }
        inds <- which(minloc$minloc == 0)
        minloc$minloc[inds] <- NA
        minloc$movrisk[inds] <- NA
        #if(length(inds) > 0){
        #    for(i in 1:length(inds)){
        #        newxi <- newx[i]
        #        tmpx <- abs(x-newxi)
        #        tmpind <- order(tmpx)[1:3]
        #        tmprisk <- apply(lrisk[tmpind, ], 2, mean)
        #        minloc$minloc[inds[i]] <- which.min(tmprisk)
        #        minloc$movrisk[inds[i]] <- min(tmprisk)
        #    }
        #    
        #}
        
        lambda <- lambda.grid[minloc$minloc]
        risk <- minloc$movrisk
        
    }else{
        n <- nrow(x)
        nn <- nrow(newx)
        if(is.null(nn)) stop("The number of col of newx must be >= 2.")
        minloc <- .Fortran("movingrisk2", as.double(lrisk), as.integer(n), as.integer(nh), 
                as.double(as.matrix(x)), as.double(as.matrix(newx)), as.integer(nn), 
                as.double(rval), movrisk=double(nn), minloc=integer(nn))
        inds <- which(minloc$minloc == 0)
        minloc$minloc[inds] <- NA
        minloc$movrisk[inds] <- NA
        
        #if(length(inds) > 0){
        #    for(i in 1:length(inds)){
        #        newxi <- newx[i,]
        #        tmpx <- as.vector(apply(t(t(x) - unlist(newxi))^2, 1, sum))
        #        tmpind <- order(tmpx)[1:6]
        #        tmprisk <- apply(lrisk[tmpind, ], 2, mean)
        #        minloc$minloc[inds[i]] <- which.min(tmprisk)
        #        minloc$movrisk[inds[i]] <- min(tmprisk)
        #    }
        #    
        #}
        lambda <- lambda.grid[minloc$minloc]
        risk <- minloc$movrisk 
        
    }
	
	return(list(lambda = lambda, risk=risk))
}


localsmooth <- function(pilot, lambda){
	y <- pilot$yM
	n <- length(y)
	
	locfit <- .Fortran("localfit", as.double(lambda), as.double(y), as.integer(n), 
					   as.double(pilot$matrices$D), as.double(pilot$QV),
					   fitted.values = double(n), diagA = double(n)) 
	
	diag_A <- locfit$diagA
	f_pilot <- locfit$fitted.values
	sigma_pilot <- sqrt(sum((y - f_pilot)^2)/(n - sum(diag_A)))
	out <- list(diagA = diag_A, fitted.values=f_pilot, sigmahat = sigma_pilot)
	return(out)
	
}



improve_local <- function(pilot, newx, rval, sparameter=NULL, maxiter=10, theta=.15, boundary=c("symmetric", "none"), termval=0.1^4){
	lambda.grid <- pilot$lambda.grid
    nh <- length(lambda.grid)
    
    ndim <- dim(pilot$x)
    if(!is.null(ndim)){
        xx <- pilot$x[pilot$uniquerows,]
        yy <- c(pilot$y)[pilot$uniquerows]
        nn <- nrow(newx)
    }else{
        xx <- pilot$x[pilot$uniquerows]
        yy <- pilot$y[pilot$uniquerows]
        nn <- length(newx)
    }
    
	n <- length(yy)
    #pilot esitmator
    smovlam <- matrix(NA, n, maxiter+1)
	predy <- matrix(NA, n, maxiter+1)
	predsigma <- rep(NA, maxiter+1)
	diagA <- matrix(NA, n, maxiter+1)
	eff.df <- rep(NA, maxiter+1)
	
	movrisk <- matrix(NA, nn, maxiter+1)
	movlambda <- matrix(NA, nn, maxiter+1)
	
	# pilot estimator
	smovlam[,1] <- pilot$p_lambda # rep(pilot$p_lambda, n)
	predy[,1] <- c(pilot$p_fitted)[pilot$uniquerows]
	predsigma[1] <- pilot$p_sigmahat
	eff.df[1] <- pilot$eff.df
	diagA[,1] <- 1/(1+pilot$p_lambda*pilot$matrice$D)
	
	# iterative moving risk and lambda
	movrisk[,1] <- NA
	movlambda[,1] <- NA # rep(pilot$p_lambda, nn)
	
	error <- Inf
	old_predlambda <- rep(pilot$p_lambda, n)
	
	iter <- 0
	############
	aa <- list()
	############
	
	while(error > termval){
	   iter <- iter + 1
		if(iter > maxiter) {
			warning(paste("The number of iteration is greater than ", maxiter, ". \n", sep=""))
			return(out)
		}
		#	cat(iter, "\n")
		
		
		elrisk <- localrisk(pilot, rmethod="estimated")

		if(iter==1) risk0 <- elrisk

		elambdafit <- movinglambda(elrisk, newx, rval, boundary[1])
		
		#aa[[iter]] <- movingriskR(elrisk, newx, rval)
		moverisk <- elambdafit$risk
		elambda <- elambdafit$lambda
		##predrisk <- cbind(predrisk, moverisk)
		movrisk[, iter+1] <- moverisk
		##movlambda <- cbind(movlambda, elambda)
		movlambda[, iter+1] <- elambda
		
        if(!is.null(ndim)){
            if(is.null(sparameter)) stop("The sparameter is NULL \n")
            fite <- fastTps(newx[!is.na(elambda),], -log(elambda[!is.na(elambda)]), lambda=sparameter, theta=theta)
        }else{
            
            tmpx <- newx[!is.na(elambda)]
            #tmpy <- -log(elambda[!is.na(elambda)])
            tmpy <- log(elambda[!is.na(elambda)])
           
            fite <- Tps(tmpx, tmpy, lambda=sparameter)
        #    cat(fite$lambda, "\n")
        #    fite <- Tps(newx, -log(elambda))
        #    fite <- sreg(newx, -log(elambda), lambda=sparameter)
		    }
        #predlambda <- c(exp(-predict(fite, xx)))
        predlambda <- c(exp(predict(fite, xx)))
		
		locfit <- .Fortran("localfit", as.double(predlambda), as.double(yy), as.integer(n), 
						   as.double(pilot$matrices$D), as.double(pilot$QV),
						   fitted.values = double(n), diagA = double(n)) 
		
		diag_A <- locfit$diagA
		f_pilot <- c(locfit$fitted.values)
		sigma_pilot <- sqrt(sum((yy - f_pilot)^2)/(n - sum(diag_A)))
        #sigma_pilot <- 0.0463065
		pilot$p_fitted <- f_pilot
		pilot$p_sigmahat <- sigma_pilot
		
######## History information
		##smovlam <- cbind(smovlam, predlambda)
		smovlam[, iter+1] <- predlambda
		##predy <- cbind(predy, f_pilot)
		predy[, iter+1] <- f_pilot
		#predsigma <- c(predsigma, sigma_pilot)
		predsigma[iter+1] <- sigma_pilot
		##diagA <- cbind(diagA, diag_A)
		diagA[,iter+1] <- diag_A
		#eff.df <- c(eff.df, sum(diag_A))
		eff.df[iter+1] <- sum(diag_A)
######## 		
		error <- sum((log(old_predlambda) - log(predlambda))^2)/n
		old_predlambda <- predlambda
		
		out <- list(smovlambda=smovlam[,1:(iter+1)], predy=predy[,1:(iter+1)], predsigma = predsigma[1:(iter+1)], 
		            movrisk=movrisk[,1:(iter+1)], movlambda = movlambda[,1:(iter+1)], 
		            finlrisk=elrisk, inilrisk=risk0, diagA=diagA[,1:(iter+1)], eff.df=eff.df[1:(iter+1)], iteration =iter, boundary=boundary[1])
	}
	#bb <<- aa
	return(out)
}



fullA <- function(pilot){
	
	yM <- pilot$y[pilot$uniquerows]
	n <- length(yM)
	
	
	nh <- length(pilot$lambda.grid)
	
	
    out <- .Fortran("iatra2", as.double(pilot$lambda.grid), as.integer(nh),
				as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
				IAall=double(n*n*nh), Tra2all=double(n*nh))[c("IAall", "Tra2all")]
    dim(out$IAall) <- c(n*n, nh)
	dim(out$Tra2all) <- c(n, nh)
	
	return(out)
}



newlocalrisk <- function(pilot, ias=NULL, tra2s=NULL, rmethod=c("estimated", "true", "corrected"), pilotf=NULL, pilotsigma=NULL){
	
	yM <- pilot$y[pilot$uniquerows]
	n <- length(yM)
	
	
	nh <- length(pilot$lambda.grid)
	
	if(rmethod[1]=="estimated"){
		
		if(is.null(pilotf)){
			if(length(pilot$p_fitted) != n){
				pilotf <- c(pilot$p_fitted)[pilot$uniquerows]
			}else{
				pilotf <- c(pilot$p_fitted)
			}
		}
		if(is.null(pilotsigma)){
			pilotsigma <- pilot$p_sigmahat
		}
		
		if ( is.null(ias) | is.null(tra2s) ){
            out <- .Fortran("localrisk1", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        } else {
            out <- .Fortran("localriskall", as.integer(n), as.integer(nh),
						as.double(pilotf), as.double(pilotsigma), as.double(ias), as.double(tra2s),
						risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        
        }
	}
	
	if(rmethod[1]=="true"){
		if(is.null(pilotf) | is.null(pilotsigma)){
			if(is.null(pilot$tf)) {
				warning("True f and variance must be...")
				return("Error")	
			}else{
				if(length(pilot$tf) != n){
					pilotf <- pilot$tf[pilot$uniquerows]
				}else{
					pilotf <- pilot$tf
				}
				
				pilotsigma <- pilot$tsigma
			}
		}
		if ( is.null(ias) | is.null(tra2s) ){
            out <- .Fortran("localrisk1", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        } else {
            out <- .Fortran("localriskall", as.integer(n), as.integer(nh),
						as.double(pilotf), as.double(pilotsigma), as.double(ias), as.double(tra2s),
						risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        
        }        
	}
	
	if(rmethod[1]=="corrected"){
		if(is.null(pilotf)){
			pilotf <- yM
		}
		if(is.null(pilotsigma)){
			pilotsigma <- pilot$p_sigmahat
		}
		
		out <- .Fortran("localrisk2", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
		dim(out) <- c(n, nh)
	}
    if(!is.null(dim(pilot$x))){
        xM = pilot$x[pilot$uniquerows, ]
    }else{
        xM = pilot$x[pilot$uniquerows]
    }
	
	return(list(risk=out, xM=xM, lambda.grid=pilot$lambda.grid))
}


newimprove_local <- function(pilot, newx, rval, sparameter=NULL, maxiter=10, theta=.15, boundary=c("symmetric", "none"), termval=0.1^4){
	lambda.grid <- pilot$lambda.grid
    nh <- length(lambda.grid)
    
    IAs <- fullA(pilot)
    
    ndim <- dim(pilot$x)
    if(!is.null(ndim)){
        xx <- pilot$x[pilot$uniquerows,]
        yy <- c(pilot$y)[pilot$uniquerows]
        nn <- nrow(newx)
    }else{
        xx <- pilot$x[pilot$uniquerows]
        yy <- pilot$y[pilot$uniquerows]
        nn <- length(newx)
    }
    
	n <- length(yy)
    #pilot esitmator
    smovlam <- matrix(NA, n, maxiter+1)
	predy <- matrix(NA, n, maxiter+1)
	predsigma <- rep(NA, maxiter+1)
	diagA <- matrix(NA, n, maxiter+1)
	eff.df <- rep(NA, maxiter+1)
	
	movrisk <- matrix(NA, nn, maxiter+1)
	movlambda <- matrix(NA, nn, maxiter+1)
	
	# pilot estimator
	smovlam[,1] <- pilot$p_lambda # rep(pilot$p_lambda, n)
	predy[,1] <- c(pilot$p_fitted)[pilot$uniquerows]
	predsigma[1] <- pilot$p_sigmahat
	eff.df[1] <- pilot$eff.df
	diagA[,1] <- 1/(1+pilot$p_lambda*pilot$matrice$D)
	
	# iterative moving risk and lambda
	movrisk[,1] <- NA
	movlambda[,1] <- NA # rep(pilot$p_lambda, nn)
	
	error <- Inf
	old_predlambda <- rep(pilot$p_lambda, n)
	
	iter <- 0
	############
	aa <- list()
	############
	
	
	while(error > termval){
	   iter <- iter + 1
		if(iter > maxiter) {
			warning(paste("The number of iteration is greater than ", maxiter, ". \n", sep=""))
			return(out)
		}
		#	cat(iter, "\n")
		
		
		elrisk <- newlocalrisk(pilot, ias=IAs$IAall, tra2s=IAs$Tra2all, rmethod="estimated")

		if(iter==1) risk0 <- elrisk

		elambdafit <- movinglambda(elrisk, newx, rval, boundary[1])
		
		#aa[[iter]] <- movingriskR(elrisk, newx, rval)
		moverisk <- elambdafit$risk
		elambda <- elambdafit$lambda
		##predrisk <- cbind(predrisk, moverisk)
		movrisk[, iter+1] <- moverisk
		##movlambda <- cbind(movlambda, elambda)
		movlambda[, iter+1] <- elambda
		
        if(!is.null(ndim)){
            if(is.null(sparameter)) stop("The sparameter is NULL \n")
            fite <- fastTps(newx[!is.na(elambda),], -log(elambda[!is.na(elambda)]), lambda=sparameter, theta=theta)
        }else{
            
            tmpx <- newx[!is.na(elambda)]
            #tmpy <- -log(elambda[!is.na(elambda)])
            tmpy <- log(elambda[!is.na(elambda)])
           
            fite <- Tps(tmpx, tmpy, lambda=sparameter)
        #    cat(fite$lambda, "\n")
        #    fite <- Tps(newx, -log(elambda))
        #    fite <- sreg(newx, -log(elambda), lambda=sparameter)
		    }
        #predlambda <- c(exp(-predict(fite, xx)))
        predlambda <- c(exp(predict(fite, xx)))
		
		locfit <- .Fortran("localfit", as.double(predlambda), as.double(yy), as.integer(n), 
						   as.double(pilot$matrices$D), as.double(pilot$QV),
						   fitted.values = double(n), diagA = double(n)) 
		
		diag_A <- locfit$diagA
		f_pilot <- c(locfit$fitted.values)
		sigma_pilot <- sqrt(sum((yy - f_pilot)^2)/(n - sum(diag_A)))
        #sigma_pilot <- 0.0463065
		pilot$p_fitted <- f_pilot
		pilot$p_sigmahat <- sigma_pilot
		
######## History information
		##smovlam <- cbind(smovlam, predlambda)
		smovlam[, iter+1] <- predlambda
		##predy <- cbind(predy, f_pilot)
		predy[, iter+1] <- f_pilot
		#predsigma <- c(predsigma, sigma_pilot)
		predsigma[iter+1] <- sigma_pilot
		##diagA <- cbind(diagA, diag_A)
		diagA[,iter+1] <- diag_A
		#eff.df <- c(eff.df, sum(diag_A))
		eff.df[iter+1] <- sum(diag_A)
######## 		
		error <- sum((log(old_predlambda) - log(predlambda))^2)/n
		old_predlambda <- predlambda
		
		out <- list(smovlambda=smovlam[,1:(iter+1)], predy=predy[,1:(iter+1)], predsigma = predsigma[1:(iter+1)], 
		            movrisk=movrisk[,1:(iter+1)], movlambda = movlambda[,1:(iter+1)], 
		            finlrisk=elrisk, inilrisk=risk0, diagA=diagA[,1:(iter+1)], eff.df=eff.df[1:(iter+1)], iteration =iter, boundary=boundary[1])
	}
	#bb <<- aa
	return(out)
}


newimprove_local2 <- function(pilot, newx, rval, sparameter=NULL, maxiter=10, theta=.15, boundary=c("symmetric", "none"), termval=0.1^4, IAs=NULL){
	lambda.grid <- pilot$lambda.grid
    nh <- length(lambda.grid)
    
    if ( is.null(IAs) ) {
    	IAs <- fullA(pilot)
	}    
    ndim <- dim(pilot$x)
    if(!is.null(ndim)){
        xx <- pilot$x[pilot$uniquerows,]
        yy <- c(pilot$y)[pilot$uniquerows]
        nn <- nrow(newx)
    }else{
        xx <- pilot$x[pilot$uniquerows]
        yy <- pilot$y[pilot$uniquerows]
        nn <- length(newx)
    }
    
	n <- length(yy)
    #pilot esitmator
    smovlam <- matrix(NA, n, maxiter+1)
	predy <- matrix(NA, n, maxiter+1)
	predsigma <- rep(NA, maxiter+1)
	diagA <- matrix(NA, n, maxiter+1)
	eff.df <- rep(NA, maxiter+1)
	
	movrisk <- matrix(NA, nn, maxiter+1)
	movlambda <- matrix(NA, nn, maxiter+1)
	
	# pilot estimator
	smovlam[,1] <- pilot$p_lambda # rep(pilot$p_lambda, n)
	predy[,1] <- c(pilot$p_fitted)[pilot$uniquerows]
	predsigma[1] <- pilot$p_sigmahat
	eff.df[1] <- pilot$eff.df
	diagA[,1] <- 1/(1+pilot$p_lambda*pilot$matrice$D)
	
	# iterative moving risk and lambda
	movrisk[,1] <- NA
	movlambda[,1] <- NA # rep(pilot$p_lambda, nn)
	
	error <- Inf
	old_predlambda <- rep(pilot$p_lambda, n)
	
	iter <- 0
	############
	aa <- list()
	############
	
	
	while(error > termval){
	   iter <- iter + 1
		if(iter > maxiter) {
			warning(paste("The number of iteration is greater than ", maxiter, ". \n", sep=""))
			return(out)
		}
		#	cat(iter, "\n")
		
		
		elrisk <- newlocalrisk(pilot, ias=IAs$IAall, tra2s=IAs$Tra2all, rmethod="estimated")

		if(iter==1) risk0 <- elrisk

		elambdafit <- movinglambda(elrisk, newx, rval, boundary[1])
		
		#aa[[iter]] <- movingriskR(elrisk, newx, rval)
		moverisk <- elambdafit$risk
		elambda <- elambdafit$lambda
		##predrisk <- cbind(predrisk, moverisk)
		movrisk[, iter+1] <- moverisk
		##movlambda <- cbind(movlambda, elambda)
		movlambda[, iter+1] <- elambda
		
        if(!is.null(ndim)){
            if(is.null(sparameter)) stop("The sparameter is NULL \n")
            fite <- fastTps(newx[!is.na(elambda),], -log(elambda[!is.na(elambda)]), lambda=sparameter, theta=theta)
        }else{
            
            tmpx <- newx[!is.na(elambda)]
            #tmpy <- -log(elambda[!is.na(elambda)])
            tmpy <- log(elambda[!is.na(elambda)])
           
            fite <- Tps(tmpx, tmpy, lambda=sparameter)
        #    cat(fite$lambda, "\n")
        #    fite <- Tps(newx, -log(elambda))
        #    fite <- sreg(newx, -log(elambda), lambda=sparameter)
		    }
        #predlambda <- c(exp(-predict(fite, xx)))
        predlambda <- c(exp(predict(fite, xx)))
		
		locfit <- .Fortran("localfit", as.double(predlambda), as.double(yy), as.integer(n), 
						   as.double(pilot$matrices$D), as.double(pilot$QV),
						   fitted.values = double(n), diagA = double(n)) 
		
		diag_A <- locfit$diagA
		f_pilot <- c(locfit$fitted.values)
		sigma_pilot <- sqrt(sum((yy - f_pilot)^2)/(n - sum(diag_A)))
        #sigma_pilot <- 0.0463065
		pilot$p_fitted <- f_pilot
		pilot$p_sigmahat <- sigma_pilot
		
######## History information
		##smovlam <- cbind(smovlam, predlambda)
		smovlam[, iter+1] <- predlambda
		##predy <- cbind(predy, f_pilot)
		predy[, iter+1] <- f_pilot
		#predsigma <- c(predsigma, sigma_pilot)
		predsigma[iter+1] <- sigma_pilot
		##diagA <- cbind(diagA, diag_A)
		diagA[,iter+1] <- diag_A
		#eff.df <- c(eff.df, sum(diag_A))
		eff.df[iter+1] <- sum(diag_A)
######## 		
		error <- sum((log(old_predlambda) - log(predlambda))^2)/n
		old_predlambda <- predlambda
		
		out <- list(smovlambda=smovlam[,1:(iter+1)], predy=predy[,1:(iter+1)], predsigma = predsigma[1:(iter+1)], 
		            movrisk=movrisk[,1:(iter+1)], movlambda = movlambda[,1:(iter+1)], 
		            finlrisk=elrisk, inilrisk=risk0, diagA=diagA[,1:(iter+1)], eff.df=eff.df[1:(iter+1)], iteration =iter, boundary=boundary[1])
	}
	#bb <<- aa
	return(out)
}



devlambda <- function(risk, xgrid = NULL, nxgrid=50, pxgrid=NULL, npxgrid=500, ntrunc=10, sparameter=0.1^8){
	x <- risk$xM
	lambda.grid <- risk$lambda.grid
	risk <- risk$risk
	
	
	if(is.null(xgrid)){
		xgrid <- seq(min(x), max(x),, nxgrid)
		xgrid <- xgrid[which(xgrid > x[order(x)[ntrunc]])]
		xgrid <- xgrid[which(xgrid < x[order(x, decreasing=TRUE)[ntrunc]])]
		xgrid <- seq(min(xgrid), max(xgrid),, nxgrid)
	}
	
	leftlambda <- rep(NA, length(xgrid))
	rightlambda <- rep(NA, length(xgrid))
	
	for(i in 1:length(xgrid)){
		left_risk <- risk[which(x < xgrid[i]) ,]
		right_risk <- risk[which(x >= xgrid[i]) ,]
		leftlambda[i] <- lambda.grid[which.min(apply(left_risk, 2, mean))]
		rightlambda[i] <- lambda.grid[which.min(apply(right_risk, 2, mean))]
	}
	
	dlambda <- abs(log(leftlambda)-log(rightlambda))
	fitdiff <- sreg(xgrid, dlambda, lambda=sparameter)
	if(is.null(pxgrid)){
		pxgrid <- seq(min(xgrid), max(xgrid),, npxgrid)
	}
	dplambda <- predict(fitdiff, pxgrid)
	devpred <- predict(fitdiff, pxgrid,  derivative=1)
	out <- list(xgrid = xgrid, dlam=dlambda, pxgrid=pxgrid, dplam=dplambda, devdlam = devpred,leftlambda=leftlambda, rightlambda=rightlambda )
	return(out)
}

newimprove_local3 <- function(pilot, newx, rval, sparameter=NULL, maxiter=10, theta=.15, boundary=c("symmetric", "none"), termval=0.1^4, IAs=NULL, xo=0.24){
	lambda.grid <- pilot$lambda.grid
    nh <- length(lambda.grid)
    
    if ( is.null(IAs) ) {
    	IAs <- fullA(pilot)
	}  
    
    ndim <- dim(pilot$x)
    if(!is.null(ndim)){
        xx <- pilot$x[pilot$uniquerows,]
        yy <- c(pilot$y)[pilot$uniquerows]
        nn <- nrow(newx)
    }else{
        xx <- pilot$x[pilot$uniquerows]
        yy <- pilot$y[pilot$uniquerows]
        nn <- length(newx)
    }
    
	n <- length(yy)
    #pilot esitmator
    smovlam <- matrix(NA, n, maxiter+1)
	predy <- matrix(NA, n, maxiter+1)
	predsigma <- rep(NA, maxiter+1)
	diagA <- matrix(NA, n, maxiter+1)
	eff.df <- rep(NA, maxiter+1)
	
	movrisk <- matrix(NA, nn, maxiter+1)
	movlambda <- matrix(NA, nn, maxiter+1)
	
	# pilot estimator
	smovlam[,1] <- pilot$p_lambda # rep(pilot$p_lambda, n)
	predy[,1] <- c(pilot$p_fitted)[pilot$uniquerows]
	predsigma[1] <- pilot$p_sigmahat
	eff.df[1] <- pilot$eff.df
	diagA[,1] <- 1/(1+pilot$p_lambda*pilot$matrice$D)
	
	# iterative moving risk and lambda
	movrisk[,1] <- NA
	movlambda[,1] <- NA # rep(pilot$p_lambda, nn)
	
	error <- Inf
	old_predlambda <- rep(pilot$p_lambda, n)
	
	iter <- 0
	############
	aa <- list()
	############
	
	riskmats <<- NULL
	biasmats <<- NULL
	varimats <<- NULL
	
	indxo <- which(abs(xx - xo) < rval)
	
	while(error > termval){
	   iter <- iter + 1
		if(iter > maxiter) {
			warning(paste("The number of iteration is greater than ", maxiter, ". \n", sep=""))
			return(out)
		}
		#	cat(iter, "\n")
		
		
		elrisk <- newlocalrisk3(pilot, ias=IAs$IAall, tra2s=IAs$Tra2all, rmethod="estimated")
		
		riskmats <<- cbind(riskmats, apply(elrisk$risk[indxo,], 2, mean))
		biasmats <<- cbind(biasmats, apply(elrisk$bias[indxo,], 2, sum))
		varimats <<- cbind(varimats, apply(elrisk$vari[indxo,], 2, sum))
		
		
		if(iter==1) risk0 <- elrisk
		
		elambdafit <- movinglambda(elrisk, newx, rval, boundary[1])
		
		#aa[[iter]] <- movingriskR(elrisk, newx, rval)
		moverisk <- elambdafit$risk
		elambda <- elambdafit$lambda
		##predrisk <- cbind(predrisk, moverisk)
		movrisk[, iter+1] <- moverisk
		##movlambda <- cbind(movlambda, elambda)
		movlambda[, iter+1] <- elambda
		
        if(!is.null(ndim)){
            if(is.null(sparameter)) stop("The sparameter is NULL \n")
            fite <- fastTps(newx[!is.na(elambda),], -log(elambda[!is.na(elambda)]), lambda=sparameter, theta=theta)
        }else{
            
            tmpx <- newx[!is.na(elambda)]
            #tmpy <- -log(elambda[!is.na(elambda)])
            tmpy <- log(elambda[!is.na(elambda)])
           
            fite <- Tps(tmpx, tmpy, lambda=sparameter)
        #    cat(fite$lambda, "\n")
        #    fite <- Tps(newx, -log(elambda))
        #    fite <- sreg(newx, -log(elambda), lambda=sparameter)
		    }
        #predlambda <- c(exp(-predict(fite, xx)))
        predlambda <- c(exp(predict(fite, xx)))
		
		locfit <- .Fortran("localfit", as.double(predlambda), as.double(yy), as.integer(n), 
						   as.double(pilot$matrices$D), as.double(pilot$QV),
						   fitted.values = double(n), diagA = double(n)) 
		
		diag_A <- locfit$diagA
		f_pilot <- c(locfit$fitted.values)
		sigma_pilot <- sqrt(sum((yy - f_pilot)^2)/(n - sum(diag_A)))
        #sigma_pilot <- 0.0463065
		pilot$p_fitted <- f_pilot
		pilot$p_sigmahat <- sigma_pilot
		
######## History information
		##smovlam <- cbind(smovlam, predlambda)
		smovlam[, iter+1] <- predlambda
		##predy <- cbind(predy, f_pilot)
		predy[, iter+1] <- f_pilot
		#predsigma <- c(predsigma, sigma_pilot)
		predsigma[iter+1] <- sigma_pilot
		##diagA <- cbind(diagA, diag_A)
		diagA[,iter+1] <- diag_A
		#eff.df <- c(eff.df, sum(diag_A))
		eff.df[iter+1] <- sum(diag_A)
######## 		
		error <- sum((log(old_predlambda) - log(predlambda))^2)/n
		old_predlambda <- predlambda
		
		out <- list(smovlambda=smovlam[,1:(iter+1)], predy=predy[,1:(iter+1)], predsigma = predsigma[1:(iter+1)], 
		            movrisk=movrisk[,1:(iter+1)], movlambda = movlambda[,1:(iter+1)], 
		            finlrisk=elrisk, inilrisk=risk0, diagA=diagA[,1:(iter+1)], eff.df=eff.df[1:(iter+1)], iteration =iter, boundary=boundary[1])
	}
	#bb <<- aa
	return(out)
}


newlocalrisk3 <- function(pilot, ias=NULL, tra2s=NULL, rmethod=c("estimated", "true", "corrected"), pilotf=NULL, pilotsigma=NULL){
	
	yM <- pilot$y[pilot$uniquerows]
	n <- length(yM)
	
	
	nh <- length(pilot$lambda.grid)
	
	if(rmethod[1]=="estimated"){
		
		if(is.null(pilotf)){
			if(length(pilot$p_fitted) != n){
				pilotf <- c(pilot$p_fitted)[pilot$uniquerows]
			}else{
				pilotf <- c(pilot$p_fitted)
			}
		}
		if(is.null(pilotsigma)){
			pilotsigma <- pilot$p_sigmahat
		}
		
		if ( is.null(ias) | is.null(tra2s) ){
            out <- .Fortran("localrisk1", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        } else {
            out <- .Fortran("localriskall2", as.integer(n), as.integer(nh),
						as.double(pilotf), as.double(pilotsigma), as.double(ias), as.double(tra2s),
						risk=double(n*nh), bias=double(n*nh), vari=double(n*nh))[c("risk", "bias", "vari")]
            dim(out$risk) <- c(n, nh)
            dim(out$bias) <- c(n, nh)
            dim(out$vari) <- c(n, nh)
        
        }
	}
	
	if(rmethod[1]=="true"){
		if(is.null(pilotf) | is.null(pilotsigma)){
			if(is.null(pilot$tf)) {
				warning("True f and variance must be...")
				return("Error")	
			}else{
				if(length(pilot$tf) != n){
					pilotf <- pilot$tf[pilot$uniquerows]
				}else{
					pilotf <- pilot$tf
				}
				
				pilotsigma <- pilot$tsigma
			}
		}
		if ( is.null(ias) | is.null(tra2s) ){
            out <- .Fortran("localrisk1", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        } else {
            out <- .Fortran("localriskall", as.integer(n), as.integer(nh),
						as.double(pilotf), as.double(pilotsigma), as.double(ias), as.double(tra2s),
						risk=double(n*nh))$risk
            dim(out) <- c(n, nh)
        
        }        
	}
	
	if(rmethod[1]=="corrected"){
		if(is.null(pilotf)){
			pilotf <- yM
		}
		if(is.null(pilotsigma)){
			pilotsigma <- pilot$p_sigmahat
		}
		
		out <- .Fortran("localrisk2", as.double(pilot$lambda.grid), as.integer(nh),
						as.double(pilot$matrices$D), as.double(pilot$QV), as.integer(n),
						as.double(pilotf), as.double(pilotsigma), risk=double(n*nh))$risk
		dim(out) <- c(n, nh)
	}
    if(!is.null(dim(pilot$x))){
        xM = pilot$x[pilot$uniquerows, ]
    }else{
        xM = pilot$x[pilot$uniquerows]
    }
	
	return(list(risk=out$risk, xM=xM, lambda.grid=pilot$lambda.grid, bias=out$bias, vari=out$vari))
}