fast.s <- function(x, y, int, N, k=2, 
        best.r=5, b=.5, cc=1.56, 
        seed) {
    #
    # int = 1 -> add a column of ones to x
    # N = cant de sub-samples
    # k = number of refining iterations in ea. 
    #     subsample
    # k = 0 means "raw-subsampling"
    # b = right hand side of the equation
    # cc = corresponding tuning constant
    # best.r = number of "best betas" to remember
    #          from the subsamples. These will be later
    #          iterated until convergence


    # the objective function, we solve loss.S(u,s,cc)=b for "s"
    loss.S <- function(u,s,cc) mean(rho(u/s,cc) )


    re.s <- function(x,y,initial.beta,initial.scale,k,conv,b,cc) {
    # does "k" IRWLS refining steps from "initial.beta"
    #
    # if "initial.scale" is present, it's used, o/w the MAD is used
    # k = number of refining steps
    # conv = 0 means "do k steps and don't check for convergence"
    # conv = 1 means "stop when convergence is detected, or the
    #                 maximum number of iterations is achieved"
    # b and cc = tuning constants of the equation
    # 

    f.w <- function(u, cc) 
    {
        # weight function = psi(u)/u
        tmp <- (1 - (u/cc)^2)^2
        tmp[ abs(u/cc) > 1 ] <- 0
        return(tmp)
    }

    n <- dim(x)[1]
    p <- dim(x)[2]
    res <- y - x %*% initial.beta
    if( missing( initial.scale ) ) 
        initial.scale <- scale <- median(abs(res))/.6745
    else
        scale <- initial.scale
    
    if( conv == 1) k <- 50
    #
    # if conv == 1 then set the max no. of iterations to 50
    # magic number alert!!!
    
    beta <- initial.beta

    lower.bound <- median(abs(res))/cc

    
    for(i in 1:k) {
        # do one step of the iterations to solve for the scale
        scale.super.old <- scale
        #lower.bound <- median(abs(res))/1.56
        scale <- sqrt( scale^2 * mean( rho( res / scale, cc ) ) / b     )
        # now do one step of IRWLS with the "improved scale"
        weights <- f.w( res/scale, cc )
        W <- matrix(weights, n, p)
        xw <- x * sqrt(W)
        yw <- y *   sqrt(weights)
        beta.1 <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
        if(any(is.na(beta.1))) { beta.1 <- initial.beta
                                    scale <- initial.scale
                                     break
        }
        if( (conv==1) )
        {
            # check for convergence
            if( norm( beta - beta.1 ) / norm(beta) < 1e-20 ) break
            # magic number alert!!!
        }
        res <- y - x %*% beta.1
        beta <- beta.1
    }
    
    res <- y - x %*% beta
    # get the residuals from the last beta
    return(list(beta.rw = beta.1, scale.rw = scale))
    
    }


    scale1 <- function(u, b, cc, initial.sc=median(abs(u))/.6745) {
        # find the scale, full iterations
        max.it <- 200
        # magic number alert
        #sc <- median(abs(u))/.6745
        sc <- initial.sc
        i <- 0
        eps <- 1e-20
        # magic number alert
        err <- 1
        while( ( (i <- i+1) < max.it ) && (err > eps) ) {
            sc2 <- sqrt( sc^2 * mean( rho( u / sc, cc ) ) / b   )
            err <- abs(sc2/sc - 1)
            sc <- sc2
        }
        return(sc)
    }





    
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    if( int == 1) {
            x <- cbind(rep(1,n), x)
            p <- p + 1
    }
    
    # save existing random seed
    if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
    
    if(!missing(seed)) set.seed(seed)

    best.betas <- matrix(0, best.r, p)
    best.scales <- rep(1e20, best.r)
    s.worst <- 1e20
    n.ref <- 1
    
    for(i in 1:N) {
        #
        # get a subsample
        #
        singular <- T
        while(singular==T) {
        indices <- sample(n,p)
        xs <- x[indices,]
        ys <- y[indices]
        beta <- our.solve(xs,ys)
        singular <- any(is.na(beta))
        }
        if(k>0) { 
            # do the refining
            tmp <- re.s(x=x,y=y,initial.beta=beta,k=k,conv=0,b=b,cc=cc)
            beta.rw <- tmp$beta.rw
            scale.rw <- tmp$scale.rw
            res.rw <- y - x %*% beta.rw
        } else { #k = 0 means "no refining"
            beta.rw <- beta
            res.rw <- y - x %*% beta.rw
            scale.rw <- median(abs(res.rw))/.6745
        }
        if( i > 1 ) { 
            # if this isn't the first iteration....
            scale.test <- loss.S(res.rw,s.worst,cc)
            if( scale.test < b ) {
                s.best <- scale1(res.rw,b,cc,scale.rw)
                ind <- order(best.scales)[best.r]
                best.scales[ind] <- s.best
                best.betas[ind,] <- beta.rw
                s.worst <- max(best.scales)
            }
        } else { # if this is the first iteration, then this is the best beta...
                best.scales[best.r]  <- scale1(res.rw,b,cc,scale.rw)
                best.betas[best.r,] <- beta.rw
        }
    }
    # do the complete refining step until convergence (conv=1) starting
    # from the best subsampling candidate (possibly refined)
    super.best.scale <- 1e20
    # magic number alert
    for(i in best.r:1) {
        tmp <- re.s(x=x,y=y,initial.beta=best.betas[i,],
                initial.scale=best.scales[i],k=0,conv=1,b=b,cc=cc)
        if(tmp$scale.rw < super.best.scale) {
            super.best.scale <- tmp$scale.rw
            super.best.beta <- tmp$beta.rw
        }   
    }
    
    # restore seed existing before call
    if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
    
    return(list(as.vector(super.best.beta), super.best.scale))

}





our.solve <- function(a,b) {
    a <- qr(a)
    da <- dim(a$qr)
    if(a$rank < (p <- da[2]))
        return(NA)
    else qr.coef(a, b)
}



rho <- function(u, cc=1.56) {
    w <- abs(u)<=cc
    v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
    v <- v*6/cc^2
    return(v)
}


norm <- function(x) sqrt( sum( x^2 ) )
