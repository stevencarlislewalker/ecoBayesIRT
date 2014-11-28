##' Construct an index for a Bayesian ordination (BO) mcmc object
##'
##' @param mcmc An mcmc object
##' @return A data frame with one row per model parameter.  Each
##' row gives the (1) parameter type, (2) object associated with 
##' that parameter (e.g. species name), and (3) the ordination 
##' axis associated with that parameter (if it is an axis-specific
##' parameter).
##' @export
BOindex <- function(mcmc) {
    namesList <- strsplit(colnames(mcmc), ".", TRUE)
    namesDF <- as.data.frame(t(sapply(namesList, "[", 1:3)),
                             stringsAsFactors = FALSE)
    names(namesDF) <- c("type", "object", "axis")
    namesDF$axis <- as.numeric(namesDF$axis)
    oneAxis <- all(is.na(namesDF$axis))
    if(oneAxis) {
        namesDF <- within(namesDF, {
            axis[type != "alpha"] <- 1
        })
    }
    return(namesDF)
}

##' Get parameters from a Bayesian ordination (getBO)
##' 
##' @param mcmc An mcmc object
##' @param Y Data matrix that matches the mcmc object
##' @param index The results of \code{\link{BOindex}}
##' @param name What to extract 
##'        (site scores, x, intercepts, a, or species scores, b)
##' @return extracted object
##' @export
getBO <- function(mcmc, Y, index, name = c("x", "a", "b")) {

    if(missing(index)) index <- BOindex(mcmc)
    x <- a <- b <- NULL
    if("a" %in% name) {
        a <- -mcmc[, index$type == "alpha"]
    }
    nAxes <- max(index$axis, na.rm = TRUE)
    if("b" %in% name) {
        b <- list()
        for(i in 1:nAxes) {
            cond <- with(index,
                         (type == "beta") &
                         (axis == i))
            b[[i]] <- mcmc[, cond]
        }
        b <- aperm(do.call(abind, c(b, list(along = 0))), c(2, 3, 1))
        dimnames(b)[[2]] <- colnames(Y)
    }
    if("x" %in% name) {
        x <- list()
        for(i in 1:nAxes) {
            cond <- with(index,
                         (type == "theta") &
                         (axis == i))
            x[[i]] <- mcmc[, cond]
        }
        x <- aperm(do.call(abind, c(x, list(along = 0))), c(2, 3, 1))
        dimnames(x)[[2]] <- rownames(Y)
    }
    
    return(list(a = a, b = b, x = x))
}

##' Rotate a Bayesian ordination
##'
##' @param BO Results of \code{\link{getBO}}
##' @return A Bayesian ordination with rotated axes
##' @export
BOrotate <- function(BO) {
    within(BO, {
        xCenter <- sweep(x, c(1, 3), apply(x, c(1, 3), mean))
        svdXCenter <- apply(xCenter, 1, svd)
        d <- dim(b)[3]
        recip <- rep(1/n, n)
        for(i in 1:dim(a)[1]) {
            a[i,] <- a[i,] + b[i,,] %*% t(x[i,,]) %*% recip
            b[i,,] <- b[i,,] %*% svdXCenter[[i]]$v %*% diag(svdXCenter[[i]]$d, d, d)
            x[i,,] <- svdXCenter[[i]]$u
        }
    })
}

##' Switch axes labels
##'
##' Attempts to switch axis labels such that more highly correlated
##' MCMC samples are labeled as the same axis.
##'
##' @param BO Results of \code{\link{getBO}}
##' @return A Bayesian ordination with switched axes
##' @export
BOswitch <- function(BO) {
    within(BO, {
        for(i in 1:nrow(a)) {
            within <- c(x[-i,,1] %*% x[i,,1], x[-i,,2] %*% x[i,,2])
            among <- c(x[,,2] %*% x[i,,1], x[,,1] %*% x[i,,2])
            if(mean(abs(among)) > mean(abs(within))) {
                x[i,,1:2] <- x[i,,2:1]
                b[i,,1:2] <- b[i,,2:1]
            }
        }
    })
}


##' Flip axes in a Bayesian ordination
##' 
##' @param BO Results of \code{\link{getBO}}
##' @param refSites Vector of reference sites for each axis
##' @return A Bayesian ordination with flipped axes
##' @export
BOflip <- function(BO, refSites) {
    .fn <- function(ii) ifelse(BO$x[, refSites[ii], ii] < 0, -1, 1)
    idntMult <- sapply(1:length(refSites), .fn)
    BO$x <- sweep(BO$x, c(1, 3), idntMult, "*")
    BO$b <- sweep(BO$b, c(1, 3), idntMult, "*")
    return(BO)
}

##' Compute the posterior of the linear predictor in a Bayesian
##' ordination
##'
##' @param BO Results of \code{\link{getBO}}
##' @return And array containing the posterior of the linear predictor
##' with dimensions relating to (1) MCMC iterations, (2) sites, and
##' (3) species
##' @export
BOeta <- function(BO) {
    with(BO, {
        nSite <- dim(BO$x)[2]
        nSpec <- dim(BO$b)[2]
        nSamp <- dim(BO$x)[1]
        outerSamp <- function(i) outer(x[i,,], b[i,,])
        A <- aperm(array(a, c(nSamp, nSpec, nSite)), c(1, 3, 2))
        XB <- aperm(simplify2array(lapply(seq_len(nSamp), outerSamp)), c(3,1,2))
        A + XB
    })
}
