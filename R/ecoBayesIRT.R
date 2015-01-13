## ecoBayesIRT: Helper functions for ecological Bayesian item response
## theory in R using MCMCpack
## 
## Copyright (C) 2014  Steven Carlisle Walker
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



##' Helper functions for ecological Bayesian item response theory in R
##' using MCMCpack
##'
##' Steps: (1) Fit a model using functions from the \code{MCMCpack}
##' package such as \code{\link{MCMCirt1d}} or
##' \code{\link{MCMCirtKd}}.  (2) Do diagnostics and convergence
##' checks on the resulting MCMC objects.  (3) Separate the three
##' types of parameters (site scores, \code{x}, species intercepts,
##' \code{a}, and species scores, \code{b}), by passing the MCMC
##' through the \code{\link{getBO}} function.  (4) Use the other
##' functions of this package to analyze the results
##' (\code{\link{BOeta}}; \code{\link{BOflip}}; \code{\link{BOp}};
##' \code{\link{BOrotate}}; \code{\link{BOswitch}}).
##'
##' @docType package
##' @name ecoBayesIRT
##' @import MCMCpack abind coda logspline plyr
##'
##' @examples
##' ## simulate data
##' nSite <- 25
##' nSpec <- 10
##' nIter <- 10000
##' set.seed(1)
##' x <- rnorm(nSite)
##' a <- rnorm(nSpec)
##' b <- rnorm(nSpec)
##' eta <- cbind(1, x) %*% rbind(a, b)
##' e <- matrix(rnorm(nSite*nSpec), nSite, nSpec)
##' p <- pnorm(eta + e)
##' Y <- matrix(rbinom(nSite*nSpec, 1, p), nSite, nSpec)
##' rownames(Y) <- letters[1:nSite]
##' colnames(Y) <- LETTERS[1:nSpec]
##'
##' ## fit model
##' Yirt <- MCMCirt1d(Y, mcmc = nIter, store.item = TRUE)
##'
##' ## getBO
##' BO <- getBO(Yirt, Y)
##'
##' ## fix identifiability
##' (refSite <- findRefSite(BO))
##' BO <- BOflip(BO, refSite)
##'
##' ## get fitted values
##' pFit <- BOp(BO)
##' pFitMean <- apply(pFit, c(2, 3), mean)
##' boxplot(pFitMean ~ Y)
##'
##' ## ordination plot
##' ordInt <- t(apply(BO$x, c(2, 3), quantile, probs = c(0.025, 0.975))[,,1])
NULL


##' Extract and replace pieces of an item response theory mcmc
##'
##' \code{extractIRT} and \code{replaceIRT} are only to be used in
##' interactive mode (i.e. not used within functions).  see
##' \url{http://adv-r.had.co.nz/Computing-on-the-language.html#calling-from-another-function}.
##' It is much safer and easier to read to use the \code{x}, \code{b},
##' and \code{a} functions, which are used to extract and replace all
##' gradient values, species slopes, and species intercepts (note
##' these functions return as arrays).  \code{getAxisParamArray} and
##' \code{setAxisParamArray} are mostly for internal use, because
##' usually \code{x}, \code{b}, and \code{a} will provide everything
##' someone might need.
##'
##' @param mcmc An \code{\link{mcmc}} object from a fitted IRT model.
##' @param condition An expression evaluating to a
##' \code{\link{logical}} vector the same length as the number of
##' columns (parameters) in \code{mcmc}.
##' @param index Optional (not really needed)
##' @param value Replacement value.
##' @param param Which parameter to obtain/replace (\code{theta},
##' \code{beta}, or \code{alpha}).
##' @param ord \code{itersIRT(mcmc)} by \code{axesIRT(mcmc)} matrix
##' with the orderings of the axes.
##' @return Another \code{\link{mcmc}} object with (1) only the
##' extracted piece (for \code{extractIRT}) or (2) the replaced values
##' (for \code{replaceIRT}).
##' @export
extractIRT <- function(mcmc, condition, index) {
    if(missing(index)) index <- BOindex(mcmc)
    mcmc[, eval(substitute(condition), index)]
}
##' @rdname extractIRT
##' @export
replaceIRT <- function(mcmc, condition, value, index) {
    if(missing(index)) index <- BOindex(mcmc)
    mcmc[, eval(substitute(condition), index)] <- value
    return(mcmc)
}
##' @rdname extractIRT
##' @export
a <- function(mcmc, index) {
    if(missing(index)) index <- BOindex(mcmc)
    condition <- with(index, type == "alpha")
    out <- mcmc[, condition]
    colnames(out) <- index[condition, "object"]
    return(out)
}
##' @rdname extractIRT
##' @export
b <- function(mcmc, index) {
    if(missing(index)) index <- BOindex(mcmc)
    getAxisParamArray(mcmc, index, "beta")
}
##' @rdname extractIRT
##' @export
x <- function(mcmc, index) {
    if(missing(index)) index <- BOindex(mcmc)
    getAxisParamArray(mcmc, index, "theta")
}
##' @rdname extractIRT
##' @export
`a<-` <- function(mcmc, value) {
    index <- BOindex(mcmc)
    mcmc[, with(index, type == "alpha")] <- value
    return(mcmc)
}
##' @rdname extractIRT
##' @export
`b<-` <- function(mcmc, value) {
    mcmc <- setAxisParamArray(mcmc, param = "beta", value = value)
}
##' @rdname extractIRT
##' @export
`x<-` <- function(mcmc, value) {
    mcmc <- setAxisParamArray(mcmc, param = "theta", value = value)
}
##' @rdname extractIRT
##' @export
getAxisParamArray <- function(mcmc, index, param = c("beta", "theta")) {
    if(missing(index)) index <- BOindex(mcmc)
    param <- param[1]
    out <- list()
    for(i in 1:axesIRT(mcmc, index)) {
        conditioni <- with(index, (type == param) & (axis == i))
        out[[i]] <- mcmc[, conditioni]
    }
    out <- aperm(do.call(abind, c(out, list(along = 0))), c(2, 3, 1))
    dimnames(out)[[2]] <- subset(index, conditioni)$object
    return(out)
}
##' @rdname extractIRT
##' @export
setAxisParamArray <- function(mcmc, index, param = c("beta", "theta"), value) {
    if(missing(index)) index <- BOindex(mcmc)
    param <- param[1]
    for(i in 1:axesIRT(mcmc, index)) {
        conditioni <- with(index, (type == param) & (axis == i))
        mcmc[, conditioni] <- value[, , i]
    }
    return(mcmc)
}
##' @rdname extractIRT
##' @export
swapAxisParamArray <- function(mcmc, index, param = c("beta", "theta"), ord) {
    ## nothing here yet
}

##' The dimensions of an item response theory mcmc
##'
##' @inheritParams extractIRT
##' 
##' @export
dimIRT <- function(mcmc, index) {
    if(missing(index)) index <- BOindex(mcmc)
    c(iters = nrow(mcmc),
      sites = length(unique(subset(index, type == "theta")$object)),
      species = length(unique(subset(index, type == "beta")$object)),
      axes = max(index$axis, na.rm = TRUE))
}
##' @rdname dimIRT
##' @export
itersIRT   <- function(mcmc, index) dimIRT(mcmc, index)["iters"]
##' @rdname dimIRT
##' @export
sitesIRT   <- function(mcmc, index) dimIRT(mcmc, index)["sites"]
##' @rdname dimIRT
##' @export
speciesIRT <- function(mcmc, index) dimIRT(mcmc, index)["species"]
##' @rdname dimIRT
##' @export
axesIRT    <- function(mcmc, index) dimIRT(mcmc, index)["axes"]



##' Follows
##'
##' Syntactic sugar for composition of two functions. Useful for apply
##' functions.
##' 
##' @param f a function of one argument
##' @param g another function with range equal to the domain of
##' \code{f}
##' @param ... additional arguments to \code{g}
##' @return The function given by \code{f} follows \code{g} of
##' \code{x}.
##' @rdname follows
##' @export
##' @examples
##' sqr <- function(x) x^2
##' min10 <- function(x) x - 10
##' (sqr %f% min10)(2)
##' div5 <- function(x) x/5
##' (sqr %f% min10 %f% div5)(2)
##' (unique %f% unlist %f% `[`)(letters[1:10], c(3, 2))
follows <- function(f, g) function(x, ...) f(g(x, ...))

##' @rdname follows
##' @export
`%f%` <- function(f, g) follows(f, g)



##' Returns a vector or array or list of values obtained by applying a
##' function to margins of two arrays or matrices.
##'
##' @param X,Y arrays (each potentially a matrix)
##' @param XMARGIN,YMARGIN vectors giving the subscripts which the
##' function will be applied over.  see \code{\link{apply}} for more
##' information. for the special case here, note that
##' \code{all(dim(X)[XMARGIN] == dim(Y)[YMARGIN])} must be
##' \code{TRUE}.
##' @param FUN the function to be applied.
##' @param ... more arguments to \code{FUN}.
##' @return a new array, or list if \code{\link{simplify2array}}
##' cannot do it's job.
##' @examples
##' dims <- c(4, 2, 3, 1, 2)
##' set.seed(1)
##' (A <- array(rnorm(prod(dims[-2])), dims[-2]))
##' (B <- array(rnorm(prod(dims[-(3:4)])), dims[-(3:4)]))
##' apply2arrays(A, B, c(1, 4), c(1, 3), `%*%`)
##' @export
apply2arrays <- function(X, Y, XMARGIN, YMARGIN, FUN, ...) {
    dX <- dim(X)
    dY <- dim(Y)
                                        # marginal dimensions
    dXmarg <- dX[XMARGIN]
    dYmarg <- dY[YMARGIN]
    if(!all(dXmarg == dYmarg)) stop("margins do not match")

                                        # dimensions to marginalize
                                        # over
    XOVER <- setdiff(1:length(dX), XMARGIN)
    YOVER <- setdiff(1:length(dY), YMARGIN)
    dXover <- dX[XOVER]
    dYover <- dY[YOVER]

                                        # put marginal dimensions
                                        # first (FIXME: maybe not best
                                        # ordering?)
    X <- aperm(X, c(XMARGIN, XOVER))
    Y <- aperm(Y, c(YMARGIN, YOVER))

                                        # collapse marginal dimensions
                                        # into one
    dim(X) <- c(prod(dXmarg), dXover)
    dim(Y) <- c(prod(dYmarg), dYover)

                                        # convert FUN into a function
                                        # that takes a single list as
                                        # an argument
    sFUN <- splat(FUN)

                                        # loop over margin
    ans <- list()
    for(i in seq_len(dim(X)[1])) {
                                        # structure (instead of drop)
                                        # required to drop first (but
                                        # not any other) dimensions
        ans[[i]] <- sFUN(list(structure(take(X, 1, i), dim = dim(X)[-1]),
                              structure(take(Y, 1, i), dim = dim(Y)[-1])
                              ), ...)
    }

                                        # convert ans to an array, put
                                        # the margin first again, and
                                        # add back the dimensions of
                                        # the margin
    ans <- simplify2array(ans) # FIXME: maybe add SIMPLIFY to argument
                               # list?
    if(is.atomic(ans)) {
        nDimAns <- length(da <- dim(ans))
        ans <- aperm(ans, c(nDimAns, 1:(nDimAns-1)))
        dim(ans) <- c(dXmarg, da[-nDimAns])
    }
    return(ans)
}



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
                                        # the names of model
                                        # parameters
                                        # (i.e. colnames(mcmc)) are
                                        # stored with dots separating
                                        # three pieces of information:
                                        # (1) type of parameter
                                        # (i.e. alpha, beta, theta),
                                        # (2) object associated with
                                        # the parameter (e.g. a
                                        # species or a site name), and
                                        # (3) the axis that the
                                        # parameter is associated
                                        # with.  namesList is a list,
                                        # each element of which gives
                                        # these three pieces of
                                        # information for each model
                                        # parameter.
    namesList <- strsplit(colnames(mcmc), ".", TRUE)

                                        # put these three pieces of
                                        # information into a data
                                        # frame, with one row per
                                        # model parameter.
    namesDF <- as.data.frame(t(sapply(namesList, "[", 1:3)),
                             stringsAsFactors = FALSE)
    names(namesDF) <- c("type", "object", "axis")
    namesDF$axis <- as.numeric(namesDF$axis)

                                        # gracefully handle the
                                        # special case of one axis.
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
    checkBO(BO)
    BO <- within(BO, {
        xCenter <- sweep(x, c(1, 3), apply(x, c(1, 3), mean))
        svdXCenter <- apply(xCenter, 1, svd)
        d <- dim(b)[3]
        n <- dim(x)[2]
        recip <- rep(1/n, n)
        for(i in 1:dim(a)[1]) {
            a[i,] <- a[i,] + b[i,,] %*% t(x[i,,]) %*% recip
            b[i,,] <- b[i,,] %*% svdXCenter[[i]]$v %*% diag(svdXCenter[[i]]$d, d, d)
            x[i,,] <- svdXCenter[[i]]$u
        }
    })
    return(BO[c("a", "b", "x")])
}

##' Switch axes labels
##'
##' Attempts to switch axis labels such that more highly correlated
##' MCMC samples are labeled as the same axis.  Only for two-axis
##' models.
##'
##' @param BO Results of \code{\link{getBO}}
##' @return A Bayesian ordination with switched axes
##' @export
BOswitch <- function(BO) {
    checkBO(BO)
    BO <- within(BO, {
        for(i in 1:nrow(a)) {
            within <- c(x[-i,,1] %*% x[i,,1], x[-i,,2] %*% x[i,,2])
            among <- c(x[,,2] %*% x[i,,1], x[,,1] %*% x[i,,2])
            if(mean(abs(among)) > mean(abs(within))) {
                x[i,,1:2] <- x[i,,2:1]
                b[i,,1:2] <- b[i,,2:1]
            }
        }
    })
    return(BO[c("a", "b", "x")])
}


##' Flip axes in a Bayesian ordination
##' 
##' @param BO Results of \code{\link{getBO}}
##' @param refSites Vector of reference sites for each axis
##' @return A Bayesian ordination with flipped axes
##' @seealso \code{\link{mcmcFlip}}
##' @export
BOflip <- function(BO, refSites) {
    checkBO(BO)
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
##' @return An array containing the posterior of the linear predictor
##' with dimensions relating to (1) MCMC iterations, (2) sites, and
##' (3) species
##' @export
BOeta <- function(BO) {
    checkBO(BO)
    eta <- with(BO, {
        nSite <- dim(BO$x)[2]
        nSpec <- dim(BO$b)[2]
        nSamp <- dim(BO$x)[1]
        outerSamp <- function(i) x[i,,] %*% t(b[i,,]) #  outer(x[i,,], b[i,,])
        A <- aperm(array(a, c(nSamp, nSpec, nSite)), c(1, 3, 2))
        XB <- aperm(simplify2array(lapply(seq_len(nSamp), outerSamp)), c(3,1,2))
        A + XB
    })
    dimnames(eta)[2:3] <- list(dimnames(BO$x)[[2]],
                               dimnames(BO$b)[[2]])
    return(eta)
}

##' Compute the posterior of the fitted values from a Bayesian
##' ordination on the probability scale
##'
##' @param BO Results of \code{\link{getBO}} (not required if
##' \code{eta} is present)
##' @param eta Results of \code{\link{BOeta}} (not required if
##' \code{BO} is present)
##' @param .seed Seed for random number generator associated with the
##' standard normal residual error term
##' @return An array containing the posterior of the fitted values on
##' the probability scale with dimensions relating to (1) MCMC
##' iterations, (2) sites, and (3) species
##' @export
BOp <- function(BO, eta, .seed = 1) {
    if(missing(eta) & (!missing(BO))) {
        eta <- BOeta(BO)
    } else if(missing(eta)) {
        stop("must specify at least one of either BO or eta")
    }
    set.seed(.seed)
    pnorm(eta + array(rnorm(prod(dim(eta))), dim(eta)))
}


checkBO <- function(BO) {
    nms <- names(BO)
    if(length(setdiff(nms, c("x", "a", "b"))) > 0) {
        stop("need x, a, and b elements to be a valid BO")
    }
    nIter <- dim(BO$x)[1]
    nSite <- dim(BO$x)[2]
    nSpec <- dim(BO$b)[2]
    nAxes <- dim(BO$x)[3]
    if(nIter != dim(BO$b)[1]) stop("inconsistent number of iterations")
    if(nSpec != dim(BO$a)[2]) stop("inconsistent number of species")
    if(nAxes != dim(BO$b)[3]) stop("inconsistent number of axes")
}

##' Find reference site
##'
##' @param BO Results of \code{\link{getBO}}
##' @export
findRefSite <- function(BO) {
    densitiesX <- apply(BO$x, 2:3, logspline)
    denAtZero <- matrix(sapply(densitiesX, dlogspline, q = 0),
                        dim(BO$x)[2:3])
    apply(denAtZero, 2, which.min)
}

##' All possible pairwise comparisons among columns in binary or
##' probability matrices
##'
##' @param y binary or probability matrix
##' @return 4-d array with comparisons
##' @export
pairComp <- function(y){
    out <- array(dim = c(2, 2, ncol(y), ncol(y)))
    out[1,1,,] <- t(y) %*% y
    out[1,2,,] <- t(y) %*% (1-y)
    out[2,1,,] <- t(1-y) %*% y
    out[2,2,,] <- t(1-y) %*% (1-y)
    dimnames(out)[1:2] <- rep(list(c("present","absent")), 2)
    dimnames(out)[3:4] <- rep(list(colnames(y)), 2)
    return(out)
}


##' Compare binary data to binned probabilistic expectations
##'
##' @param y vector of binary data
##' @param p vector of probabilities
##' @param length number of bins
##' @return \code{data.frame} with observed proportions and expected
##' proportions in evenly spaced bins.
##' @export
fitProbBins <- function(y, p, length = 11) {
    obsExp <- data.frame(y, p)
    incr <- 1/(length-1)
    mid1 <- incr/2
    mid2 <- 1 - mid1
    obsProp <- tapply(obsExp$y, cut(obsExp$p, seq(0, 1, length = length)), mean)
    expProp <- seq(mid1, mid2, by = incr)
    data.frame(obsProp = obsProp,
               expProp = expProp)
}

##' Differentiate between specific forms of uni- and bi-modality
##'
##' If there is one mode on either side of zero, then this function
##' attempts to return the location of those two modes.  If there is
##' only one mode, the function attempts to return the location of
##' that mode.  The algorithm is based on optimizing over the negative
##' numbers, then over the positive numbers, and seeing if either is
##' on the boundary.  This is not well tested.
##'
##' @param den output of \code{\link{logspline}}
##' @param tol tolerance when comparing modes with zero
##' @return A vector with the location of one or two modes
##' @export
##' @examples
##' set.seed(1)
##' x <- rnorm(1000, c(-2, 2), 1)
##' den <- logspline(x)
##' plot(den)
##' abline(v = findModes(den))
findModes <- function(den, tol = 1e-2) {
                                        # find region with 'most' of
                                        # the probability
    importantRegion <- qlogspline(c(0.01, 0.99), den)
                                        # if the important region is
                                        # either all positive or all
                                        # negative, just optimize and
                                        # return a single mode
    if(all(importantRegion > 0) || all(importantRegion < 0)) {
        return(unname(optimize(dlogspline, importantRegion,
                               maximum = TRUE, fit = den)$maximum))
    }
                                        # compute two halves of the
                                        # real line
    interval <- list(c(den$range[1], 0),
                     c(0, den$range[2]))
                                        # optimize the density over
                                        # both halves
    opts <- lapply(interval, optimize,
                   f = dlogspline, maximum = TRUE,
                   fit = den)
                                        # extract the locations of the
                                        # maxima
    maxs <- unlist(lapply(opts, "[", "maximum"))
                                        # discard modes at zero
                                        # (because they are probably
                                        # not local maxima over the
                                        # entire real line)
    keeps <- !sapply(mapply(all.equal, maxs, c(0, 0),
                            SIMPLIFY = FALSE,
                            MoreArgs = list(tolerance = tol)),
                     isTRUE)
    return(unname(maxs[keeps]))
}



##' Flip axes in an mcmc object for a Bayesian ordination
##' 
##' @param mcmc An \code{\link{mcmc}} object
##' @param refSites Vector of reference sites for each axis
##' @return An \code{\link{mcmc}} object for a Bayesian ordination
##' with flipped axes
##' @seealso \code{\link{BOflip}}
##' @export
mcmcFlip <- function(mcmc, refSites) {
    index <- BOindex(mcmc)
    nAxes <- max(index$axis, na.rm = TRUE)
    for(i in 1:nAxes) {
        thetaInds <- with(index, (type == "theta") & (axis == i))
         betaInds <- with(index, (type ==  "beta") & (axis == i))
        idntMult <- ifelse(mcmc[, which(thetaInds)[refSites[i]]] < 0, -1, 1)
        mcmc[, thetaInds] <- sweep(mcmc[, thetaInds], 1, idntMult, "*")
        mcmc[,  betaInds] <- sweep(mcmc[,  betaInds], 1, idntMult, "*")
    }
    return(mcmc)
}

## mcmcFlipOrd <- function(mcmc, multOrd) {
##     index <- BOindex(mcmc)
##     nAxes <- max(index$axis, na.rm = TRUE)
##     nIter <- nrow(mcmc)
##     nSite <- with(index, type == "theta")
##     for(i in 1:nAxes) {
##         thetaInds <- with(index, (type == "theta") & (axis == i))
##          betaInds <- with(index, (type ==  "beta") & (axis == i))        
##         mcmc[, thetaInds] <- sweep(mcmc[, thetaInds], 1, multOrd[i,,2], "*")
##         mcmc[,  betaInds] <- sweep(mcmc[,  betaInds], 1, multOrd[i,,2], "*")
##     }
##     for(i in 1:nIter) {
##         thetaInds <- with(index, (type == "theta")
##     }
##     return(mcmc)
## }

## getOrdMult <- function(BO, eta, etaMean) {
##     nAxes <- dim(BO$x)[3]
##                                         # rotate (FIXME:  svd twice!)
##     if(nAxes > 1) BO <- BOrotate(BO)
    
##                                         # fill in anything that's
##                                         # missing
##     print("calculating eta stuff")
##     if(missing(eta)) eta <- BOeta(BO)
##     if(missing(etaMean)) etaMean <- apply(eta, 2:3, mean)

##                                         # calculate SVDs and
##                                         # correlations between SVDs of
##                                         # the posterior and the
##                                         # posterior mean
##     print("calculating means and cors")
##     .ufn <- function(xx) svd(scale(xx, scale = FALSE))$u[, 1:nAxes, drop = FALSE]
##     etaSvd <- alply(eta, 1, .ufn)
##     etaSvdMean <- .ufn(etaMean)
##     etaSvdCors <- aaply(BO$x, 1, cor, y = etaSvdMean)

##                                         # find the best ordering of
##                                         # the gradients in the
##                                         # posterior
##     print("calculating orderings")
##     .multfn <- function(xx) {
##         if(nAxes == 1L) dim(xx) <- c(length(xx), 1)
##         ord <- mult <- numeric(nrow(xx))
##         for(i in 1:nrow(xx)) {
##             ord[i] <- which.max(abs(xx[i, ]))
##             mult[i] <- sign(xx[i, ord[i]])
##             xx[, ord[i]] <- 0L
##         }
##         cbind(ord = ord, mult = mult)
##     }
##     if(nAxes == 1) {
##         ordMult <- laply(etaSvdCors, .multfn)
##     } else {
##         ordMult <- aaply(etaSvdCors, 1, .multfn)
##     }
##     return(ordMult)
## }

##' Fix 2d-ordinations
##'
##' @param BO 2d Bayesian ordination
##' @param eta optional results of \code{BOeta(BO)}
##' @param etaMean optional posterior mean of \code{eta}
##' @export
BO2dFix <- function(BO, eta, etaMean) {
    nAxes <- dim(BO$x)[3]
                                        # rotate (FIXME:  svd twice!)
    if(nAxes > 1) BO <- BOrotate(BO)
    
                                        # fill in anything that's
                                        # missing
    print("calculating eta stuff")
    if(missing(eta)) eta <- BOeta(BO)
    if(missing(etaMean)) etaMean <- apply(eta, 2:3, mean)

                                        # calculate SVDs and
                                        # correlations between SVDs of
                                        # the posterior and the
                                        # posterior mean
    print("calculating means and cors")
    .ufn <- function(xx) svd(scale(xx, scale = FALSE))$u[, 1:nAxes, drop = FALSE]
    etaSvd <- alply(eta, 1, .ufn)
    etaSvdMean <- .ufn(etaMean)
    etaSvdCors <- aaply(BO$x, 1, cor, y = etaSvdMean)

                                        # find the best ordering of
                                        # the gradients in the
                                        # posterior
    print("calculating orderings")
    .multfn <- function(xx) {
        if(nAxes == 1L) dim(xx) <- c(length(xx), 1)
        ord <- mult <- numeric(nrow(xx))
        for(i in 1:nrow(xx)) {
            ord[i] <- which.max(abs(xx[i, ]))
            mult[i] <- sign(xx[i, ord[i]])
            xx[, ord[i]] <- 0L
        }
        data.frame(ord = ord, mult = mult)
    }
    if(nAxes == 1) {
        ordMult <- llply(etaSvdCors, .multfn)
    } else {
        ordMult <- alply(etaSvdCors, 1, .multfn)
    }
    
    .fixfn <- function(i, BOelement, ordMult) {
        with(ordMult[[i]], {
            flip <- sweep(BOelement[i,,,drop = FALSE], 3, mult, "*")
            return(flip[, , ord])
        })
    }
    BO$x <- laply(seq_along(ordMult), .fixfn, BOelement = BO$x, ordMult = ordMult)
    BO$b <- laply(seq_along(ordMult), .fixfn, BOelement = BO$b, ordMult = ordMult)
    return(BO)
}

##' Two-dimensional quantile plot
##'
##' @param x numeric vector
##' @param y numeric vector the same length as \code{y}
##' @param p probability
##' @param ... other arguments to \code{\link{polygon}}
##' @return plots a convex hull that contains a fraction \code{p} of
##' the \code{x}-\code{y} points.
##' @export
quantile2d <- function(x, y, p, ...) {
    hull <- chull(x, y)
    xy <- cbind(x, y)
    mahal <- mahalanobis(xy, apply(xy, 2, mean), cov(xy))
    xyOrd <- xy[order(mahal), ][1:round(length(x)*p), ]
    hullOrd <- chull(xyOrd[, 1], xyOrd[, 2])
    polygon(xyOrd[hullOrd, 1], xyOrd[hullOrd, 2], ...)
}

##' Plot correlations between 2d indirect and 1d direct gradients
##'
##' @param xIndirect posterior distribution of indirect gradient
##' values (MCMC iterations-by-sites-by-2 axes array)
##' @param xDirect direct gradient (vector with same length as second
##' dimension of \code{xIndirect})
##' @param .seed random seed for generation of the null correlations
##' @param plotArgs list of arguments for setting up the plotting region
##' @param nullPtArgs list of arguments for plotting correlations from
##' a null distribution
##' @param xPtArgs list of arguments for plotting correlations with direct gradient
##' @param nullContArgs list of arguments for plotting null distribution contour
##' @param xContArgs list of arguments for plotting contour
##' @export
plot2dCors <- function(xIndirect, xDirect, .seed = 1,
                       plotArgs = list(las = 1, 
                           xlim = c(-1, 1), ylim = c(-1, 1),
                           xlab = "", ylab = ""),
                       nullPtArgs = list(col = rgb(1, 0, 0, 0.1)),
                       xPtArgs = list(col = rgb(0, 0, 1, 0.1)),
                       nullContArgs = list(lwd = 3, cont = 95),
                       xContArgs = list(lwd = 3, cont = 95)) {

    message("computing correlations...")
    xCors <- apply(xIndirect, 1, cor, xDirect)
    .nullfn <- function(i, xI, xD) {
        apply(xI[i,,], 2, cor, xD[sample(dim(xI)[2])])
    }
    set.seed(.seed)
    nullCors <- sapply(1:dim(xIndirect)[1], .nullfn,
                       xI = xIndirect,
                       xD = xDirect)
    
    message("computing densities...")
    xDen <- kde(t(xCors))
    nullDen <- kde(t(nullCors))
    
    message("plotting...")
    plotArgs <- c(list(x = 0, y = 0, type = "n"), plotArgs)
    nullPtArgs = c(list(x = t(nullCors)), nullPtArgs)
    xPtArgs = c(list(x = t(xCors)), xPtArgs)
    nullContArgs = c(list(x = nullDen, add = TRUE), nullContArgs)
    xContArgs = c(list(x = xDen, add = TRUE), xContArgs)
    do.call(plot, plotArgs)
    abline(v = 0, h = 0, col = grey(0.5))
    do.call(points, xPtArgs)
    do.call(points, nullPtArgs)
    do.call(plot, xContArgs)
    do.call(plot, nullContArgs)
    invisible()
}
