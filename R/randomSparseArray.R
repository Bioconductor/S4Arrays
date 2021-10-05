### =========================================================================
### randomSparseArray() and poissonSparseArray()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### randomSparseArray()
###

### Returns an SVT_SparseArray object of type "double".
randomSparseArray <- function(dim=1L, density=0.05)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (!isSingleNumber(density) || density < 0 || density > 1)
        stop(wmsg("'density' must be a number >= 0 and <= 1"))

    ## Start with an empty sparse array.
    ans <- new_SVT_SparseArray(dim, type="integer")

    ## Add the nonzero values to it.
    ans_len <- length(ans)
    nzcount <- as.integer(ans_len * density)
    Lindex <- sample.int(ans_len, nzcount)
    nzvals <- signif(rnorm(nzcount), 2)
    if (nzcount <= .Machine$integer.max) {
        ## Using an Mindex seems to be slightly faster (4%-5%) than using an
        ## Lindex but we can only do this when the resulting Mindex matrix
        ## has < 2^31 rows.
        ans[Lindex2Mindex(Lindex, dim(ans))] <- nzvals
    } else {
        ans[Lindex] <- nzvals
    }

    ans
}

randomSparseMatrix <- function(nrow=1L, ncol=1L, density=0.05)
{
    if (!isSingleNumber(nrow) || !isSingleNumber(ncol))
        stop(wmsg("'nrow' and 'ncol' must be single integers"))
    randomSparseArray(c(nrow, ncol), density)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### poissonSparseArray()
###

### Replacement for rpois() when 'n' is big and 'lambda' is small.
### For example:
###     .sparse_rpois(3e9, 0.005)  # takes about 1 min. and uses < 1G of RAM
###     rpois(3e9, 0.005)          # takes about 55 sec. and uses 12G of RAM
.sparse_rpois <- function(n, lambda, chunksize=5e6L)
{
    if (n == 0L)
        return(list(integer(0), integer(0)))
    nzidx_env <- new.env(parent=emptyenv())
    nzvals_env <- new.env(parent=emptyenv())
    offset <- 0  # double to avoid integer overflow when n >= 2^31
    k <- 1L
    while (offset < n) {
        nn <- n - offset
        if (nn > chunksize)
            nn <- chunksize
        vals <- rpois(nn, lambda)
        nzidx <- which(vals != 0L)
        key <- sprintf("%04d", k)
        assign(key, offset + nzidx, envir=nzidx_env)
        assign(key, vals[nzidx], envir=nzvals_env)
        offset <- offset + nn
        k <- k + 1L
    }
    nzidx <- as.list(nzidx_env, all.names=TRUE, sorted=TRUE)
    nzidx <- unlist(nzidx, recursive=FALSE, use.names=FALSE)
    nzvals <- as.list(nzvals_env, all.names=TRUE, sorted=TRUE)
    nzvals <- unlist(nzvals, recursive=FALSE, use.names=FALSE)
    list(nzidx, nzvals)
}

### Returns an SVT_SparseArray object of type "integer".
poissonSparseArray <- function(dim=1L, lambda=0.05)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (!isSingleNumber(lambda) || lambda < 0 || lambda > 1)
        stop(wmsg("'lambda' must be a number >= 0 and <= 1"))

    ## Start with an empty sparse array.
    ans <- new_SVT_SparseArray(dim, type="integer")

    ## Add the nonzero values to it.
    ans_len <- length(ans)
    srp <- .sparse_rpois(ans_len, lambda)
    ans[srp[[1L]]] <- srp[[2L]]

    ans
}

poissonSparseMatrix <- function(nrow=1L, ncol=1L, lambda=0.05)
{
    if (!isSingleNumber(nrow) || !isSingleNumber(ncol))
        stop(wmsg("'nrow' and 'ncol' must be single integers"))
    poissonSparseArray(c(nrow, ncol), lambda)
}

