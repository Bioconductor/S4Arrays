### =========================================================================
### Manipulation of array selections
### -------------------------------------------------------------------------
###

### Note that Lindex2Mindex() and Mindex2Lindex() do not allow NAs at the
### moment (they trigger an error) even though R base allows them in matrix
### subsripts.
### TODO: Allow NAs in the input of Lindex2Mindex() and Mindex2Lindex().
### We could either always allow them or add an argument (e.g. 'allowNAs')
### that the user would set to TRUE to allow them.

### Like base::arrayInd() but faster and accepts a matrix for 'dim' (with 1
### row per element in 'Lindex').
Lindex2Mindex <- function(Lindex, dim, use.names=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    ## 'dim' can be a matrix so it's important to use storage.mode()
    ## instead of as.integer().
    if (storage.mode(dim) == "double")
        storage.mode(dim) <- "integer"
    ## 'Lindex' and 'dim' will be fully checked at the C level.
    .Call2("C_Lindex2Mindex", Lindex, dim, use.names, PACKAGE="S4Arrays")
}

Mindex2Lindex <- function(Mindex, dim, use.names=FALSE, as.integer=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.integer))
        stop("'as.integer' must be TRUE or FALSE")
    ## 'dim' and/or 'Mindex' can be matrices so it's important to use
    ## storage.mode() instead of as.integer(). Also, unlike as.integer(),
    ## this preserves the names/dimnames.
    if (storage.mode(dim) == "double")
        storage.mode(dim) <- "integer"
    if (storage.mode(Mindex) == "double")
        storage.mode(Mindex) <- "integer"
    ## 'Mindex' and 'dim' will be fully checked at the C level.
    .Call2("C_Mindex2Lindex", Mindex, dim, use.names,
                              as.integer, PACKAGE="S4Arrays")
}

### NOT exported but used in the SparseArray and DelayedArray packages!
Mindex_order <- function(Mindex)
{
    stopifnot(is.matrix(Mindex), ncol(Mindex) != 0L)
    cols <- lapply(ncol(Mindex):1, function(j) Mindex[ , j])
    do.call(order, cols)
}

