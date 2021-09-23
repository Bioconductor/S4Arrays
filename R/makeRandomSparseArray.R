### =========================================================================
### makeRandomSparseArray()
### -------------------------------------------------------------------------


makeRandomSparseArray <- function(dim=1L, density=0.05)
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
    #TODO!

    ans
}

makeRandomSparseMatrix <- function(nrow=1L, ncol=1L, density=0.05)
{
    if (!isSingleNumber(nrow) || !isSingleNumber(ncol))
        stop(wmsg("'nrow' and 'ncol' must be single integers"))
    makeRandomSparseArray(c(nrow, ncol), density)
}

