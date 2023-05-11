### =========================================================================
### Dim tuning utilities
### -------------------------------------------------------------------------
###
### Dim tuning is the act of adding and/or dropping ineffective dimensions
### to/from an array-like object. The exact actions to perform on the
### dimensions of the object are described via the 'dim_tuner' argument.
### See src/dim_tuning_utils.c or more information.


### NOT exported but used in the SparseArray package!
tune_dims <- function(dim, dim_tuner)
{
    stopifnot(is.integer(dim),
              is.integer(dim_tuner))
    .Call2("C_tune_dims", dim, dim_tuner, PACKAGE="S4Arrays")
}

### NOT exported but used in the SparseArray package!
tune_dimnames <- function(dimnames, dim_tuner)
{
    stopifnot(is.null(dimnames) || is.list(dimnames),
              is.integer(dim_tuner))
    .Call2("C_tune_dimnames", dimnames, dim_tuner, PACKAGE="S4Arrays")
}

