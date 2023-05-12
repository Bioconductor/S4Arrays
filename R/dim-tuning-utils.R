### =========================================================================
### Dim tuning utilities
### -------------------------------------------------------------------------
###
### "Dim tuning" is the act of adding and/or dropping "ineffective
### dimensions" to/from an array-like object, typically via the drop()
### and/or dim() setter. The exact transformation to operate on the vector
### of dimensions of the object can be precisely described by supplying
### a 'dim_tuner' vector.
### See src/dim_tuning_utils.c for additional information.


### NOT exported but used in the DelayedArray packages!
normalize_dim_replacement_value <- function(value, x_dim)
{
    if (is.null(value))
        stop(wmsg("you can't do that, sorry"))
    if (!is.numeric(value))
        stop(wmsg("the supplied dim vector must be numeric"))
    if (length(value) == 0L)
        stop(wmsg("the supplied dim vector cannot be empty"))
    if (!is.integer(value))
        value <- as.integer(value)
    if (S4Vectors:::anyMissingOrOutside(value, 0L))
        stop(wmsg("the supplied dim vector cannot contain negative ",
                  "or NA values"))
    prod1 <- prod(value)
    prod2 <- prod(x_dim)
    if (prod1 != prod2)
        stop(wmsg("the supplied dims [product ", prod1, "] do not match ",
                  "the length of object [", prod2, "]"))
    value
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The tune_Array_dims() low-level generic
###
### Array derivatives (e.g. SparseArray or DelayedArray objects) only need
### to implement a tune_Array_dims() method to have drop() and the dim()
### setter work out-of-the-box.
###
### Note that a "dim tuning" operation does NOT change the length of the
### object (which is prod(dim(x))) or alter its content, and should always
### be reversible (except when it drops ineffective dimensions with dimnames
### on them). To revert a "dim tuning" operation, simply tune again
### with '-dim_tuner' (i.e. with minus 'dim_tuner'). More precisely, for
### tune_Array_dims(), 'x2' should always be identical to 'x' here:
###
###     y <- tune_Array_dims(x, dim_tuner)
###     x2 <- tune_Array_dims(y, -dim_tuner)
###     identical(x2, x)  # should be TRUE
###
### This should be the case for any array-like object 'x' (with no dimnames
### on its ineffective dimensions) and any 'dim_tuner' vector compatible
### with 'dim(x)'.

setGeneric("tune_Array_dims", signature="x",
    function(x, dim_tuner) standardGeneric("tune_Array_dims")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tune_dims() and tune_dimnames()
###
### NOT exported but used by the tune_Array_dims() method for SVT_SparseArray
### objects defined in the SparseArray package.

tune_dims <- function(dim, dim_tuner)
{
    stopifnot(is.integer(dim),
              is.integer(dim_tuner))
    .Call2("C_tune_dims", dim, dim_tuner, PACKAGE="S4Arrays")
}

tune_dimnames <- function(dimnames, dim_tuner)
{
    stopifnot(is.null(dimnames) || is.list(dimnames),
              is.integer(dim_tuner))
    .Call2("C_tune_dimnames", dimnames, dim_tuner, PACKAGE="S4Arrays")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop() method
###

### Expected to be semantically equivalent to 'drop(as.array(x))'.
### Will work out-of-the-box on any Array derivative that supports
### tune_Array_dims() and as.array(). Note that the latter is used
### only if 'x' has at most one effective dimension.
.drop_Array <- function(x)
{
    is_effective <- dim(x) != 1L
    if (sum(is_effective) <= 1L)
        return(drop(as.array(x)))  # ordinary vector
    dim_tuner <- -as.integer(!is_effective)
    tune_Array_dims(x, dim_tuner)
}

setMethod("drop", "Array", .drop_Array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() setter
###

.diff_dims <- function(old_dim, new_dim, x_class)
{
    stopifnot(is.integer(old_dim), is.integer(new_dim))

    cannot_map_msg <- c(
        "Cannot map the supplied dim vector to the current dimensions of ",
        "the object. On a ", x_class, " object, the dim() setter can only ",
        "be used to drop and/or add \"ineffective dimensions\" (i.e. ",
        "dimensions equal to 1) to the object."
    )
    can_map_effective_dimensions <- function(effdim_idx1, effdim_idx2) {
        if (length(effdim_idx1) != length(effdim_idx2))
            return(FALSE)
        if (length(effdim_idx1) == 0L)
            return(TRUE)
        all(old_dim[effdim_idx1] == new_dim[effdim_idx2])
    }

    ## Get index of old and new effective dimensions.
    effdim_idx1 <- which(old_dim != 1L)
    effdim_idx2 <- which(new_dim != 1L)

    if (!can_map_effective_dimensions(effdim_idx1, effdim_idx2))
        stop(wmsg(cannot_map_msg))

    compute_dim_tuner <- function(effdim_idx1, effdim_idx2) {
        idx1 <- c(effdim_idx1, length(old_dim) + 1L)
        idx2 <- c(effdim_idx2, length(new_dim) + 1L)
        diffs1 <- S4Vectors:::diffWithInitialZero(idx1)
        diffs2 <- S4Vectors:::diffWithInitialZero(idx2)
        deltas <- pmax(diffs1, diffs2)
        nonzero_runlengths <- deltas - 1L
        ans_len <- sum(nonzero_runlengths) + length(effdim_idx1)
        ans <- integer(ans_len)
        offsets <- c(0L, head(cumsum(deltas), n=-1L))
        for (k in seq_along(idx1)) {
            d <- diffs2[[k]] - diffs1[[k]]
            op <- ifelse(d > 0L, 1L, -1L)
            ans[offsets[[k]] + seq_len(abs(d))] <- rep.int(op, abs(d))
        }
        ans
    }

    compute_dim_tuner(effdim_idx1, effdim_idx2)
}

.set_Array_dim <- function(x, value)
{
    x_dim <- dim(x)
    value <- normalize_dim_replacement_value(value, x_dim)
    dim_tuner <- .diff_dims(x_dim, value, class(x))
    ans <- tune_Array_dims(x, dim_tuner)
    stopifnot(identical(dim(ans), value))  # sanity check
    ans
}

setReplaceMethod("dim", "Array", .set_Array_dim)

