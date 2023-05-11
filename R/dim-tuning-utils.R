### =========================================================================
### Dim tuning utilities
### -------------------------------------------------------------------------
###
### Dim tuning is the act of adding and/or dropping ineffective dimensions
### to/from an array-like object. The exact actions to perform on the
### dimensions of the object are described via the 'dim_tuner' argument.
### See src/dim_tuning_utils.c for more information.


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

### NOT exported but used in the SparseArray and DelayedArray packages!
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

### NOT exported but used in the SparseArray package!
make_dim_tuner_from_old2new_dims <- function(old_dim, new_dim, x_class)
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

