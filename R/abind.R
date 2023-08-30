### =========================================================================
### Bind multidimensional arrays along an arbitrary dimension
### -------------------------------------------------------------------------


.normarg_along <- function(along)
{
    if (!isSingleNumber(along))
        stop(wmsg("'along' must be a single number"))
    if (!is.integer(along))
        along <- as.integer(along)
    if (along <= 0L)
        stop(wmsg("'along' must be a positive integer"))
    along
}

### Return a matrix with one row per dim and one column per object if the
### objects are "bindable". Otherwise return a string describing why they
### are not. This design allows the function to be used in the context of
### a validity method.
get_dims_to_bind <- function(objects, along)
{
    stopifnot(is.list(objects))
    along <- .normarg_along(along)
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    if (any(ndims < along))
        stop(wmsg("the array-like objects to bind must have at least ",
                  along, " dimensions for this binding operation"))
    ndim <- ndims[[1L]]
    if (!all(ndims == ndim))
        return(paste0("all the objects to bind must have ",
                      "the same number of dimensions"))
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return("the objects to bind have no dimensions")
    dims <- matrix(tmp, nrow=ndim)
    tmp <- dims[-along, , drop=FALSE]
    if (!all(tmp == tmp[ , 1L]))
        return("the objects to bind have incompatible dimensions")
    dims
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_dims_along() and combine_dimnames_along()
###
### Compute the expected dims and dimnames of abind()'s result.
### Both combine_dims_along() and combine_dimnames_along() generalize what
### rbind/cbind do in the 2D case to the multidimensional case.
###
### A NOTE ABOUT CRAN PACKAGE abind 1.4-5: By default, abind::abind() does
### NOT follow the same rules as rbind() and cbind() for propagation of the
### dimnames. This is despite its man page (?abind::abind) claiming that it
### does (see documentation of the 'use.first.dimnames' argument), when in
### fact it has it backward! Very misleading!
### The abind0() function defined below in this file is a simple wrapper
### around abind::abind() that fixes that.

### Combine the dims the rbind/cbind way.
combine_dims_along <- function(dims, along)
{
    stopifnot(is.matrix(dims))
    along <- .normarg_along(along)
    stopifnot(along <= nrow(dims))
    ans_dim <- dims[ , 1L]
    ans_dim[[along]] <- sum(dims[along, ])
    ans_dim
}

### Assume all the arrays in 'objects' have the same number of dimensions.
.combine_dimnames <- function(objects)
{
    lapply(seq_along(dim(objects[[1L]])),
        function(n) {
            for (object in objects) {
                dn <- dimnames(object)[[n]]
                if (!is.null(dn))
                    return(dn)
            }
            NULL
        })
}

### Combine the dimnames the rbind/cbind way.
combine_dimnames_along <- function(objects, dims, along)
{
    stopifnot(is.matrix(dims),
              isSingleInteger(along), along >= 1L, along <= nrow(dims))
    dimnames <- .combine_dimnames(objects)
    along_names <- lapply(objects, function(object) dimnames(object)[[along]])
    along_names_lens <- lengths(along_names)
    if (any(along_names_lens != 0L)) {
        fix_idx <- which(along_names_lens != dims[along, ])
        along_names[fix_idx] <- lapply(dims[along, fix_idx], character)
    }
    along_names <- unlist(along_names, use.names=FALSE)
    if (!is.null(along_names))
        dimnames[[along]] <- along_names
    simplify_NULL_dimnames(dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simple_abind()
###

### 'objects' is assumed to be a list of vector-like objects.
### 'nblock' is assumed to be a single integer value (stored as a numeric)
### that is a common divisor of the object lengths.
.intertwine_blocks <- function(objects, nblock, ans_dim)
{
    x0 <- unlist(lapply(objects, `[`, 0L), recursive=FALSE, use.names=FALSE)
    objects_lens <- lengths(objects)
    if (all(objects_lens == 0L))
        return(set_dim(x0, ans_dim))

    idx <- which(vapply(objects,
        function(object) { typeof(object) != typeof(x0) },
        logical(1),
        USE.NAMES=FALSE))
    if (length(idx) != 0L)
        objects[idx] <- lapply(objects[idx], `storage.mode<-`, typeof(x0))

    .Call2("C_abind", objects, nblock, ans_dim, PACKAGE="S4Arrays")
}

### A stripped-down version of abind::abind().
### Some differences:
###   (a) Treatment of dimnames: simple_abind() treatment of dimnames is
###       consistent with base::rbind() and base::cbind(). This is not the
###       case for abind::abind() which does some strange things with the
###       dimnames.
###   (b) Performance: simple_abind() is much faster than abind::abind()
###       (between 3x and 15x). Also note that in the 'along=1L' and 'along=2L'
###       cases, it's generally as fast (and most of the time faster) than
###       base::rbind() and base::cbind().
###       For example, with 'm <- matrix(1:30000000, nrow=5000)',
###       'simple_abind(m, m, m, along=1L)' is 14x faster than
###       'abind::abind(m, m, m, along=1L)' and 11x faster than
###       'base::rbind(m, m, m)'.
###   (c) abind::abind() is broken on matrices of type "list".
simple_abind <- function(..., along)
{
    along <- .normarg_along(along)
    objects <- S4Vectors:::delete_NULLs(list(...))
    if (length(objects) == 0L)
        return(NULL)

    ## Check dim compatibility.
    dims <- get_dims_to_bind(objects, along)
    if (is.character(dims))
        stop(wmsg(dims))
    if (length(objects) == 1L)
        return(objects[[1L]])

    ## Perform the binding.
    nblock <- prod(dims[-seq_len(along), 1L])  # numeric that can be >
                                               # .Machine$integer.max
    ans <- .intertwine_blocks(objects, nblock, combine_dims_along(dims, along))

    ## Combine and set the dimnames.
    set_dimnames(ans, combine_dimnames_along(objects, dims, along))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### simple_abind2()
###
### A wrapper to simple_abind() that adds some of the capabilities originally
### found in abind::abind():
### 1. Not all arrays (supplied via 'objects') are required to have the same
###    number of dimensions. If N is the number of dimensions of the arrays
###    with the most dimensions, then artifical (a.k.a. ineffective) dimensions
###    are added to the arrays with less dimensions so that they also have
###    N dimensions.
### 2. 'along' can be any integer value between 1 and N+1. When set to N+1,
###    then one more artifical dimension is added to all arrays.
### 3. simple_abind2() supports the 'rev.along' argument which can be any
###    value between 0 and N.
### 4. Arguments 'along' and 'rev.along' can be omitted, in which case binding
###    happens along the N-th dimension.

.normarg_rev_along <- function(rev.along)
{
    if (!isSingleNumber(rev.along))
        stop(wmsg("'rev.along' must be a single number"))
    if (!is.integer(rev.along))
        rev.along <- as.integer(rev.along)
    if (rev.along < 0L)
        stop(wmsg("'rev.along' must be a non-negative integer"))
    rev.along
}

### To be re-used by other abind() methods e.g. by the method for SparseArray
### objects.
get_along <- function(N, along=NULL, rev.along=NULL)
{
    if (is.null(along) && is.null(rev.along))
        return(N)
    if (is.null(rev.along)) {
        along <- .normarg_along(along)
        if (along > N + 1L)
            stop(wmsg("'along' must be <= N + 1, where N is the number of ",
                      "dimensions of the arrays with the most dimensions"))
    } else {
        rev.along <- .normarg_rev_along(rev.along)
        if (rev.along > N)
            stop(wmsg("'rev.along' must be <= N, where N is the number of ",
                      "dimensions of the arrays with the most dimensions"))
        along <- N + 1L - rev.along
    }
    along
}

### 'ndim' is expected to be >= number of dimensions of the arrays with
### the most dimensions.
### To be re-used by other abind() methods e.g. by the method for SparseArray
### objects.
add_missing_dims <- function(objects, ndim)
{
    stopifnot(is.list(objects), isSingleInteger(ndim))
    lapply(objects,
        function(object) {
            object_dim <- dim(object)
            object_ndim <- length(object_dim)
            if (object_ndim >= ndim)
                return(object)
            set_dim(object, c(object_dim, rep.int(1L, ndim - object_ndim)))
        }
    )
}

simple_abind2 <- function(objects, along=NULL, rev.along=NULL)
{
    stopifnot(is.list(objects))
    objects <- S4Vectors:::delete_NULLs(objects)
    if (length(objects) == 0L)
        return(NULL)
    ndims <- vapply(objects, function(object) length(dim(object)), integer(1))
    N <- max(ndims)
    along <- get_along(N, along=along, rev.along=rev.along)
    ans_ndim <- max(N, along)
    objects <- add_missing_dims(objects, ans_ndim)
    do.call(simple_abind, c(objects, list(along=along)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### abind0()
###
### A simple wrapper around abind::abind() that uses 'use.first.dimnames=TRUE'
### by default to correct mishandling of the dimnames.
### See "A NOTE ABOUT CRAN PACKAGE abind 1.4-5" above in this file.

abind0 <- function(..., use.first.dimnames=TRUE)
{
    abind::abind(..., use.first.dimnames=use.first.dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The abind() generic and default method
###

### Like in the original abind::abind(), the default for 'along' is the last
### dimension. However, here in the generic function, the default value is
### NULL instead of 'N', but it means the same thing.
setGeneric("abind", signature="...",
    function(..., along=NULL, rev.along=NULL) standardGeneric("abind")
)

### All argument names of abind::abind(), ignoring the ellipsis.
.ORGINAL_ABIND_ARGNAMES <- setdiff(names(formals(abind::abind)), "...")

### Return the supplied objects in an ordinary list. Any argument that is
### not a "recognized" argument of the original abind::abind() is considered
### an object.
.extract_objects <- function(...)
{
    objects <- list(...)
    objects[.ORGINAL_ABIND_ARGNAMES] <- NULL
    unname(S4Vectors:::delete_NULLs(objects))
}

### Return the list of supplied arguments that belong to the original
### abind::abind() interface.
.extract_crazy_args <- function(...)
{
    dots <- list(...)
    dots[intersect(names(dots), .ORGINAL_ABIND_ARGNAMES)]
}

### One object in 'objects' **must** be an Array derivative, but not all!
### Try to re-dispatch on an existing abind() method defined for some Array
### derivative. If such method is found, then it should return an Array
### derivative.
.abind_as_Array <- function(objects, is_Array, is_array,
                            along=NULL, rev.along=NULL)
{
    stopifnot(is.list(objects),
              is.logical(is_Array), length(is_Array) == length(objects),
              is.logical(is_array), length(is_array) == length(objects))
    if (!all(is_Array | is_array))
        stop(wmsg("all objects passed to abind() must be array-like ",
                  "objects when some of them are Array derivatives"))

    ## Avoid infinite recursion if re-dispatch fails to find an abind()
    ## method defined for some Array derivatives.
    if (all(is_Array)) {
        x1 <- objects[[1L]]
        stop(wmsg("no abind() method found for ", class(x1)[[1L]], " objects"))
    }

    ## Coerce all arrays or matrices to the class of the first Array
    ## derivative found in 'objects'.
    x1 <- objects[[which(is_Array)[[1L]]]]
    objects[is_array] <- lapply(objects[is_array], S4Vectors:::coerce2, x1)

    ## Call the abind() generic in the hope that method dispatch will find
    ## a method defined for objects of class 'class(x1)'.
    do.call(S4Arrays::abind, c(objects, list(along=along, rev.along=rev.along)))
}

### .default_abind() tries to dispatch on the following functions, in this
### order:
### 1. Dispatch on the abind() generic: If the objects to bind are supplied
###    in a list then .default_abind() will dispatch on the abind() generic.
###    However this time the objects are passed to abind() as individual
###    arguments. This gives the abind() generic the opportunity to dispatch
###    on the appropriate S4 method (e.g. on the abind() method for SparseArray
###    objects if all the objects are SparseArray derivatives). If no specific
###    S4 method is found then .default_abind() will be called again but note
###    that infinite recursion is naturally avoided.
### 2. Dispatch on .abind_as_Array(): If at least one of the objects to
###    bind is an Array derivative then .default_abind() will dispatch
###    on .abind_as_Array().
### 3. Dispatch on simple_abind2(): When all the objects to bind are ordinary
###    arrays (including matrices) and no "crazy arguments" are supplied,
###    then .default_abind() will dispatch on simple_abind2(). Note that
###    a "crazy argument" is any formal argument located after the ellipsis
###    in abind::abind() except 'along' and 'rev.along'.
### 4. Finally, if all the above fails, then dispatch on abind::abind()
###    wrapper abind0().
.default_abind <- function(..., along=NULL, rev.along=NULL)
{
    objects <- .extract_objects(...)
    if (length(objects) == 0L)
        return(NULL)
    crazy_args <- .extract_crazy_args(...)
    all_extra_args <- c(list(along=along, rev.along=rev.along), crazy_args)
    ok1 <- vapply(objects, is.list, logical(1))
    ok2 <- vapply(objects, is.array, logical(1))
    ok3 <- vapply(objects, is, logical(1), "Array")
    if (any(ok1 & !(ok2 | ok3))) {
        ## 1. Dispatch on the abind() generic.
        if (length(objects) != 1L)
            stop(wmsg("when supplying a list, all the objects ",
                      "to bind must be supplied via the list"))
        x <- objects[[1L]]
        return(do.call(S4Arrays::abind, c(x, all_extra_args)))
    }
    if (any(ok3)) {
        ## 2. Dispatch on .abind_as_Array().
        if (length(crazy_args) != 0L) {
            argnames <- paste0("'", names(crazy_args), "'", collapse=",")
            stop(wmsg("unsupported argument(s) when calling abind() ",
                      "on Array derivatives: ", argnames))
        }
        return(.abind_as_Array(objects, ok3, ok2,
                               along=along, rev.along=rev.along))
    }
    if (all(ok2) && length(crazy_args) == 0L) {
        ## 3. Dispatch on simple_abind2(), which is significantly faster than
        ##    abind::abind().
        return(simple_abind2(objects, along=along, rev.along=rev.along))
    }
    ## 4. Dispatch on abind::abind() wrapper abind0().
    all_extra_args <- S4Vectors:::delete_NULLs(all_extra_args)
    do.call(abind0, c(objects, all_extra_args))
}

setMethod("abind", "ANY", .default_abind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bind arrays along their 1st or 2nd dimension
###
### TODO: No need for these functions to be generics. They should just be
### simple wrappers to 'abind(..., along=1L)' and 'abind(..., along=2L)',
### respectively.

setGeneric("arbind", function(...) standardGeneric("arbind"))
setGeneric("acbind", function(...) standardGeneric("acbind"))

setMethod("arbind", "ANY", function(...) abind(..., along=1L))
setMethod("acbind", "ANY", function(...) abind(..., along=2L))

