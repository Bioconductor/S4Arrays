### =========================================================================
### ArrayGrid objects
### -------------------------------------------------------------------------


### An ArrayGrid object represents a grid on top of an array (called "the
### reference array" or "the underlying array"). The ArrayGrid class is a
### virtual class with 3 concrete subclasses, ArbitraryArrayGrid and
### RegularArrayGrid, for representing an arbitrarily-spaced or a
### regularly-spaced grid, respectively. The API we implement on these objects
### is divided into 3 groups of methods:
###   1) One special method:
###      - refdim(x): Return the dimensions of the reference array.
###   2) Methods from the array API (an ArrayGrid object can be seen as an
###      array of ArrayViewport objects that cover the reference array without
###      overlapping with each others):
###      - dim(x): Return the number of grid elements (i.e. viewports) along
###        each dimension of the reference array.
###      - x[[i_1, i_2, ..., i_n]]: Multidimensional double bracket
###        subsetting. Return an ArrayViewport object.
###   3) Methods from the list API:
###      - length(): Return the total number of grid elements.
###      - x[[i]]: Linear double bracket subsetting. Return an ArrayViewport
###        object.
### Groups 2) and 3) give these objects 2 semantics: array-like and list-like.
### Note that length() and "linear double bracket subsetting" are consistent
### with dim() and "multidimensional double bracket subsetting", respectively.
### So the array-like and list-like semantics are compatible.
###

setClass("ArrayGrid",
    contains=c("Array", "List"),
    representation("VIRTUAL"),
    prototype(elementType="ArrayViewport")
)

setClass("DummyArrayGrid",
    contains="ArrayGrid",
    representation(
        refdim="integer"    # Dimensions of the reference array.
    ),
    prototype(elementType="DummyArrayViewport")
)

setClass("ArbitraryArrayGrid",
    contains="ArrayGrid",
    representation(
        tickmarks="list"    # A list of integer vectors, one along each
                            # dimension of the reference array,
                            # representing the tickmarks along that
                            # dimension. Each integer vector must be sorted
                            # in ascending order.
    ),
    prototype(elementType="SafeArrayViewport")
)

setClass("RegularArrayGrid",
    contains="ArrayGrid",
    representation(
        refdim="integer",   # Dimensions of the reference array.
        spacings="integer"  # Grid spacing along each dimension.
    ),
    prototype(elementType="SafeArrayViewport")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

.prod2 <- function(x)
{
    p <- prod(x)
    ## We only use .prod2() in a context where 'x' is guaranteed to contain
    ## non-negative values so 'p' will always be >= 0.
    if (p <= .Machine$integer.max) as.integer(p) else p
}

.get_DummyArrayGrid_spacings_along <- function(x, along) x@refdim[[along]]

.get_ArbitraryArrayGrid_spacings_along <- function(x, along)
    S4Vectors:::diffWithInitialZero(x@tickmarks[[along]])

.get_ArbitraryArrayGrid_max_spacings <- function(x)
{
    vapply(seq_along(x@tickmarks),
        function(along)
            max(0L, .get_ArbitraryArrayGrid_spacings_along(x, along)),
        integer(1),
        USE.NAMES=FALSE
    )
}

### Get length of biggest viewport in ArbitraryArrayGrid object 'x'.
.get_ArbitraryArrayGrid_maxlength <- function(x)
{
    .prod2(.get_ArbitraryArrayGrid_max_spacings(x))
}

### Get length of biggest viewport in RegularArrayGrid object 'x'.
.get_RegularArrayGrid_maxlength <- function(x) .prod2(x@spacings)

get_RegularArrayGrid_dim <- function(refdim, spacings)
{
    ans <- refdim %/% spacings + (refdim %% spacings != 0L)
    ans[is.na(ans)] <- 1L
    ans
}

.get_RegularArrayGrid_spacings_along <- function(x, along)
{
    D <- x@refdim[[along]]
    if (D == 0L)
        return(0L)
    spacing <- x@spacings[[along]]
    ans <- rep.int(spacing, D %/% spacing)
    r <- D %% spacing
    if (r != 0L)
        ans <- c(ans, r)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_DummyArrayGrid <- function(x)
{
    msg <- validate_dim_slot(x, "refdim")
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("DummyArrayGrid", .validate_DummyArrayGrid)

.valid_tickmarks <- function(tm)
{
    is.integer(tm) && !S4Vectors:::anyMissingOrOutside(tm, 0L) && isSorted(tm)
}
.validate_ArbitraryArrayGrid <- function(x)
{
    x_tickmarks <- x@tickmarks
    if (!is.list(x_tickmarks))
        return("'tickmarks' slot must be a list")
    ok <- vapply(x_tickmarks, .valid_tickmarks, logical(1), USE.NAMES=FALSE)
    if (!all(ok))
        return(paste0("each list element in 'tickmarks' slot must be a ",
                      "sorted integer vector of non-negative values"))
    x_maxlen <- .get_ArbitraryArrayGrid_maxlength(x)
    if (x_maxlen > .Machine$integer.max)
        return(paste0("grid is too coarse (all grid elements must have a ",
                      "length <= .Machine$integer.max)"))
    TRUE
}
setValidity2("ArbitraryArrayGrid", .validate_ArbitraryArrayGrid)

.validate_RegularArrayGrid <- function(x)
{
    msg <- validate_dim_slot(x, "refdim")
    if (!isTRUE(msg))
        return(msg)
    msg <- validate_dim_slot(x, "spacings")
    if (!isTRUE(msg))
        return(msg)
    x_spacings <- x@spacings
    x_refdim <- x@refdim
    if (length(x_spacings) != length(x_refdim))
        return("'spacings' and 'refdim' slots must have the same length")
    if (!all(x_spacings <= x_refdim))
        return(paste0("values in 'spacings' slot must be <= their ",
                      "corresponding value in 'refdim' slot"))
    if (any(x_spacings == 0L & x_refdim != 0L))
        return(paste0("values in 'spacings' slot cannot be 0 unless their ",
                      "corresponding value in 'refdim' slot is also 0"))
    x_maxlen <- .get_RegularArrayGrid_maxlength(x)
    if (x_maxlen > .Machine$integer.max)
        return(paste0("grid is too coarse (all grid elements must have a ",
                      "length <= .Machine$integer.max)"))
    TRUE
}
setValidity2("RegularArrayGrid", .validate_RegularArrayGrid)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("refdim", "DummyArrayGrid", function(x) x@refdim)

setMethod("refdim", "ArbitraryArrayGrid",
    function(x)
    {
        mapply(
            function(tm, tm_len) if (tm_len == 0L) 0L else tm[[tm_len]],
            x@tickmarks,
            lengths(x@tickmarks),
            USE.NAMES=FALSE
        )
    }
)

setMethod("refdim", "RegularArrayGrid", function(x) x@refdim)

setMethod("dim", "DummyArrayGrid", function(x) rep.int(1L, length(refdim(x))))

setMethod("dim", "ArbitraryArrayGrid", function(x) lengths(x@tickmarks))

setMethod("dim", "RegularArrayGrid",
    function(x) get_RegularArrayGrid_dim(refdim(x), x@spacings)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

DummyArrayGrid <- function(refdim)
{
    if (!is.numeric(refdim))
        stop(wmsg("'refdim' must be an integer vector"))
    if (!is.integer(refdim))
        refdim <- as.integer(refdim)
    new2("DummyArrayGrid", refdim=refdim, check=TRUE)
}

ArbitraryArrayGrid <- function(tickmarks)
{
    if (!is.list(tickmarks))
        stop(wmsg("'tickmarks' must be a list"))
    new2("ArbitraryArrayGrid", tickmarks=tickmarks, check=TRUE)
}

### Note that none of the dimensions of an RegularArrayGrid object can be 0,
### even when some dimensions of the reference array are 0 (in which case,
### the corresponding dimensions of the grid object are set to 1). As a
### consequence, an RegularArrayGrid object always contains at least 1 grid
### element. Each dimension of the first grid element is always equal to the
### spacing along that dimension i.e. for any RegularArrayGrid object,
### 'dim(grid[[1]])' is identical to 'spacings'.
### If 'spacings' is omitted, return a grid with a single grid element
### covering the whole reference array.
RegularArrayGrid <- function(refdim, spacings=refdim)
{
    if (!is.numeric(refdim))
        stop(wmsg("'refdim' must be an integer vector"))
    if (!is.integer(refdim))
        refdim <- as.integer(refdim)
    if (!is.numeric(spacings))
        stop(wmsg("'spacings' must be an integer vector"))
    if (!is.integer(spacings))
        spacings <- as.integer(spacings)
    if (length(refdim) != length(spacings))
        stop(wmsg("'refdim' and 'spacings' must have the same length"))
    new2("RegularArrayGrid", refdim=refdim, spacings=spacings, check=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [[
###

setMethod("getArrayElement", "DummyArrayGrid",
    function(x, subscripts)
    {
        stopifnot(is.integer(subscripts),
                  identical(subscripts, rep.int(1L, length(subscripts))))
        DummyArrayViewport(refdim(x))
    }
)

setMethod("getArrayElement", "ArbitraryArrayGrid",
    function(x, subscripts)
    {
        stopifnot(is.integer(subscripts))
        x_refdim <- refdim(x)
        ans_end <- mapply(`[[`, x@tickmarks, subscripts, USE.NAMES=FALSE)
        ans_width <- mapply(
            function(along, i)
                .get_ArbitraryArrayGrid_spacings_along(x, along)[[i]],
            seq_along(x_refdim),
            subscripts,
            USE.NAMES=FALSE
        )
        ans_ranges <- IRanges(end=ans_end, width=ans_width)
        ArrayViewport(x_refdim, ans_ranges)
    }
)

setMethod("getArrayElement", "RegularArrayGrid",
    function(x, subscripts)
    {
        stopifnot(is.integer(subscripts))
        x_refdim <- refdim(x)
        ans_offset <- (subscripts - 1L) * x@spacings
        ans_end <- pmin(ans_offset + x@spacings, refdim(x))
        ans_ranges <- IRanges(start=ans_offset + 1L, end=ans_end)
        ArrayViewport(x_refdim, ans_ranges)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dims(), lengths(), and maxlength()
###

### NOT exported.
setGeneric("get_spacings_along", signature="x",
    function(x, along) standardGeneric("get_spacings_along")
)
setMethod("get_spacings_along", "DummyArrayGrid",
    .get_DummyArrayGrid_spacings_along
)
setMethod("get_spacings_along", "ArbitraryArrayGrid",
    .get_ArbitraryArrayGrid_spacings_along
)
setMethod("get_spacings_along", "RegularArrayGrid",
    .get_RegularArrayGrid_spacings_along
)

### Equivalent to 't(vapply(x, dim, refdim(x)))' but faster.
setMethod("dims", "ArrayGrid",
    function(x)
    {
        ans <- as.matrix(get_spacings_along(x, 1L))
        x_ndim <- length(refdim(x))
        if (x_ndim >= 2L) {
            for (along in 2:x_ndim) {
                spacings_along <- get_spacings_along(x, along)
                ans <- cbind(
                    S4Vectors:::rep.int_along_ROWS(ans, length(spacings_along)),
                    rep(spacings_along, each=nrow(ans))
                )
            }
        }
        ans
    }
)

### Equivalent to 'vapply(x, length, integer(1))' or to 'rowProds(dims(x))'
### but faster.
### The sum of the hyper-volumes of all the grid elements should be equal
### to the hyper-volume of the reference array.
### More concisely: sum(lengths(x)) should be equal to 'prod(refdim(x))'.
setMethod("lengths", "ArrayGrid",
    function(x, use.names=TRUE)
    {
        ans <- get_spacings_along(x, 1L)
        x_ndim <- length(refdim(x))
        if (x_ndim >= 2L) {
            for (along in 2:x_ndim) {
                spacings_along <- get_spacings_along(x, along)
                ans <- ans * rep(spacings_along, each=length(ans))
            }
        }
        ans
    }
)

setMethod("lengths", "DummyArrayGrid",
    function(x, use.names=TRUE) .prod2(refdim(x))
)

### Equivalent to 'max(lengths(x))' except that when 'x' is an ArrayGrid
### object this can be computed without computing 'lengths(x)' first so is
### very efficient.
setGeneric("maxlength", function(x) standardGeneric("maxlength"))
setMethod("maxlength", "ANY", function(x) max(lengths(x)))
setMethod("maxlength", "ArbitraryArrayGrid", .get_ArbitraryArrayGrid_maxlength)
setMethod("maxlength", "RegularArrayGrid", .get_RegularArrayGrid_maxlength)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

### S3/S4 combo for as.character.ArrayGrid
.as.character.ArrayGrid <- function(x, collapse=",", with.brackets=FALSE)
{
    data <- vapply(x,
        function(viewport)
            make_string_from_ArrayViewport(viewport,
                                           collapse=collapse,
                                           with.brackets=with.brackets),
        character(1),
        USE.NAMES=FALSE
    )
    array(data, dim(x))
}
as.character.ArrayGrid <- function(x, ...) .as.character.ArrayGrid(x, ...)
setMethod("as.character", "ArrayGrid", .as.character.ArrayGrid)

setMethod("show", "ArrayGrid",
    function(object)
    {
        if (!is(object, "DummyArrayGrid")) {
            dim_in1string <- paste0(dim(object), collapse=" x ")
            cat(dim_in1string, " ")
        }
        refdim_in1string <- paste0(refdim(object), collapse=" x ")
        cat(class(object), " object on a ", refdim_in1string, " array:\n",
            sep="")
        ## Turn 'object' into a character array.
        print(as.character(object, with.brackets=TRUE),
              quote=FALSE, right=TRUE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Like aperm2() defined in this package, the methods below extend
### base::aperm() by allowing dropping and/or adding ineffective dimensions.
###
### NOTE: The methods below use normarg_perm() which doesn't perform a
### full check of 'perm' e.g. it allows duplicates in it (see normarg_perm()
### in aperm2.R). Instead we call replaceSlots() with 'check=TRUE' to trigger
### validation of the modified ArrayGrid object with the expectation that it
### will fail if for example 'perm' contains NAs or values > length(dim(a)).

.aperm.DummyArrayGrid <- function(a, perm)
{
    perm <- normarg_perm(perm, dim(a))
    ans_refdim <- a@refdim[perm]
    ans_refdim[is.na(perm)] <- 1L
    ## See above NOTE for why it's important to call replaceSlots()
    ## with 'check=TRUE'.
    BiocGenerics:::replaceSlots(a, refdim=ans_refdim, check=TRUE)
}
### S3/S4 combo for aperm.DummyArrayGrid
aperm.DummyArrayGrid <-
    function(a, perm, ...) .aperm.DummyArrayGrid(a, perm, ...)
setMethod("aperm", "DummyArrayGrid", aperm.DummyArrayGrid)

.aperm.ArbitraryArrayGrid <- function(a, perm)
{
    perm <- normarg_perm(perm, dim(a))
    ans_tickmarks <- a@tickmarks[perm]
    ans_tickmarks[is.na(perm)] <- list(1L)
    ## See above NOTE for why it's important to call replaceSlots()
    ## with 'check=TRUE'.
    BiocGenerics:::replaceSlots(a, tickmarks=ans_tickmarks, check=TRUE)
}
### S3/S4 combo for aperm.ArbitraryArrayGrid
aperm.ArbitraryArrayGrid <-
    function(a, perm, ...) .aperm.ArbitraryArrayGrid(a, perm, ...)
setMethod("aperm", "ArbitraryArrayGrid", aperm.ArbitraryArrayGrid)

.aperm.RegularArrayGrid <- function(a, perm)
{
    perm <- normarg_perm(perm, dim(a))
    ans_refdim <- a@refdim[perm]
    ans_refdim[is.na(perm)] <- 1L
    ans_spacings <- a@spacings[perm]
    ans_spacings[is.na(perm)] <- 1L
    ## See above NOTE for why it's important to call replaceSlots()
    ## with 'check=TRUE'.
    BiocGenerics:::replaceSlots(a, refdim=ans_refdim,
                                   spacings=ans_spacings,
                                   check=TRUE)
}
### S3/S4 combo for aperm.RegularArrayGrid
aperm.RegularArrayGrid <-
    function(a, perm, ...) .aperm.RegularArrayGrid(a, perm, ...)
setMethod("aperm", "RegularArrayGrid", aperm.RegularArrayGrid)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### downsample() generic and methods
###
### Reduce the "resolution" of a grid by the specified ratio.
### Act as an endomorphism.
###

setGeneric("downsample", signature="x",
    function(x, ratio=1L) standardGeneric("downsample")
)

.normarg_ratio <- function(ratio, x_dim)
{
    if (!is.numeric(ratio))
        stop(wmsg("'ratio' must be an integer vector"))
    if (!is.integer(ratio))
        ratio <- as.integer(ratio)
    ndim <- length(x_dim)
    if (length(ratio) != 1L && length(ratio) != ndim)
        stop(wmsg("'length(ratio)' must be 1 or the sane as 'length(dim(x))'"))
    if (S4Vectors:::anyMissingOrOutside(ratio, 0L))
        stop(wmsg("'ratio' cannot contain negative or NA values"))
    if (length(ratio) != ndim)
        ratio <- rep.int(ratio, ndim)
    if (any(ratio == 0L & x_dim != 0L))
        stop(wmsg("values in 'ratio' cannot be 0 unless their ",
                  "corresponding dimension in 'x' is also 0"))
    ratio
}

setMethod("downsample", "ArbitraryArrayGrid",
    function(x, ratio=1L)
    {
        ratio <- .normarg_ratio(ratio, dim(x))
        ans_tickmarks <- mapply(
            function(tm, tm_len, r) {
                if (tm_len == 0L)
                    return(integer(0))
                tm[seq2(tm_len, r)]
            },
            x@tickmarks,
            lengths(x@tickmarks),
            ratio,
            SIMPLIFY=FALSE,
            USE.NAMES=FALSE
        )
        ArbitraryArrayGrid(ans_tickmarks)
    }
)

setMethod("downsample", "RegularArrayGrid",
    function(x, ratio=1L)
    {
        ratio <- .normarg_ratio(ratio, dim(x))
        x_refdim <- refdim(x)
        ## We turn 'ratio' into a double vector to prevent a potential
        ## integer overflow.
        ans_spacings <- pmin(x@spacings * as.double(ratio), x_refdim)
        RegularArrayGrid(x_refdim, ans_spacings)
    }
)

