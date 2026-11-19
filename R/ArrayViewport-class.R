### =========================================================================
### ArrayViewport objects
### -------------------------------------------------------------------------


setClass("ArrayViewport",
    contains="Array",
    representation(
        "VIRTUAL",
        refdim="integer"  # Dimensions of "the reference array" i.e. the
                          # array on top of which the viewport is defined
                          # (a.k.a. "the underlying array").
    )
)

### Represents a viewport that covers the whole reference array i.e. that
### has the same dimensions as the reference array.
### IMPORTANT NOTE: Unlike SafeArrayViewport objects below, the length
### of a DummyArrayViewport object (i.e. the number of array elements
### in the viewport) is allowed to be > .Machine$integer.max. This makes
### these objects **unsafe** to use in some contexts. For example a
### DummyArrayViewport object can be used to extract a sparse block
### that is too big to be handled safely:
###   library(Matrix)
###   m0 <- sparseMatrix(i=1, j=1, x=8, dims=c(5e4, 5e4))
###   block0 <- read_block(m0, DummyArrayViewport(x), as.sparse=TRUE)
###   m <- as(block0, "dgCMatrix")
###   m + 1
###   #Error in asMethod(object) :
###   #  Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105
setClass("DummyArrayViewport", contains="ArrayViewport")

### Represents an arbitrary viewport on the reference array with the
### restriction that the length of the viewport (i.e. the number of
### array elements in the viewport) must be <= .Machine$integer.max.
### Unlike DummyArrayViewport objects above, this restriction on their
### length makes these objects always **safe** to use.
### Note that we don't extend the IRanges class because we don't want
### to inherit the full Ranges API (most operations in that API would
### do the wrong thing on SafeArrayViewport objects).
setClass("SafeArrayViewport",
    contains="ArrayViewport",
    representation(
        ranges="IRanges"  # One range per dimension.
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_ArrayViewport <- function(x)
{
    msg <- validate_dim_slot(x, "refdim")
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("ArrayViewport", .validate_ArrayViewport)

.validate_ranges_slot <- function(x)
{
    x_ranges <- x@ranges
    x_refdim <- x@refdim
    if (length(x_ranges) != length(x_refdim))
        return("'ranges' and 'refdim' slots must have the same length")

    ## Check that the viewport is contained in the reference array.
    x_start <- start(x_ranges)
    x_end <- end(x_ranges)
    if (!(all(x_start >= 1L) && all(x_end <= x_refdim)))
        return(paste0("object represents a viewport that is not ",
                      "within the bounds of the reference array"))

    ## A viewport cannot be longer than 2^31 - 1.
    x_dim <- width(x_ranges)
    if (prod(x_dim) > .Machine$integer.max)
        return("a viewport cannot be longer than .Machine$integer.max")

    TRUE
}
.validate_SafeArrayViewport <- function(x)
{
    msg <- .validate_ranges_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("SafeArrayViewport", .validate_SafeArrayViewport)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setGeneric("refdim", function(x) standardGeneric("refdim"))

setMethod("refdim", "ArrayViewport", function(x) x@refdim)

.dims2ranges <- function(dim) IRanges(rep.int(1L, length(dim)), dim)

setMethod("ranges", "DummyArrayViewport", function(x) .dims2ranges(refdim(x)))
setMethod("ranges", "SafeArrayViewport", function(x) x@ranges)

setMethod("start", "ArrayViewport", function(x) start(ranges(x)))

setMethod("width", "ArrayViewport", function(x) width(ranges(x)))

setMethod("end", "ArrayViewport", function(x) end(ranges(x)))

### 'width(x)' and 'dim(x)' are synonyms.
setMethod("dim", "ArrayViewport", function(x) width(ranges(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

DummyArrayViewport <- function(refdim)
    new2("DummyArrayViewport", refdim=refdim, check=TRUE)

### If 'ranges' is omitted, returns a viewport that covers the whole
### reference array.
ArrayViewport <- function(refdim, ranges=NULL)
{
    if (is.null(ranges))
        ranges <- .dims2ranges(refdim)
    new2("SafeArrayViewport", refdim=refdim, ranges=ranges, check=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

setMethod("classNameForDisplay", "ArrayViewport", function(x) "ArrayViewport")

make_string_from_ArrayViewport <- function(viewport, dimnames=NULL,
                                           as.2Dslice=FALSE,
                                           collapse=",", with.brackets=FALSE)
{
    if (!isTRUEorFALSE(as.2Dslice))
        stop("'as.2Dslice' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.brackets))
        stop("'with.brackets' must be TRUE or FALSE")

    viewport_ranges <- ranges(viewport)
    ans <- as.character(viewport_ranges)

    ## Place "blank" subscripts.
    viewport_dim <- dim(viewport)
    viewport_refdim <- refdim(viewport)
    useblank <- viewport_dim == viewport_refdim
    ndim <- length(viewport_dim)
    if (as.2Dslice && ndim >= 3L)
        useblank[3:ndim] <- FALSE
    ans[useblank] <- ""

    if (!is.null(dimnames)) {
        stopifnot(is.list(dimnames), length(dimnames) == ndim)
        usename_idx <- which(!useblank &
                             viewport_dim == 1L &
                             lengths(dimnames) != 0L)
        ans[usename_idx] <- unlist(mapply(`[`,
                                          dimnames[usename_idx],
                                          start(viewport_ranges)[usename_idx],
                                          SIMPLIFY=FALSE,
                                          USE.NAMES=FALSE))
    }
    if (with.brackets)
        ans <- lapply(ans, function(x) if (x == "") " " else x)
    ans <- paste0(ans, collapse=collapse)
    if (with.brackets)
        ans <- paste0("[", ans, "]")
    ans
}

setMethod("show", "ArrayViewport",
    function(object)
    {
        dim_in1string <- paste0(dim(object), collapse=" x ")
        refdim_in1string <- paste0(refdim(object), collapse=" x ")
        cat(dim_in1string, " ", classNameForDisplay(object), " object ",
            "on a ", refdim_in1string, " array: ", sep="")
        s <- make_string_from_ArrayViewport(object, with.brackets=TRUE)
        cat(s, "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeNindexFromArrayViewport()
###

makeNindexFromArrayViewport <- function(viewport, expand.RangeNSBS=FALSE)
{
    viewport_ranges <- ranges(viewport)
    viewport_dim <- dim(viewport)
    viewport_refdim <- refdim(viewport)
    ndim <- length(viewport_dim)
    Nindex <- vector("list", length=ndim)
    is_not_missing <- viewport_dim < viewport_refdim
    if (expand.RangeNSBS) {
        expand_idx <- which(is_not_missing)
    } else {
        viewport_starts <- start(viewport_ranges)
        viewport_ends <- end(viewport_ranges)
        is_width1 <- viewport_dim == 1L
        expand_idx <- which(is_not_missing & is_width1)
        RangeNSBS_idx <- which(is_not_missing & !is_width1)
        Nindex[RangeNSBS_idx] <- lapply(RangeNSBS_idx,
            function(i) {
                range_start <- viewport_starts[[i]]
                range_end <- viewport_ends[[i]]
                upper_bound <- viewport_refdim[[i]]
                new2("RangeNSBS", subscript=c(range_start, range_end),
                                  upper_bound=upper_bound,
                                  check=FALSE)
            }
        )
    }
    if (length(expand_idx) != 0L)
        Nindex[expand_idx] <- as.list(as(viewport_ranges[expand_idx],
                                         "CompressedIntegerList"))
    Nindex
}

