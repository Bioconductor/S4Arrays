### =========================================================================
### write_block()
### -------------------------------------------------------------------------


### 'sink' must be a **writable** array-like object, typically a
### RealizationSink concrete subclass (RealizationSink is a virtual class
### defined in the DelayedArray package), but not necessarily. See default
### write_block() method below.
### Note that for now dispatch is only on the first argument ('sink') but
### we could change that in the future to also dispatch on the third
### argument ('block') if the need arises.
### Must return the modified 'sink'.
setGeneric("write_block", signature="sink",
    function(sink, viewport, block)
    {
        sink_dim <- dim(sink)
        if (is.null(sink_dim))
            stop(wmsg("the first argument to write_block() must be an ",
                      "array-like object (i.e. it must have dimensions)"))
        stopifnot(is(viewport, "ArrayViewport"),
                  identical(refdim(viewport), sink_dim),
                  identical(dim(block), dim(viewport)))
        standardGeneric("write_block")
    }
)

### Based on replace_by_Nindex() which is based on subassignment ('[<-'),
### so work on any array-like object 'sink' that supports subassignment.
### Thanks to this method, write_block() will work out-of-the-box on an
### ordinary array and other in-memory array-like object that supports
### subassignment (e.g. SparseArray object from the SparseArray package
### or sparseMatrix derivative from the Matrix package).
setMethod("write_block", "ANY",
    function(sink, viewport, block)
    {
        if (is.array(sink)) {
            ## Subassignment of an ordinary array only works if the right
            ## value is also an ordinary array.
            if (!is.array(block))
                block <- as.array(block)
        } else if (is(sink, "sparseMatrix")) {
            ## Subassignment of a sparseMatrix derivative (e.g. dgCMatrix
            ## object) only works if the right value is a sparseMatrix
            ## derivative or ordinary array.
            if (!(is(block, "sparseMatrix") || is.array(block))) {
                if (is_sparse(block))
                    block <- as(block, "sparseMatrix")
                else
                    block <- as.array(block)
            }
        }
        Nindex <- makeNindexFromArrayViewport(viewport)
        replace_by_Nindex(sink, Nindex, block)
    }
)

