### =========================================================================
### type()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Default type() method
###

### Perform 'extract_array(x, list(integer(0), ..., integer(0)))'.
### 'x' is **trusted** to be an array-like object.
extract_empty_array <- function(x)
{
    index <- rep.int(list(integer(0)), length(dim(x)))
    extract_array(x, index)
}

### The default type() method (the generic is defined in BiocGenerics)
### implements the 'typeof(as.array(x))' semantic. It will work out-of-the-box
### on any array-like object that supports extract_array().
setMethod("type", "ANY",
    function(x)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("the default type() method only supports ",
                      "array-like objects"))
        type(extract_empty_array(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() method for DataFrame objects
###

.get_DataFrame_type <- function(x)
{
    ## Make sure that this remains consistent with .extract_DataFrame_slice0()
    ## defined in R/extract_array.R
    x0 <- x[0L, , drop=FALSE]
    df0 <- as.data.frame(x0)
    if (ncol(df0) != ncol(x))
        stop(wmsg("the type() method for DataFrame objects only suports ",
                  "objects for which 'as.data.frame(x)' preserves the ",
                  "number of columns"))
    type(df0)
}

### Equivalent to 'typeof(as.matrix(as.data.frame(x)))' (granted
### that 'as.data.frame(x)' preserves the number of columns) but with an
### almost-zero memory footprint (it avoids the cost of turning 'x' first
### into a data frame then into a matrix).
setMethod("type", "DataFrame", .get_DataFrame_type)

