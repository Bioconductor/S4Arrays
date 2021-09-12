### =========================================================================
### readSparseCSV()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readSparseCSV()
###

.scan_first_two_lines <- function(filepath, sep=",")
{
    con <- file(filepath, "r")
    on.exit(close(con))
    line1 <- scan(con, what=character(), sep=sep, nlines=1L, quiet=TRUE)
    line2 <- scan(con, what=character(), sep=sep, nlines=1L, quiet=TRUE)
    list(line1, line2)
}

.looks_like_a_name <- function(x)
    nzchar(x) & is.na(suppressWarnings(as.numeric(x)))

### Guess by looking at the first 2 lines in the file.
### Not ready (this is a mess!)
.guess_dimnames_and_ncol <- function(line1, line2, rownames=NA, colnames=NA)
{
    if (!(is.logical(rownames) && length(rownames) == 1L))
        stop(wmsg("'rownames' must be a single logical value"))
    if (!(is.logical(colnames) && length(colnames) == 1L))
        stop(wmsg("'colnames' must be a single logical value"))
    n1 <- length(line1)
    n2 <- length(line2)
    if (n2 == 0L) {
        if (n1 == 0L)
            stop(wmsg("invalid file: first two lines are empty"))
        if (isTRUE(rownames) && isTRUE(colnames))
            stop(wmsg("file does not seem to contain both rownames and ",
                      "colnames (2nd line is empty)"))
        if (isTRUE(rownames))
            return(list(c(TRUE, FALSE), n1 - 1L))
        if (isTRUE(colnames))
            return(list(c(FALSE, TRUE), n1))
        return(list(c(FALSE, FALSE), n1))
    }
    if (n1 == n2 - 1L) {
        if (isFALSE(rownames) || isFALSE(colnames))
            stop(wmsg("file seems to contain both rownames and colnames ",
                      "(1st line contains one less item than 2nd line)"))
        return(list(c(TRUE, TRUE), n1))
    }
    if (n1 != n2)
        stop(wmsg("invalid file: nb of items in 2nd line is not equal to ",
                  "n2 or to n2-1 where n2 is the nb of items in 1st line"))
    item1_looks_like_a_name <- .looks_like_a_name(line1[[1L]])
    if (!item1_looks_like_a_name) {
        if (is.na(colnames))
            colnames <- any(.looks_like_a_name(line1[-1L]))
        if (is.na(rownames))
            rownames <- .looks_like_a_name(line2[[1L]])
        return(list(c(rownames, colnames), n2))
    }
    ## File contains either rownames or colnames, but not both.
    if (isTRUE(rownames) && isTRUE(colnames))
        stop(wmsg("file does not seem to contain both rownames and colnames ",
                  "(1st item in the file looks like a name)"))
    if (isFALSE(rownames) && isFALSE(colnames))
        stop(wmsg("file seems to contain either rownames or colnames ",
                  "(1st item in the file looks like a name)"))
    if (!is.na(rownames))
        return(c(rownames, !rownames))
    if (!is.na(colnames))
        return(c(!colnames, colnames))
    if (.looks_like_a_name(line2[[1L]]))
        return(c(TRUE, FALSE))
    if (any(.looks_like_a_name(line1[-1L])))
        return(c(FALSE, TRUE))
    ## File maybe is more likely to have rownames than colnames but who knows,
    ## this is just a random guess.
    c(TRUE, FALSE)
}

readSparseCSV <- function(filepath, sep=",", transpose=FALSE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string"))
    if (!isTRUEorFALSE(transpose))
        stop(wmsg("'transpose' must be TRUE or FALSE"))
    first_two_lines <- .scan_first_two_lines(filepath, sep=sep)
    line1 <- first_two_lines[[1L]]
    line2 <- first_two_lines[[2L]]
    n1 <- length(line1)
    n2 <- length(line2)
    if (n1 < 2L)
        stop(wmsg("first line in the file must contain ",
                  "at least 2 items (found ", n1, ")"))
    if (n1 != n2)
        stop(wmsg("first two lines in the file must contain ",
                  "the same number of items"))
    #dimnames_and_ncol <- .guess_dimnames_and_ncol(line1, line2,
    #                                              rownames, colnames)
    #rownames <- dimnames_and_ncol[[1L]][1L]
    #colnames <- dimnames_and_ncol[[1L]][2L]
    #ncol <- dimnames_and_ncol[[2L]]

    #filexp <- open_input_files(filepath)[[1L]]
    con <- file(filepath, "r")
    on.exit(close(con))

    C_ans <- .Call2("C_readSparseCSV", con, sep, PACKAGE="S4Arrays")
    if (transpose) {
        ans_rownames <- line1[-1L]
        ans_colnames <- C_ans[[1L]]
        ans_nzcoo1 <- C_ans[[3L]]
        ans_nzcoo2 <- C_ans[[2L]]
    } else {
        ans_rownames <- C_ans[[1L]]
        ans_colnames <- line1[-1L]
        ans_nzcoo1 <- C_ans[[2L]]
        ans_nzcoo2 <- C_ans[[3L]]
    }
    ans_nzdata <- C_ans[[4L]]
    ans_dim <- c(length(ans_rownames), length(ans_colnames))
    ans_dimnames <- list(ans_rownames, ans_colnames)

    ## Construct dgCMatrix object.
    CsparseMatrix(ans_dim, ans_nzcoo1, ans_nzcoo2, ans_nzdata,
                  dimnames=ans_dimnames)
}

readSparseTable <- function(...)
{
    .Deprecated("readSparseCSV")
    readSparseCSV(...)
}

