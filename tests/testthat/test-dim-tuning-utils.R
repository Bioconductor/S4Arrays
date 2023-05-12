
test_that("tune_dims()", {
    tune_dims <- S4Arrays:::tune_dims
    dim <- c(A=1L, B=4L, C=15L)

    ## No-op.
    dim_tuner <- c(0L, 0L, 0L)
    expect_identical(tune_dims(dim, dim_tuner), dim)

    ## Drop ineffective dimensions.
    dim_tuner <- c(-1L, 0L, 0L)
    current <- tune_dims(dim, dim_tuner)
    expect_identical(current, dim[-1L])
    expect_identical(tune_dims(current, -dim_tuner),
                     setNames(dim, c("", "B", "C")))

    ## A simple "pure R" version of tune_dims() that doesn't validate
    ## the 'dim_tuner' vector.
    simple_tune_dims <- function(dim, dim_tuner) {
        is_keep_or_add <- dim_tuner >= 0L
        is_keep_or_drop <- dim_tuner <= 0L
        idx <- rep.int(NA_integer_, sum(is_keep_or_add))
        i1 <- which(dim_tuner[is_keep_or_add] == 0L)
        i2 <- which(dim_tuner[is_keep_or_drop] == 0L)
        idx[i1] <- i2
	ans <- dim[idx]
        na_idx <- which(is.na(ans))
        ans[na_idx] <- 1L
        names(ans)[na_idx] <- ""
        ans
    }

    ## Add ineffective dimensions.
    dim_tuner <- c(1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L)
    current <- tune_dims(dim, dim_tuner)
    expected <- simple_tune_dims(dim, dim_tuner)
    expect_identical(current, expected)
    expect_identical(tune_dims(current, -dim_tuner), dim)

    ## Add and drop ineffective dimensions.
    dim_tuner <- c(-1L, 0L, 1L, 0L, 1L, 1L)
    current <- tune_dims(dim, dim_tuner)
    expected <- simple_tune_dims(dim, dim_tuner)
    expect_identical(current, expected)
    expect_identical(tune_dims(current, -dim_tuner),
                     setNames(dim, c("", "B", "C")))

    ## Trying to drop effective dimensions.
    expect_error(tune_dims(dim, c(0L, -1L, 0L)), "internal error")
    expect_error(tune_dims(dim, c(1L, 0L, -1L, 0L)), "internal error")

    ## Invalid 'dim_tuner' vector.
    dim <- c(1L, 1L, 1L)
    expect_error(tune_dims(dim, c(0L, 0L, 2L)), "internal error")
    expect_error(tune_dims(dim, c(0L, 0L)), "internal error")
    expect_error(tune_dims(dim, c(-1L, 0L)), "internal error")
    expect_error(tune_dims(dim, c(-1L, 0L, 1L)), "internal error")
    expect_error(tune_dims(dim, c(0L, 0L, 0L, 0L)), "internal error")
    expect_error(tune_dims(dim, c(0L, -1L, 0L, 0L, 1L)), "internal error")
    expect_error(tune_dims(dim, c(-1L, -1L, -1L)), "internal error")
})

test_that("tune_dimnames()", {
    tune_dimnames <- S4Arrays:::tune_dimnames
    dimnames <- list(NULL, "B", "C")

    expect_identical(tune_dimnames(dimnames, c(0L, 0L, 0L)), dimnames)
    expect_identical(tune_dimnames(dimnames, c(0L, 0L, -1L)), dimnames[-3L])
    expect_identical(tune_dimnames(dimnames, c(0L, -1L, -1L)), NULL)
    expect_identical(tune_dimnames(dimnames, c(1L, 0L, -1L, 0L, 1L)),
                     list(NULL, NULL, "C", NULL))
    expect_identical(tune_dimnames(dimnames, c(1L, 0L, -1L, -1L)), NULL)

    expect_error(tune_dimnames(dimnames, c(0L, 0L, -1L, 0L)), "internal error")
})

