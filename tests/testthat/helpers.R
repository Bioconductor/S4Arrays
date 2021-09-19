check_SparseArray_object <- function(object, expected_class, a0)
{
    expect_true(class(object) == expected_class)
    expect_true(validObject(object))
    expect_identical(dim(object), dim(a0))
    expect_identical(dimnames(object), dimnames(a0))
    expect_identical(type(object), type(a0))
    a <- as.array(object)
    expect_identical(a, a0)
}

