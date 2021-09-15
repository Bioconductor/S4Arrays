/****************************************************************************
 *                      Combining SparseArray objects                       *
 ****************************************************************************/
#include "SparseArray_combine.h"

#include "SVT_SparseArray_class.h"


static SEXP check_and_combine_object_dims(SEXP objects, int along0,
		int *dims_along)
{
	SEXP object, dim, ans_dim;
	int nb_objects, n;

	object = VECTOR_ELT(objects, 0);
	dim = GET_SLOT(object, install("dim"));
	if (along0 < 0 || along0 >= LENGTH(dim))
		error("'along' must be >= 1 and <= the number "
		      "of dimensions of the objects to bind");
	dims_along[0] = INTEGER(dim)[along0];

	ans_dim = PROTECT(duplicate(dim));
	nb_objects = LENGTH(objects);
	for (n = 1; n < nb_objects; n++) {
		object = VECTOR_ELT(objects, n);
		dim = GET_SLOT(object, install("dim"));
		if (XLENGTH(dim) != XLENGTH(ans_dim)) {
			UNPROTECT(1);
			error("all the objects to bind must have "
			      "the same number of dimensions");
		}
		INTEGER(ans_dim)[along0] +=
			dims_along[n] = INTEGER(dim)[along0];
	}
	UNPROTECT(1);
	return ans_dim;
}

static SEXP *prepare_SVTs_buf(SEXP objects, int ndim, int along0)
{
	int nb_objects, n;
	SEXP *SVTs_buf, object;

	nb_objects = LENGTH(objects);
	SVTs_buf = (SEXP *) R_alloc(nb_objects * (ndim - along0), sizeof(SEXP));
	for (n = 0; n < nb_objects; n++) {
		object = VECTOR_ELT(objects, n);
		SVTs_buf[n] = GET_SLOT(object, install("SVT"));
	}
	return SVTs_buf;
}

static int all_NULLs(SEXP *SVTs, int nb_objects)
{
	int n;

	for (n = 0; n < nb_objects; n++)
		if (SVTs[n] != R_NilValue)
			return 0;
	return 1;
}

static SEXP concatenate_SVTs(SEXP *SVTs, int nb_objects,
		const int *dims_along, int sum_dims_along,
		SEXPTYPE ans_Rtype)
{
	int n, offset, SVT_len, k;
	SEXP SVT, ans;

	for (n = 0; n < nb_objects; n++) {
		SVT = SVTs[n];
		if (SVT != R_NilValue) {
			if (LENGTH(SVT) != dims_along[n])
				error("input object %s is an invalid "
				      "SVT_SparseArray", n + 1);
		}
	}

	ans = PROTECT(NEW_LIST(sum_dims_along));

	offset = 0;
	for (n = 0; n < nb_objects; n++) {
		SVT = SVTs[n];
		if (SVT == R_NilValue) {
			offset += dims_along[n];
			continue;
		}
		SVT_len = LENGTH(SVT);
		for (k = 0; k < SVT_len; k++, offset++)
			SET_VECTOR_ELT(ans, offset, VECTOR_ELT(SVT, k));
	}

	UNPROTECT(1);
	/* Sanity check (should never fail). */
	if (offset != sum_dims_along)
		error("S4Arrays internal error in concatenate_SVTs():\n"
		      "  offset != sum_dims_along");
	return ans;
}

/* All SVT objects in 'SVTs' are expected to have their last dimension (a.k.a.
   rightmost dimension, a.k.a. outer dimension) equal to 'd'. */
static int collect_SVTs_kth_elt(SEXP *SVTs, int nb_objects, int k, int d,
		SEXP *subSVTs_buf)
{
	int n;
	SEXP SVT, subSVT;

	for (n = 0; n < nb_objects; n++) {
		SVT = SVTs[n];
		/* 'SVT' must be NULL or a list of length 'd'. */
		if (SVT == R_NilValue) {
			subSVT = R_NilValue;
		} else {
			if (!isVectorList(SVT))  // IS_LIST() is broken
				return -1;
			if (LENGTH(SVT) != d)
				return -1;
			subSVT = VECTOR_ELT(SVT, k);
		}
		subSVTs_buf[n] = subSVT;
	}
	return 0;
}

/* Recursive.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_combine_SVTs(SEXP *SVTs, int nb_objects,
		const int *ans_dim, int ndim, int along0,
		const int *dims_along, SEXPTYPE ans_Rtype)
{
	SEXP *subSVTs_buf, ans, ans_elt;
	int ans_len, is_empty, k, ret;

	if (all_NULLs(SVTs, nb_objects))
		return R_NilValue;

	if (along0 == 0)
		return R_NilValue; /* not ready yet */

	if (along0 == ndim - 1)
		return concatenate_SVTs(SVTs, nb_objects,
					dims_along, ans_dim[along0],
					ans_Rtype);

	subSVTs_buf = SVTs + nb_objects;
	ans_len = ans_dim[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (k = 0; k < ans_len; k++) {
		ret = collect_SVTs_kth_elt(SVTs, nb_objects, k, ans_len,
					   subSVTs_buf);
		if (ret < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in "
			      "REC_combine_SVTs():\n"
			      "  collect_SVTs_kth_elt() returned an error");
		}
		ans_elt = REC_combine_SVTs(subSVTs_buf, nb_objects,
					   ans_dim, ndim - 1, along0,
					   dims_along, ans_Rtype);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, k, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_abind_SVT_SparseArray_objects(SEXP objects, SEXP along, SEXP ans_type)
{
	SEXPTYPE ans_Rtype;
	int along0, nb_objects, *dims_along, ndim;
	SEXP ans_dim, *SVTs_buf, ans_SVT, ans;

	if (!isVectorList(objects))  // IS_LIST() is broken
		error("'objects' must be a list of SVT_SparseArray objects");

	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	if (!IS_INTEGER(along) || XLENGTH(along) != 1)
		error("'along' must be a single positive integer");
	along0 = INTEGER(along)[0] - 1;

	nb_objects = LENGTH(objects);
	if (nb_objects == 0)
		error("'objects' cannot be an empty list");

	dims_along = (int *) R_alloc(nb_objects, sizeof(int));
	ans_dim = PROTECT(
		check_and_combine_object_dims(objects, along0, dims_along)
	);
	ndim = LENGTH(ans_dim);

	SVTs_buf = prepare_SVTs_buf(objects, ndim, along0);
	ans_SVT = REC_combine_SVTs(SVTs_buf, nb_objects,
				   INTEGER(ans_dim), ndim, along0,
				   dims_along, ans_Rtype);
	if (ans_SVT != R_NilValue)
		PROTECT(ans_SVT);

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	if (ans_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 1, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return ans;
}

