/****************************************************************************
 *             Low-level manipulation of ArraySelection objects             *
 ****************************************************************************/
#include "ArraySelection_class.h"


static R_xlen_t get_SelectionTree_length_rec(SEXP selection)
{
	R_xlen_t ans;
	int list_len, k;
	SEXP list_elt;

	if (isNull(selection))
		return 0;

	if (IS_INTEGER(selection))
		return XLENGTH(selection);

	if (!isVectorList(selection))  // IS_LIST() is broken
		error("S4Arrays internal error in "
		      "get_SelectionTree_length_rec():\n"
		      "  invalid SelectionTree object");
	ans = 0;
	list_len = LENGTH(selection);
	for (k = 0; k < list_len; k++) {
		list_elt = VECTOR_ELT(selection, k);
		ans += get_SelectionTree_length_rec(list_elt);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_SelectionTree_length(SEXP x_refdim, SEXP x_selection)
{
	R_xlen_t selection_len;

	selection_len = get_SelectionTree_length_rec(x_selection);
	if (selection_len > INT_MAX)
		return ScalarReal((double) selection_len);
	return ScalarInteger((int) selection_len);
}

static int from_SelectionTree_to_matrix_rec(SEXP selection,
		int *out, int out_nrow, int out_ncol, int *out_offset,
		int *rowbuf, int rowbuf_offset)
{
	int list_len, k, ret, *p, j;
	SEXP list_elt;
	R_xlen_t nrow;

	if (isNull(selection))
		return 0;
	if (rowbuf_offset > 0) {
		if (!isVectorList(selection))  // IS_LIST() is broken
			return -1;
		list_len = LENGTH(selection);
		for (k = 0; k < list_len; k++) {
			list_elt = VECTOR_ELT(selection, k);
			rowbuf[rowbuf_offset] = k + 1;
			ret = from_SelectionTree_to_matrix_rec(list_elt,
					out, out_nrow, out_ncol, out_offset,
					rowbuf, rowbuf_offset - 1);
			if (ret < 0)
				return ret;
		}
		return 0;
	}
	if (!IS_INTEGER(selection))
		return -1;
	nrow = XLENGTH(selection);
	if (nrow > INT_MAX)
		return -1;
	for (k = 0; k < (int) nrow; k++) {
		rowbuf[0] = INTEGER(selection)[k];
		/* Copy 'rowbuf' to 'out'. */
		p = out + *out_offset;
		for (j = 0; j < out_ncol; j++) {
			*p = rowbuf[j];
			p += out_nrow;
		}
		(*out_offset)++;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SelectionTree_to_matrix(SEXP x_refdim, SEXP x_selection)
{
	R_xlen_t selection_len;
	int ans_nrow, ans_ncol, *rowbuf, out_offset, ret;
	SEXP ans;

	selection_len = get_SelectionTree_length_rec(x_selection);
	if (selection_len > INT_MAX)
		error("SelectionTree object is too long "
		      "to be turned into an ordinary matrix");

	ans_nrow = (int) selection_len;
	ans_ncol = LENGTH(x_refdim);
	rowbuf = (int *) R_alloc(ans_ncol, sizeof(int));
	ans = PROTECT(allocMatrix(INTSXP, ans_nrow, ans_ncol));

	out_offset = 0;
	ret = from_SelectionTree_to_matrix_rec(x_selection,
			INTEGER(ans), ans_nrow, ans_ncol,
			&out_offset,
			rowbuf, ans_ncol - 1);
	if (ret < 0) {
		UNPROTECT(1);
		error("S4Arrays internal error in "
		      "C_from_SelectionTree_to_matrix():\n"
		      "  invalid SelectionTree object");
	}

	/* Should never happen. */
	if (out_offset != ans_nrow) {
		UNPROTECT(1);
		error("S4Arrays internal error in "
		      "C_from_SelectionTree_to_matrix():\n"
		      "  *out_offset != ans_nrow");
	}

	UNPROTECT(1);
	return ans;
}

