/****************************************************************************
 *             Low-level manipulation of ArraySelection objects             *
 ****************************************************************************/
#include "ArraySelection_class.h"


/* Recursive. */
static R_xlen_t count_tree_leaves_REC(SEXP tree)
{
	R_xlen_t ans;
	int tree_len, k;
	SEXP subtree;

	if (isNull(tree))
		return 0;

	if (IS_INTEGER(tree))
		return XLENGTH(tree);

	if (!isVectorList(tree))  // IS_LIST() is broken
		error("S4Arrays internal error in "
		      "count_tree_leaves_REC():\n"
		      "  invalid SelectionTree object");
	ans = 0;
	tree_len = LENGTH(tree);
	for (k = 0; k < tree_len; k++) {
		subtree = VECTOR_ELT(tree, k);
		ans += count_tree_leaves_REC(subtree);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_SelectionTree_length(SEXP x_refdim, SEXP x_selection)
{
	R_xlen_t selection_len;

	selection_len = count_tree_leaves_REC(x_selection);
	if (selection_len > INT_MAX)
		return ScalarReal((double) selection_len);
	return ScalarInteger((int) selection_len);
}


/****************************************************************************
 * Going back and forth between SelectionTree and matrix objects
 */

/* Recursive. */
static int dump_tree_as_matrix_REC(SEXP tree,
		int *out, int out_nrow, int out_ncol, int *out_offset,
		int *rowbuf, int rowbuf_offset)
{
	int tree_len, k, ret, *p, j;
	SEXP subtree;
	R_xlen_t nrow;

	if (isNull(tree))
		return 0;
	if (rowbuf_offset > 0) {
		if (!isVectorList(tree))  // IS_LIST() is broken
			return -1;
		tree_len = LENGTH(tree);
		for (k = 0; k < tree_len; k++) {
			subtree = VECTOR_ELT(tree, k);
			rowbuf[rowbuf_offset] = k + 1;
			ret = dump_tree_as_matrix_REC(subtree,
					out, out_nrow, out_ncol, out_offset,
					rowbuf, rowbuf_offset - 1);
			if (ret < 0)
				return ret;
		}
		return 0;
	}
	if (!IS_INTEGER(tree))
		return -1;
	nrow = XLENGTH(tree);
	if (nrow > INT_MAX)
		return -1;
	for (k = 0; k < (int) nrow; k++) {
		rowbuf[0] = INTEGER(tree)[k];
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

	selection_len = count_tree_leaves_REC(x_selection);
	if (selection_len > INT_MAX)
		error("SelectionTree object is too long "
		      "to be turned into an ordinary matrix");

	ans_nrow = (int) selection_len;
	ans_ncol = LENGTH(x_refdim);
	rowbuf = (int *) R_alloc(ans_ncol, sizeof(int));
	ans = PROTECT(allocMatrix(INTSXP, ans_nrow, ans_ncol));

	out_offset = 0;
	ret = dump_tree_as_matrix_REC(x_selection,
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

static SEXP make_depth1_tree(const int *in, int in_len, int refdim1)
{
	SEXP ans;
	int i, v;

	ans = PROTECT(NEW_INTEGER(in_len));
	for (i = 0; i < in_len; i++) {
		v = in[i];
		if (v < 1 || v > refdim1) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
		INTEGER(ans)[i] = v;
	}
	UNPROTECT(1);
	return ans;
}

static int grow_tree(const int *in,
		int in_nrow, int in_ncol, int in_offset,
		SEXP tree, const int *refdim)
{
	const int *p;
	int j, k;
	SEXP subtree;

	p = in + in_offset;
	if (*p < 1  || *p > refdim[0])
		return -1;

	if (in_ncol >= 3) {
		p += in_nrow * in_ncol;
		for (j = in_ncol - 2; j >= 1; j--) {
			p -= in_nrow;
			k = *p - 1;
			if (k < 0 || k >= LENGTH(tree))
				return -1;
			subtree = VECTOR_ELT(tree, k);
			if (j == 1)
				break;
			/* 'subtree' is NULL or a list. */
			if (subtree == R_NilValue) {
				subtree = PROTECT(NEW_LIST(refdim[j]));
				SET_VECTOR_ELT(tree, k, subtree);
				UNPROTECT(1);
			}
			tree = subtree;
		}
		/* 'subtree' is NULL or an integer vector of counts. */
		if (subtree == R_NilValue) {
			subtree = PROTECT(NEW_INTEGER(refdim[j]));
			SET_VECTOR_ELT(tree, k, subtree);
			UNPROTECT(1);
		}
		tree = subtree;
	}

	p = in + in_offset + in_nrow;
	k = *p - 1;
	if (k < 0 || k >= LENGTH(tree))
		return -1;
	INTEGER(tree)[k]++;
	return 0;
}

static SEXP make_terminal_list(SEXP x)
{
	int x_len, k, ans_elt_len;
	SEXP ans, ans_elt;

	x_len = LENGTH(x);
	ans = PROTECT(NEW_LIST(x_len));
	for (k = 0; k < x_len; k++) {
		ans_elt_len = INTEGER(x)[k];
		if (ans_elt_len != 0) {
			ans_elt = PROTECT(NEW_INTEGER(ans_elt_len));
			SET_VECTOR_ELT(ans, k, ans_elt);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}

static int add_matrix_row_to_tree(const int *in,
		int in_nrow, int in_ncol, int in_offset,
		SEXP tree)
{
	const int *p;
	int j, k;
	SEXP subtree;

	if (in_ncol >= 3) {
		p = in + in_offset + in_nrow * in_ncol;
		for (j = in_ncol - 2; j >= 1; j--) {
			p -= in_nrow;
			k = *p - 1;
			subtree = VECTOR_ELT(tree, k);
			if (j == 1)
				break;
			tree = subtree;
		}
		/* 'subtree' is an integer vector of counts or a list. */
		if (IS_INTEGER(subtree)) {
			subtree = PROTECT(make_terminal_list(subtree));
			SET_VECTOR_ELT(tree, k, subtree);
			UNPROTECT(1);
		}
		tree = subtree;
	}

	p = in + in_offset + in_nrow;
	k = *p - 1;
	subtree = VECTOR_ELT(tree, k);
	for (k = 0; k < LENGTH(subtree); k++) {
		if (INTEGER(subtree)[k] == 0)
			break;
	}
	p = in + in_offset;
	INTEGER(subtree)[k] = *p;
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_matrix_to_SelectionTree(SEXP m, SEXP refdim)
{
	SEXP m_dim, ans;
	int m_nrow, m_ncol, ans_len, i, ret;

	m_dim = GET_DIM(m);
	if (LENGTH(m_dim) != 2)
		error("'m' must be a matrix");

	m_nrow = INTEGER(m_dim)[0];
	m_ncol = INTEGER(m_dim)[1];
	if (m_ncol != LENGTH(refdim))
		error("ncol(m) != length(refdim)");

	if (m_nrow == 0)
		return R_NilValue;

	if (m_ncol == 1)
		return make_depth1_tree(INTEGER(m), m_nrow, INTEGER(refdim)[0]);

	ans_len = INTEGER(refdim)[m_ncol - 1];

	/* 1st pass: Grow the branches of the tree but don't add any
	   leaves to it, only count them. */
	if (m_ncol == 2) {
		ans = PROTECT(NEW_INTEGER(ans_len));
	} else {
		ans = PROTECT(NEW_LIST(ans_len));
	}
	for (i = 0; i < m_nrow; i++) {
		ret = grow_tree(INTEGER(m),
				m_nrow, m_ncol, i,
				ans, INTEGER(refdim));
		if (ret < 0) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
	}

	/* 2nd pass: Add the leaves to the tree. */
	if (m_ncol == 2)
		ans = PROTECT(make_terminal_list(ans));
	for (i = 0; i < m_nrow; i++) {
		ret = add_matrix_row_to_tree(INTEGER(m),
					     m_nrow, m_ncol, i,
					     ans);
		if (ret < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in "
			      "C_from_matrix_to_SelectionTree():\n"
			      "  add_matrix_row_to_tree() "
			      "returned an error");
		}
	}

	if (m_ncol == 2)
		UNPROTECT(1);
	UNPROTECT(1);
	return ans;
}

