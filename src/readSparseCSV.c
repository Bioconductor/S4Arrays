/****************************************************************************
 *                     Workhorse behind readSparseCSV()                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "readSparseCSV.h"

#include "S4Vectors_interface.h"
#include "XVector_interface.h"

#include "leaf_vector_utils.h"

#include <R_ext/Connections.h>  /* for R_ReadConnection() */

#include <string.h>  /* for memcpy */

#define	IOBUF_SIZE 8000002


/****************************************************************************
 * filexp_gets2(): A version of filexp_gets() that also works on connections
 * Copied from rtracklayer/src/readGFF.c
 */

Rconnection getConnection(int n);  /* not in <R_ext/Connections.h>, why? */

static char con_buf[250000];
static int con_buf_len, con_buf_offset;

static void init_con_buf()
{
	con_buf_len = con_buf_offset = 0;
	return;
}

static int filexp_gets2(SEXP filexp, char *buf, int buf_size, int *EOL_in_buf)
{
	Rconnection con;
	int buf_offset;
	char c;

	if (TYPEOF(filexp) == EXTPTRSXP)
		return filexp_gets(filexp, buf, buf_size, EOL_in_buf);
	buf_offset = *EOL_in_buf = 0;
	while (buf_offset < buf_size - 1) {
		if (con_buf_offset == con_buf_len) {
			con = getConnection(asInteger(filexp));
			con_buf_len = (int) R_ReadConnection(con,
					con_buf,
					sizeof(con_buf) / sizeof(char));
			if (con_buf_len == 0)
				break;
			con_buf_offset = 0;
		}
		c = con_buf[con_buf_offset++];
		buf[buf_offset++] = c;
		if (c == '\n') {
			*EOL_in_buf = 1;
			break;
		}
	}
	buf[buf_offset] = '\0';
	if (buf_offset == 0)
		return 0;
	if (con_buf_len == 0 || *EOL_in_buf)
		return 2;
	return 1;
}


/****************************************************************************
 * Stuff shared by C_readSparseCSV_as_SVT_SparseMatrix() and
 * C_readSparseCSV_as_COO_SparseMatrix() below
 */

static char errmsg_buf[200];

static char get_sep_char(SEXP sep)
{
	SEXP sep0;

	if (!(IS_CHARACTER(sep) && LENGTH(sep) == 1))
		error("'sep' must be a single character");
	sep0 = STRING_ELT(sep, 0);
	if (sep0 == NA_STRING || length(sep0) != 1)
		error("'sep' must be a single character");
	return CHAR(sep0)[0];
}

static inline void IntAE_fast_append(IntAE *ae, int val)
{
	/* We don't use IntAE_get_nelt() for maximum speed. */
	if (ae->_nelt == ae->_buflength)
		IntAE_extend(ae, increase_buflength(ae->_buflength));
	ae->elts[ae->_nelt++] = val;
	return;
}

static inline void CharAEAE_fast_append(CharAEAE *aeae, CharAE *ae)
{
	/* We don't use CharAEAE_get_nelt() for maximum speed. */
	if (aeae->_nelt == aeae->_buflength)
		CharAEAE_extend(aeae, increase_buflength(aeae->_buflength));
	aeae->elts[aeae->_nelt++] = ae;
	return;
}

static void load_csv_rowname(const char *data, int data_len,
			     CharAEAE *csv_rownames_buf)
{
	CharAE *ae = new_CharAE(data_len);
	memcpy(ae->elts, data, data_len);
	/* We don't use CharAE_set_nelt() for maximum speed. */
	ae->_nelt = data_len;
	CharAEAE_fast_append(csv_rownames_buf, ae);
}


/****************************************************************************
 * C_readSparseCSV_as_SVT_SparseMatrix()
 */

static void load_csv_data1(const char *data, int data_len,
		int off, IntAE *offs_buf, IntAE *vals_buf)
{
	int val;

	data_len = delete_trailing_LF_or_CRLF(data, data_len);
	if (data_len == 0 || (val = as_int(data, data_len)) == 0)
		return;
	IntAE_fast_append(offs_buf, off);
	IntAE_fast_append(vals_buf, val);
	return;
}

static void load_csv_row_to_SVT_bufs(const char *line, char sep,
		CharAEAE *csv_rownames_buf, IntAE *offs_buf, IntAE *vals_buf)
{
	int col_idx, i, data_len;
	const char *data;
	char c;

	col_idx = i = 0;
	data = line;
	data_len = 0;
	while ((c = line[i++])) {
		if (c != sep) {
			data_len++;
			continue;
		}
		if (col_idx == 0) {
			load_csv_rowname(data, data_len, csv_rownames_buf);
		} else {
			load_csv_data1(data, data_len,
				       col_idx - 1, offs_buf, vals_buf);
		}
		col_idx++;
		data = line + i;
		data_len = 0;
	}
	load_csv_data1(data, data_len,
		       col_idx - 1, offs_buf, vals_buf);
	return;
}

static const char *read_sparse_csv_as_SVT_SparseMatrix(
		SEXP filexp, char sep,
		CharAEAE *csv_rownames_buf,
		IntAEAE *offs_bufs, IntAEAE *vals_bufs)
{
	int row_idx, lineno, ret_code, EOL_in_buf;
	char buf[IOBUF_SIZE];
	IntAE *offs_buf, *vals_buf;

	if (TYPEOF(filexp) != EXTPTRSXP)
		init_con_buf();
	row_idx = 1;
	for (lineno = 1;
	     (ret_code = filexp_gets2(filexp, buf, IOBUF_SIZE, &EOL_in_buf));
	     lineno += EOL_in_buf)
	{
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (!EOL_in_buf) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, "
				 "line is too long", lineno);
			return errmsg_buf;
		}
		if (lineno == 1)
			continue;
		offs_buf = new_IntAE(0, 0, 0);
		vals_buf = new_IntAE(0, 0, 0);
		load_csv_row_to_SVT_bufs(buf, sep,
					 csv_rownames_buf, offs_buf, vals_buf);
		IntAEAE_insert_at(offs_bufs, row_idx - 1, offs_buf);
		IntAEAE_insert_at(vals_bufs, row_idx - 1, vals_buf);
		row_idx++;
	}
	return NULL;
}

/* 'offs_buf' and 'vals_buf' are **asumed** to have the same length and this
   length is **assumed** to be != 0. We don't check this! */
static SEXP make_leaf_vector_from_bufs(const IntAE *offs_buf,
				       const IntAE *vals_buf)
{
	SEXP ans_offs, ans_vals, ans;

	ans_offs = PROTECT(new_INTEGER_from_IntAE(offs_buf));
	ans_vals = PROTECT(new_INTEGER_from_IntAE(vals_buf));
	ans = _new_leaf_vector(ans_offs, ans_vals);  // unprotected!
	UNPROTECT(2);
	return ans;
}

/* 'offs_bufs' and 'vals_bufs' are **asumed** to have the same shape.
   We don't check this! */
static SEXP make_SVT_from_bufs(const IntAEAE *offs_bufs,
			       const IntAEAE *vals_bufs)
{
	int SVT_len, is_empty, i, lv_len;
	SEXP ans, ans_elt;
	const IntAE *offs_buf, *vals_buf;

	SVT_len = IntAEAE_get_nelt(offs_bufs);
	ans = PROTECT(NEW_LIST(SVT_len));
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		offs_buf = offs_bufs->elts[i];
		lv_len = IntAE_get_nelt(offs_buf);
		if (lv_len == 0)
			continue;
		vals_buf = vals_bufs->elts[i];
		ans_elt = make_leaf_vector_from_bufs(offs_buf, vals_buf);
		PROTECT(ans_elt);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
		is_empty = 0;
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   filexp: A "file external pointer" (see src/io_utils.c in the XVector
 *           package). TODO: Support connections (see src/readGFF.c in the
 *           rtracklayer package for how to do that).
 * Returns 'list(csv_rownames, SVT)'.
 */
SEXP C_readSparseCSV_as_SVT_SparseMatrix(SEXP filexp, SEXP sep)
{
	CharAEAE *csv_rownames_buf;
	IntAEAE *offs_bufs, *vals_bufs;
	const char *errmsg;
	SEXP ans, ans_elt;

	csv_rownames_buf = new_CharAEAE(0, 0);
	offs_bufs = new_IntAEAE(0, 0);
	vals_bufs = new_IntAEAE(0, 0);

	errmsg = read_sparse_csv_as_SVT_SparseMatrix(
				filexp, get_sep_char(sep),
				csv_rownames_buf, offs_bufs, vals_bufs);
	if (errmsg != NULL)
		error("reading file: %s", errmsg);

	ans = PROTECT(NEW_LIST(2));

	ans_elt = PROTECT(new_CHARACTER_from_CharAEAE(csv_rownames_buf));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(make_SVT_from_bufs(offs_bufs, vals_bufs));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_readSparseCSV_as_COO_SparseMatrix()
 */

static void load_csv_data2(const char *data, int data_len,
		int row_idx, int col_idx,
		IntAE *nzcoo1_buf, IntAE *nzcoo2_buf, IntAE *nzvals_buf)
{
	int val;

	data_len = delete_trailing_LF_or_CRLF(data, data_len);
	if (data_len == 0 || (val = as_int(data, data_len)) == 0)
		return;
	IntAE_fast_append(nzcoo1_buf, row_idx);
	IntAE_fast_append(nzcoo2_buf, col_idx);
	IntAE_fast_append(nzvals_buf, val);
	return;
}

static void load_csv_row_to_COO_bufs(const char *line, char sep, int row_idx,
		CharAEAE *csv_rownames_buf,
		IntAE *nzcoo1_buf, IntAE *nzcoo2_buf, IntAE *nzvals_buf)
{
	int col_idx, i, data_len;
	const char *data;
	char c;

	col_idx = i = 0;
	data = line;
	data_len = 0;
	while ((c = line[i++])) {
		if (c != sep) {
			data_len++;
			continue;
		}
		if (col_idx == 0) {
			load_csv_rowname(data, data_len, csv_rownames_buf);
		} else {
			load_csv_data2(data, data_len,
				       row_idx, col_idx,
				       nzcoo1_buf, nzcoo2_buf,
				       nzvals_buf);
		}
		col_idx++;
		data = line + i;
		data_len = 0;
	}
	load_csv_data2(data, data_len,
		       row_idx, col_idx,
		       nzcoo1_buf, nzcoo2_buf,
		       nzvals_buf);
	return;
}

static const char *read_sparse_csv_as_COO_SparseMatrix(
		SEXP filexp, char sep,
		CharAEAE *csv_rownames_buf,
		IntAE *nzcoo1_buf, IntAE *nzcoo2_buf,
		IntAE *nzvals_buf)
{
	int row_idx, lineno, ret_code, EOL_in_buf;
	char buf[IOBUF_SIZE];

	if (TYPEOF(filexp) != EXTPTRSXP)
		init_con_buf();
	row_idx = 1;
	for (lineno = 1;
	     (ret_code = filexp_gets2(filexp, buf, IOBUF_SIZE, &EOL_in_buf));
	     lineno += EOL_in_buf)
	{
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (!EOL_in_buf) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, "
				 "line is too long", lineno);
			return errmsg_buf;
		}
		if (lineno == 1)
			continue;
		load_csv_row_to_COO_bufs(buf, sep, row_idx,
				csv_rownames_buf,
				nzcoo1_buf, nzcoo2_buf, nzvals_buf);
		row_idx++;
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   filexp: A "file external pointer" (see src/io_utils.c in the XVector
 *           package). TODO: Support connections (see src/readGFF.c in the
 *           rtracklayer package for how to do that).
 * Returns 'list(csv_rownames, nzcoo1, nzcoo2, nzvals)'.
 */
SEXP C_readSparseCSV_as_COO_SparseMatrix(SEXP filexp, SEXP sep)
{
	CharAEAE *csv_rownames_buf;
	IntAE *nzcoo1_buf, *nzcoo2_buf, *nzvals_buf;
	const char *errmsg;
	SEXP ans, ans_elt;

	csv_rownames_buf = new_CharAEAE(0, 0);
	nzcoo1_buf = new_IntAE(0, 0, 0);
	nzcoo2_buf = new_IntAE(0, 0, 0);
	nzvals_buf = new_IntAE(0, 0, 0);

	errmsg = read_sparse_csv_as_COO_SparseMatrix(
				filexp, get_sep_char(sep),
				csv_rownames_buf,
				nzcoo1_buf, nzcoo2_buf,
				nzvals_buf);
	if (errmsg != NULL)
		error("reading file: %s", errmsg);

	ans = PROTECT(NEW_LIST(4));

	ans_elt = PROTECT(new_CHARACTER_from_CharAEAE(csv_rownames_buf));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(new_INTEGER_from_IntAE(nzcoo1_buf));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(new_INTEGER_from_IntAE(nzcoo2_buf));
	SET_VECTOR_ELT(ans, 2, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(new_INTEGER_from_IntAE(nzvals_buf));
	SET_VECTOR_ELT(ans, 3, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

