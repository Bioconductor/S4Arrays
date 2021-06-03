/****************************************************************************
 *                    Workhorse behind readSparseTable()                    *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "read_sparse_table.h"

#include "S4Vectors_interface.h"
#include "XVector_interface.h"

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
 * C_read_sparse_table()
 */

static char errmsg_buf[200];

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

static void load_rowname(const char *data, int data_len,
			 CharAEAE *rownames_buf)
{
	CharAE *ae = new_CharAE(data_len);
	memcpy(ae->elts, data, data_len);
	/* We don't use CharAE_set_nelt() for maximum speed. */
	ae->_nelt = data_len;
	CharAEAE_fast_append(rownames_buf, ae);
}

static void load_table_data(const char *data, int data_len,
		int row_idx, int col_idx,
		IntAE *nzindex1_buf, IntAE *nzindex2_buf, IntAE *nzdata_buf)
{
	int val;

	data_len = delete_trailing_LF_or_CRLF(data, data_len);
	val = as_int(data, data_len);
	if (val != 0) {
		IntAE_fast_append(nzindex1_buf, row_idx);
		IntAE_fast_append(nzindex2_buf, col_idx);
		IntAE_fast_append(nzdata_buf, val);
	}
	return;
}

static void load_table_row(const char *line, char sep, int row_idx,
		CharAEAE *rownames_buf,
		IntAE *nzindex1_buf, IntAE *nzindex2_buf, IntAE *nzdata_buf)
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
			load_rowname(data, data_len, rownames_buf);
		} else {
			load_table_data(data, data_len,
					row_idx, col_idx,
					nzindex1_buf, nzindex2_buf,
					nzdata_buf);
		}
		col_idx++;
		data = line + i;
		data_len = 0;
	}
	load_table_data(data, data_len,
			row_idx, col_idx,
			nzindex1_buf, nzindex2_buf,
			nzdata_buf);
	return;
}

static const char *read_sparse_table(SEXP filexp, char sep,
		CharAEAE *rownames_buf,
		IntAE *nzindex1_buf, IntAE *nzindex2_buf,
		IntAE *nzdata_buf)
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
		load_table_row(buf, sep, row_idx,
			       rownames_buf,
			       nzindex1_buf, nzindex2_buf, nzdata_buf);
		row_idx++;
	}
	return NULL;
}

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

/* Args:
 *   filexp: A "file external pointer" (see src/io_utils.c in the XVector
 *           package). TODO: Support connections (see src/readGFF.c in the
 *           rtracklayer package for how to do that).
 * Return 'list(rownames, nzindex1, nzindex2, nzdata)'.
 */
SEXP C_read_sparse_table(SEXP filexp, SEXP sep)
{
	CharAEAE *rownames_buf;
	IntAE *nzindex1_buf;
	IntAE *nzindex2_buf;
	IntAE *nzdata_buf;
	const char *errmsg;
	SEXP ans, ans_elt;

	rownames_buf = new_CharAEAE(0, 0);
	nzindex1_buf = new_IntAE(0, 0, 0);
	nzindex2_buf = new_IntAE(0, 0, 0);
	nzdata_buf = new_IntAE(0, 0, 0);

	errmsg = read_sparse_table(filexp, get_sep_char(sep),
				   rownames_buf,
				   nzindex1_buf, nzindex2_buf,
				   nzdata_buf);
	if (errmsg != NULL)
		error("reading file: %s", errmsg);

	ans = PROTECT(NEW_LIST(4));

	ans_elt = PROTECT(new_CHARACTER_from_CharAEAE(rownames_buf));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(new_INTEGER_from_IntAE(nzindex1_buf));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(new_INTEGER_from_IntAE(nzindex2_buf));
	SET_VECTOR_ELT(ans, 2, ans_elt);
	UNPROTECT(1);

	ans_elt = PROTECT(new_INTEGER_from_IntAE(nzdata_buf));
	SET_VECTOR_ELT(ans, 3, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

