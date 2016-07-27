#include <R.h>
#include <Rdefines.h>
#include "MSnbase.h"
#include <stdlib.h> /* for qsort() */

/* Integer ordering code comes from S4Vectors, sorting of double
   values bases on the same principle, i.e. ensuring stable sorting.*/

static const int *aa;
static const double *dd;

/* --- .Call ENTRY POINT --- */
SEXP Integer_order(SEXP x, SEXP decreasing)
{
  int ans_length;
  SEXP ans;

  ans_length = LENGTH(x);
  PROTECT(ans = NEW_INTEGER(ans_length));
  _get_order_of_int_array(INTEGER(x), ans_length,
			  LOGICAL(decreasing)[0], INTEGER(ans), 1);
  UNPROTECT(1);
  return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP Double_order(SEXP x, SEXP decreasing)
{
  int ans_length;
  SEXP ans;

  ans_length = LENGTH(x);
  PROTECT(ans = NEW_INTEGER(ans_length));
  _get_order_of_double_array(REAL(x), ans_length,
			     LOGICAL(decreasing)[0], INTEGER(ans), 1);
  UNPROTECT(1);
  return ans;
}

/* --- .Call ENTRY POINT --- */
/* Simple test function to sort a numeric vector; is in fact
   slower than the R implementation */
SEXP sort_numeric(SEXP x, SEXP decreasing)
{
  int ans_length, i;
  SEXP ans;
  ans_length = LENGTH(x);
  int idx[ans_length];
  PROTECT(ans = NEW_NUMERIC(ans_length));
  // Use a out_shift of 0 to keep it C-like 0-base index.
  _get_order_of_double_array(REAL(x), ans_length,
			     LOGICAL(decreasing)[0], idx, 0);
  for (i = 0; i < ans_length; i++) {
    REAL(ans)[i] = REAL(x)[(idx[i])];
  }
  UNPROTECT(1);
  return ans;
}


static int compar_aa_for_stable_asc_order(const void *p1, const void *p2)
{
	int i1, i2, ret;

	i1 = *((const int *) p1);
	i2 = *((const int *) p2);
	ret = aa[i1] - aa[i2];
	if (ret != 0)
		return ret;
	/* Break tie by position so the ordering is "stable". */
	return i1 - i2;
}

/* We cannot just define compar_aa_for_stable_desc_order(p1, p2) to be
 * compar_aa_for_stable_asc_order(p2, p1) because of the tie-break
 * by position. */
static int compar_aa_for_stable_desc_order(const void *p1, const void *p2)
{
	int i1, i2, ret;

	i1 = *((const int *) p1);
	i2 = *((const int *) p2);
	ret = aa[i2] - aa[i1];
	if (ret != 0)
		return ret;
	/* Break tie by position so the ordering is "stable". */
	return i1 - i2;
}

static int compar_double_dd_for_stable_asc_order(const void *p1, const void *p2)
{
  int i1, i2, ret;

  i1 = *((const int *) p1);
  i2 = *((const int *) p2);
  if (dd[i1] > dd[i2])
    return 1;
  else if (dd[i1] < dd[i2])
    return -1;
  // return 0;
  // or i1 - i2?
  return i1 - i2;
}

static int compar_double_dd_for_stable_desc_order(const void *p1, const void *p2)
{
  int i1, i2, ret;

  i1 = *((const int *) p1);
  i2 = *((const int *) p2);
  if (dd[i1] > dd[i2])
    return -1;
  else if (dd[i1] < dd[i2])
    return 1;
  // return 0;
  // or i1 - i2?
  return i1 - i2;
}

/* Compare float/double*/
static int compar_double_asc_order(const void *p1, const void *p2)
{
  double d1, d2, ret;
  d1 = *((const double *) p1);
  d2 = *((const double *) p2);
  if (d1 < d2) return -1;
  else if (d1 > d2) return 1;
  return 0;
}

static int compar_double_desc_order(const void *p1, const void *p2)
{
  double d1, d2, ret;
  d1 = *((const double *) p1);
  d2 = *((const double *) p2);
  if (d2 < d1) return -1;
  else if (d2 > d1) return 1;
  return 0;
}


void _get_order_of_int_array(const int *x, int nelt,
		int desc, int *out, int out_shift)
{
	int i, (*compar)(const void *, const void *);

	aa = x - out_shift;
	for (i = 0; i < nelt; i++)
		out[i] = i + out_shift;
	compar = desc ? compar_aa_for_stable_desc_order :
			compar_aa_for_stable_asc_order;
	qsort(out, nelt, sizeof(int), compar);
	return;
}

void _get_order_of_double_array(const double *x, int nelt,
		int desc, int *out, int out_shift)
{
	int i, (*compar)(const void *, const void *);

	dd = x - out_shift;
	for (i = 0; i < nelt; i++)
		out[i] = i + out_shift;
	compar = desc ? compar_double_dd_for_stable_desc_order :
			compar_double_dd_for_stable_asc_order;
	qsort(out, nelt, sizeof(int), compar);
	return;
}

