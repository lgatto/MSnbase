#include <R.h>
#include <Rinternals.h>
#include <string.h>

/* SEXP Spectrum1_constructor( */
/* 			   SEXP msLevel, */
/* 			   SEXP peaksCount, */
/* 			   SEXP rt, */
/* 			   SEXP acquisitionNum, */
/* 			   SEXP scanIndex, */
/* 			   SEXP tic, */
/* 			   SEXP mz, */
/* 			   SEXP intensity, */
/* 			   SEXP fromFile, */
/* 			   SEXP centroided, */
/* 			   SEXP smoothed, */
/* 			   SEXP polarity, */
/* 			   SEXP check */
/* 			   ); */

/* SEXP _logical_Rle_constructor( */
/* 			      const int *values, */
/* 			      int nvalues, */
/* 			      const int *lengths, */
/* 			      int buflength */
/* ); */

/* SEXP Rle_constructor( */
/* 		     SEXP values, */
/* 		     SEXP lengths, */
/* 		     SEXP check, */
/* 		     SEXP buflength */
/* ); */

void _get_order_of_int_array(
	const int *x,
	int nelt,
	int desc,
	int *out,
	int out_shift
);

void _get_order_of_double_array(
	const double *x,
	int nelt,
	int desc,
	int *out,
	int out_shift
);
