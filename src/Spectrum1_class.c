#include <R.h>
#include <Rdefines.h>
#include "Spectrum1_class.h"

/************************************************************
 * The Spectrum1 constructor in C.
 *
 */

static SEXP _new_Spectrum1(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP polarity)
{
  SEXP classdef, ans;

  PROTECT(classdef = MAKE_CLASS("Spectrum1"));
  PROTECT(ans = NEW_OBJECT(classdef));
  /* Set the slot values */
  /* SET_SLOT(ans, install("msLevel", msLevel)); */
  /* SET_SLOT(ans, install("peaksCount", peaksCount)); */
  /* SET_SLOT(ans, install("rt", rt)); */
  /* SET_SLOT(ans, install("acquisitionNum", acquisitionNum)); */
  /* SET_SLOT(ans, install("scanIndex", scanIndex)); */
  /* SET_SLOT(ans, install("tic", tic)); */
  /* SET_SLOT(ans, install("mz", mz)); */
  /* SET_SLOT(ans, install("intensity", intensity)); */
  /* SET_SLOT(ans, install("fromFile", fromFile)); */
  /* SET_SLOT(ans, install("centroided", centroided)); */
  /* SET_SLOT(ans, install("polarity", polarity)); */
  SET_SLOT(ans, install("msLevel"), msLevel);
  SET_SLOT(ans, install("peaksCount"), peaksCount);
  SET_SLOT(ans, install("rt"), rt);
  SET_SLOT(ans, install("acquisitionNum"), acquisitionNum);
  SET_SLOT(ans, install("scanIndex"), scanIndex);
  SET_SLOT(ans, install("tic"), tic);
  SET_SLOT(ans, install("mz"), mz);
  SET_SLOT(ans, install("intensity"), intensity);
  SET_SLOT(ans, install("fromFile"), fromFile);
  SET_SLOT(ans, install("centroided"), centroided);
  SET_SLOT(ans, install("polarity"), polarity);

  UNPROTECT(2);
  return ans;
}

SEXP Spectrum1_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP polarity, SEXP check)
{
  //int nvalues;
  SEXP ans;
  //const int *lengths_p;

  //nvalues = LENGTH(intensity);
  if(LOGICAL(check)[0]) {
    if( LENGTH(mz) != LENGTH(intensity)) {
      error("'length(intensity)' != 'length(mz)'");
    }
  }

  PROTECT(ans = _new_Spectrum1(msLevel, peaksCount, rt, acquisitionNum,
			      scanIndex, tic, mz, intensity, fromFile,
			      centroided, polarity));
  UNPROTECT(1);
  return(ans);
  return R_NilValue;
}


/* static SEXP _new_Rle(SEXP values, SEXP lengths) */
/* { */
/*   SEXP classdef, ans; */

/*   PROTECT(classdef = MAKE_CLASS("Rle")); */
/*   PROTECT(ans = NEW_OBJECT(classdef)); */
/*   SET_SLOT(ans, install("values"), values); */
/*   SET_SLOT(ans, install("lengths"), lengths); */
/*   UNPROTECT(2); */
/*   return ans; */
/* } */

/* SEXP _logical_Rle_constructor(const int *values, int nvalues, */
/* 			      const int *lengths, int buflength) */
/* { */
/*   int nrun, *buf_lengths; */
/*   int *buf_values; */
/*   SEXP ans_lengths, ans_values, ans; */

/*   if (buflength > nvalues) */
/*     buflength = nvalues; */
/*   if (buflength == 0) { */
/*     /\* 1st pass: compute only the nb of runs *\/ */
/*     buf_values = NULL; */
/*     buf_lengths = NULL; */
/*   } else { */
/*     buf_values = (int *) R_alloc(buflength, sizeof(int)); */
/*     buf_lengths = (int *) R_alloc(buflength, sizeof(int)); */
/*   } */
/*   nrun = compute_int_runs(values, nvalues, lengths, */
/* 			  buf_values, buf_lengths); */
/*   PROTECT(ans_values = NEW_LOGICAL(nrun)); */
/*   PROTECT(ans_lengths = NEW_INTEGER(nrun)); */
/*   if (buflength == 0) { */
/*     /\* 2nd pass: fill 'ans_values' and 'ans_lengths' *\/ */
/*     compute_int_runs(values, nvalues, lengths, */
/* 		     LOGICAL(ans_values), INTEGER(ans_lengths)); */
/*   } else { */
/*     memcpy(LOGICAL(ans_values), buf_values, nrun * sizeof(int)); */
/*     memcpy(INTEGER(ans_lengths), buf_lengths, nrun * sizeof(int)); */
/*   } */
/*   PROTECT(ans = _new_Rle(ans_values, ans_lengths)); */
/*   UNPROTECT(3); */
/*   return ans; */
/* } */


/* SEXP Rle_constructor(SEXP values, SEXP lengths, SEXP check, SEXP buflength) */
/* { */
/*   int nvalues, buflength0; */
/*   const int *lengths_p; */

/*   nvalues = LENGTH(values); */
/*   if (LOGICAL(check)[0] && LENGTH(lengths) > 0) { */
/*     if (LENGTH(lengths) != nvalues) */
/*       error("'length(lengths)' != 'length(values)'"); */
/*     _sum_non_neg_ints(INTEGER(lengths), LENGTH(lengths), */
/* 		      "lengths"); */
/*   } */
/*   lengths_p = LENGTH(lengths) > 0 ? INTEGER(lengths) : NULL; */
/*   buflength0 = INTEGER(buflength)[0]; */
/*   return _logical_Rle_constructor(LOGICAL(values), nvalues, */
/* 				  lengths_p, buflength0); */
/*   return R_NilValue; */
/* } */
