#include <R.h>
#include <Rdefines.h>
#include "MSnbase.h"
#include <stdio.h>

/************************************************************
 * The Spectrum1 constructor in C.
 *
 */

static SEXP _new_Spectrum1(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP smoothed, SEXP polarity)
{
  SEXP classdef, ans;

  PROTECT(classdef = MAKE_CLASS("Spectrum1"));
  PROTECT(ans = NEW_OBJECT(classdef));
  /* Set the slot values */
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
  SET_SLOT(ans, install("smoothed"), smoothed);
  SET_SLOT(ans, install("polarity"), polarity);

  UNPROTECT(2);
  return ans;
}

SEXP Spectrum1_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP smoothed, SEXP polarity,
			   SEXP check)
{
  //int nvalues;
  SEXP ans;
  //const int *lengths_p;

  //nvalues = LENGTH(intensity);
  if (LOGICAL(check)[0]) {
    if ( LENGTH(mz) != LENGTH(intensity)) {
      error("'length(intensity)' != 'length(mz)'");
    }
  }

  PROTECT(ans = _new_Spectrum1(msLevel, peaksCount, rt, acquisitionNum,
			       scanIndex, tic, mz, intensity, fromFile,
			       centroided, smoothed, polarity));
  UNPROTECT(1);
  return(ans);
  return R_NilValue;
}

/* Constructor that automatically sorts M/Z and intensity values by M/Z.
 */
SEXP Spectrum1_constructor_mz_sorted(SEXP msLevel, SEXP peaksCount,
				     SEXP rt, SEXP acquisitionNum,
				     SEXP scanIndex, SEXP tic,
				     SEXP mz, SEXP intensity,
				     SEXP fromFile, SEXP centroided,
				     SEXP smoothed, SEXP polarity,
				     SEXP check)
{
  //int nvalues;
  SEXP ans, oMz, oInt;
  int nvals, i;
  //const int *lengths_p;
  nvals = LENGTH(mz);

  //nvalues = LENGTH(intensity);
  if (LOGICAL(check)[0]) {
    if ( nvals != LENGTH(intensity)) {
      error("'length(intensity)' != 'length(mz)'");
    }
  }
  // initialize new M/Z, intensity vectors.
  PROTECT(oMz = NEW_NUMERIC(nvals));
  PROTECT(oInt = NEW_NUMERIC(nvals));
  int idx[nvals];
  _get_order_of_double_array(REAL(mz), nvals, 0, idx, 0);
  // sort M/Z and intensity by M/Z
  for (i = 0; i < nvals; i++) {
    REAL(oMz)[i] = REAL(mz)[idx[i]];
    REAL(oInt)[i] = REAL(intensity)[idx[i]];
  }
  PROTECT(ans = _new_Spectrum1(msLevel, peaksCount, rt, acquisitionNum,
			       scanIndex, tic, oMz, oInt, fromFile,
			       centroided, smoothed, polarity));
  UNPROTECT(3);
  return(ans);
  return R_NilValue;
}

/*************************************************************
 * Multi_Spectrum1_constructor
 *
 * Simple C-function to create a list of Spectrum1 object based on the
 * submitted values.
 * All arguments (except mz, intensity) are supposed to be vectors of
 * length equal to the number of spectra. Argument nvalues defining the
 * number of M/Z and intensity values per spectrum.
 * This function is called by the R-function: Spectra1, argument checking is supposed
 * to take place there.
 */
SEXP Multi_Spectrum1_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt,
				 SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
				 SEXP mz, SEXP intensity, SEXP fromFile,
				 SEXP centroided, SEXP smoothed, SEXP polarity,
				 SEXP nvalues, SEXP check)
{
  int n = length(nvalues);
  int currentN = 0;
  int startN = 0;
  SEXP cMz, cIntensity;
  SEXP out = PROTECT(allocVector(VECSXP, n));
  double *pRt, *pTic;
  int *pPeaksCount, *pAcquisitionNum, *pScanIndex, *pPolarity,
    *pFromFile, *pCentroided, *pSmoothed, *pNvalues;

  pRt = REAL(rt);
  pTic = REAL(tic);
  pPeaksCount = INTEGER(peaksCount);
  pAcquisitionNum = INTEGER(acquisitionNum);
  pScanIndex = INTEGER(scanIndex);
  pFromFile = INTEGER(fromFile);
  pPolarity = INTEGER(polarity);
  pNvalues = INTEGER(nvalues);
  pCentroided = LOGICAL(centroided);
  pSmoothed = LOGICAL(smoothed);

  for (int i = 0; i < n; i++) {
    // Creating the mz and intensity vectors.
    currentN = pNvalues[i];
    PROTECT(cMz = NEW_NUMERIC(currentN));
    PROTECT(cIntensity = NEW_NUMERIC(currentN));
    for (int j = 0; j < currentN; j++) {
      REAL(cMz)[j] = REAL(mz)[startN + j];
      REAL(cIntensity)[j] = REAL(intensity)[startN + j];
    }

    SET_VECTOR_ELT(out, i, _new_Spectrum1(msLevel,
					  ScalarInteger(pPeaksCount[i]),
					  ScalarReal(pRt[i]),
					  ScalarInteger(pAcquisitionNum[i]),
					  ScalarInteger(pScanIndex[i]),
					  ScalarReal(pTic[i]),
					  cMz,
					  cIntensity,
					  ScalarInteger(pFromFile[i]),
					  ScalarLogical(pCentroided[i]),
					  ScalarLogical(pSmoothed[i]),
					  ScalarInteger(pPolarity[i])));
    UNPROTECT(2);
    startN = startN + currentN;
  }

  UNPROTECT(1);
  return(out);
}


/*************************************************************
 * Multi_Spectrum1_constructor_mz_sorted
 *
 * Same as Multi_Spectrum1_constructor but ensures that M/Z-intensity pairs
 * are ordered by M/Z.
 */
SEXP Multi_Spectrum1_constructor_mz_sorted(SEXP msLevel, SEXP peaksCount,
					   SEXP rt, SEXP acquisitionNum,
					   SEXP scanIndex, SEXP tic,
					   SEXP mz, SEXP intensity,
					   SEXP fromFile, SEXP centroided,
					   SEXP smoothed, SEXP polarity,
					   SEXP nvalues, SEXP check)
{
  int n = length(nvalues);
  int currentN = 0;
  int startN = 0;
  SEXP cMz, cIntensity, currentMz;
  SEXP out = PROTECT(allocVector(VECSXP, n));
  double *pRt, *pTic;
  int *pPeaksCount, *pAcquisitionNum, *pScanIndex, *pPolarity,
    *pFromFile, *pCentroided, *pSmoothed, *pNvalues;

  pRt = REAL(rt);
  pTic = REAL(tic);
  pPeaksCount = INTEGER(peaksCount);
  pAcquisitionNum = INTEGER(acquisitionNum);
  pScanIndex = INTEGER(scanIndex);
  pFromFile = INTEGER(fromFile);
  pPolarity = INTEGER(polarity);
  pNvalues = INTEGER(nvalues);
  pCentroided = LOGICAL(centroided);
  pSmoothed = LOGICAL(smoothed);

  for (int i = 0; i < n; i++) {
    currentN = pNvalues[i];
    // Creating the mz and intensity vectors.
    PROTECT(currentMz = NEW_NUMERIC(currentN));
    PROTECT(cMz = NEW_NUMERIC(currentN));
    PROTECT(cIntensity = NEW_NUMERIC(currentN));
    // Fill the vector with the M/Z values of the current spectrum
    for (int j = 0; j < currentN; j++) {
      REAL(currentMz)[j] = REAL(mz)[startN + j];
    }
    // Get the order of the M/Z values.
    int idx[currentN];
    _get_order_of_double_array(REAL(currentMz), currentN, 0, idx, 0);
    // Sort the M/Z and intensity values for the current spectrum by M/Z
    for (int j = 0; j < currentN; j++) {
      REAL(cMz)[j] = REAL(currentMz)[idx[j]];
      REAL(cIntensity)[j] = REAL(intensity)[startN + (idx[j])];
    }
    // Create the new Spectrum1
    SET_VECTOR_ELT(out, i, _new_Spectrum1(msLevel,
					  ScalarInteger(pPeaksCount[i]),
					  ScalarReal(pRt[i]),
					  ScalarInteger(pAcquisitionNum[i]),
					  ScalarInteger(pScanIndex[i]),
					  ScalarReal(pTic[i]),
					  cMz,
					  cIntensity,
					  ScalarInteger(pFromFile[i]),
					  ScalarLogical(pCentroided[i]),
					  ScalarLogical(pSmoothed[i]),
					  ScalarInteger(pPolarity[i])));
    UNPROTECT(3);
    startN = startN + currentN;
  }

  UNPROTECT(1);
  return(out);
}

