#include <R.h>
#include <Rdefines.h>
#include "MSnbase.h"
#include <stdio.h>

/************************************************************
 * The Spectrum2 constructor in C.
 *
 */

static SEXP _new_Spectrum2(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP smoothed, SEXP polarity, SEXP merged,
			   SEXP precScanNum, SEXP precursorMz,
			   SEXP precursorIntensity, SEXP precursorCharge,
			   SEXP collisionEnergy)
{
  SEXP classdef, ans;

  PROTECT(classdef = MAKE_CLASS("Spectrum2"));
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
  SET_SLOT(ans, install("merged"), merged);
  SET_SLOT(ans, install("precScanNum"), precScanNum);
  SET_SLOT(ans, install("precursorMz"), precursorMz);
  SET_SLOT(ans, install("precursorIntensity"), precursorIntensity);
  SET_SLOT(ans, install("precursorCharge"), precursorCharge);
  SET_SLOT(ans, install("collisionEnergy"), collisionEnergy);

  UNPROTECT(2);
  return ans;
}

SEXP Spectrum2_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP smoothed, SEXP polarity, SEXP merged,
			   SEXP precScanNum, SEXP precursorMz,
			   SEXP precursorIntensity, SEXP precursorCharge,
			   SEXP collisionEnergy, SEXP check)
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

  PROTECT(ans = _new_Spectrum2(msLevel, peaksCount, rt, acquisitionNum,
			       scanIndex, tic, mz, intensity, fromFile,
			       centroided, smoothed, polarity, merged, precScanNum,
			       precursorMz, precursorIntensity, precursorCharge,
			       collisionEnergy));
  UNPROTECT(1);
  return(ans);
  return R_NilValue;
}

/*************************************************************
 * Multi_Spectrum2_constructor
 *
 * Simple C-function to create a list of Spectrum2 object based on the
 * submitted values.
 * All arguments (except mz, intensity) are supposed to be vectors of
 * length equal to the number of spectra. Argument nvalues defining the
 * number of M/Z and intensity values per spectrum.
 * This function is called by the R-function: Spectra1, argument checking is supposed
 * to take place there.
 */
SEXP Multi_Spectrum2_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt,
				 SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
				 SEXP mz, SEXP intensity, SEXP fromFile,
				 SEXP centroided, SEXP smoothed, SEXP polarity, SEXP merged,
				 SEXP precScanNum, SEXP precursorMz,
				 SEXP precursorIntensity, SEXP precursorCharge,
				 SEXP collisionEnergy, SEXP nvalues,
				 SEXP check)
{
  int n = length(nvalues);
  int currentN = 0;
  int startN = 0;
  SEXP cMz, cIntensity;
  SEXP out = PROTECT(allocVector(VECSXP, n));
  double *pRt, *pTic, *pMerged, *pPrecursorMz, *pPrecursorIntensity,
      *pCollisionEnergy;
  int *pPeaksCount, *pAcquisitionNum, *pScanIndex, *pPolarity,
      *pFromFile, *pCentroided, *pSmoothed, *pNvalues, *pPrecScanNum,
      *pPrecursorCharge, *pMsLevel;

  pRt = REAL(rt);
  pTic = REAL(tic);
  pPrecursorMz = REAL(precursorMz);
  pPrecursorIntensity = REAL(precursorIntensity);
  pMerged = REAL(merged);
  pCollisionEnergy = REAL(collisionEnergy);
  pMsLevel = INTEGER(msLevel);
  pPeaksCount = INTEGER(peaksCount);
  pAcquisitionNum = INTEGER(acquisitionNum);
  pScanIndex = INTEGER(scanIndex);
  pFromFile = INTEGER(fromFile);
  pPolarity = INTEGER(polarity);
  pNvalues = INTEGER(nvalues);
  pPrecScanNum = INTEGER(precScanNum);
  pPrecursorCharge = INTEGER(precursorCharge);
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

    SET_VECTOR_ELT(out, i, _new_Spectrum2(ScalarInteger(pMsLevel[i]),
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
					  ScalarInteger(pPolarity[i]),
					  ScalarReal(pMerged[i]),
					  ScalarInteger(pPrecScanNum[i]),
					  ScalarReal(pPrecursorMz[i]),
					  ScalarReal(pPrecursorIntensity[i]),
					  ScalarInteger(pPrecursorCharge[i]),
					  ScalarReal(pCollisionEnergy[i])));
    UNPROTECT(2);
    startN = startN + currentN;
  }

  UNPROTECT(1);
  return(out);
}

