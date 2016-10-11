#include "MSnbase.h"

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

  PROTECT(classdef = R_getClassDef("Spectrum1"));
  PROTECT(ans = R_do_new_object(classdef));
  /* Set the slot values */
  PROTECT(R_do_slot_assign(ans, install("msLevel"), msLevel));
  PROTECT(R_do_slot_assign(ans, install("peaksCount"), peaksCount));
  PROTECT(R_do_slot_assign(ans, install("rt"), rt));
  PROTECT(R_do_slot_assign(ans, install("acquisitionNum"), acquisitionNum));
  PROTECT(R_do_slot_assign(ans, install("scanIndex"), scanIndex));
  PROTECT(R_do_slot_assign(ans, install("tic"), tic));
  PROTECT(R_do_slot_assign(ans, install("mz"), mz));
  PROTECT(R_do_slot_assign(ans, install("intensity"), intensity));
  PROTECT(R_do_slot_assign(ans, install("fromFile"), fromFile));
  PROTECT(R_do_slot_assign(ans, install("centroided"), centroided));
  PROTECT(R_do_slot_assign(ans, install("smoothed"), smoothed));
  PROTECT(R_do_slot_assign(ans, install("polarity"), polarity));

  UNPROTECT(14);
  //  UNPROTECT(2);
  return ans;
}

/* Constructor that adds the classVersion of Spectrum1 and Spectrum */
static SEXP _new_versioned_Spectrum1(SEXP msLevel, SEXP peaksCount, SEXP rt,
				     SEXP acquisitionNum, SEXP scanIndex,
				     SEXP tic, SEXP mz, SEXP intensity,
				     SEXP fromFile, SEXP centroided,
				     SEXP smoothed, SEXP polarity,
				     SEXP versions)
{
  SEXP classdef, ans, vers_classdef, vers;

  PROTECT(classdef = R_getClassDef("Spectrum1"));
  PROTECT(ans = R_do_new_object(classdef));
  /* Set the slot values */
  PROTECT(ans = R_do_slot_assign(ans, install("msLevel"), msLevel));
  PROTECT(ans = R_do_slot_assign(ans, install("peaksCount"), peaksCount));
  PROTECT(ans = R_do_slot_assign(ans, install("rt"), rt));
  PROTECT(ans = R_do_slot_assign(ans, install("acquisitionNum"), acquisitionNum));
  PROTECT(ans = R_do_slot_assign(ans, install("scanIndex"), scanIndex));
  PROTECT(ans = R_do_slot_assign(ans, install("tic"), tic));
  PROTECT(ans = R_do_slot_assign(ans, install("mz"), mz));
  PROTECT(ans = R_do_slot_assign(ans, install("intensity"), intensity));
  PROTECT(ans = R_do_slot_assign(ans, install("fromFile"), fromFile));
  PROTECT(ans = R_do_slot_assign(ans, install("centroided"), centroided));
  PROTECT(ans = R_do_slot_assign(ans, install("smoothed"), smoothed));
  PROTECT(ans = R_do_slot_assign(ans, install("polarity"), polarity));

  // Create the class version; version is supposed to be a NAMED LIST!!!
  PROTECT(vers_classdef = R_getClassDef("Versions"));
  PROTECT(vers = R_do_new_object(vers_classdef));
  PROTECT(vers = R_do_slot_assign(vers, install(".Data"), versions));
  // Assign the version
  PROTECT(ans = R_do_slot_assign(ans, install(".__classVersion__"), vers));

  UNPROTECT(18);
  return ans;
}

/* Constructor function settings also the classVersions. 'versions' has to be
 * a named list with names corresponding to the class name(s).
 */
SEXP Spectrum1_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP smoothed, SEXP polarity,
			   SEXP check, SEXP versions) {
  //int nvalues;
  SEXP ans;
  if (LOGICAL(check)[0]) {
    if (LENGTH(mz) != LENGTH(intensity)) {
      error("'length(intensity)' != 'length(mz)'");
    }
  }
  PROTECT(ans = _new_versioned_Spectrum1(msLevel, peaksCount, rt, acquisitionNum,
					 scanIndex, tic, mz, intensity, fromFile,
					 centroided, smoothed, polarity,
					 versions));
  UNPROTECT(1);
  return(ans);
}


/* Constructor that automatically sorts M/Z and intensity values by M/Z.
 * Adds in addition the class version(s) (see issue #163).
 */
SEXP Spectrum1_constructor_mz_sorted(SEXP msLevel, SEXP peaksCount,
				     SEXP rt, SEXP acquisitionNum,
				     SEXP scanIndex, SEXP tic,
				     SEXP mz, SEXP intensity,
				     SEXP fromFile, SEXP centroided,
				     SEXP smoothed, SEXP polarity,
				     SEXP check, SEXP versions)
{
  //int nvalues;
  SEXP ans, oMz, oInt;
  int nvals, i, calcTic;
  double theTic;

  calcTic = 0;
  theTic = asReal(tic);
  if (theTic == 0)
    calcTic = 1;
  //const int *lengths_p;
  nvals = LENGTH(mz);

  // Create pointers.
  double *pOint, *pOmz, *pIntensity, *pMz;

  //nvalues = LENGTH(intensity);
  if (LOGICAL(check)[0]) {
    if ( nvals != LENGTH(intensity)) {
      error("'length(intensity)' != 'length(mz)'");
    }
  }
  // initialize new M/Z, intensity vectors.
  PROTECT(oMz = allocVector(REALSXP, nvals));
  PROTECT(oInt = allocVector(REALSXP, nvals));
  // create pointers.
  pOmz = REAL(oMz);
  pOint = REAL(oInt);
  pMz = REAL(mz);
  pIntensity = REAL(intensity);
  int idx[nvals];
  _get_order_of_double_array(pMz, nvals, 0, idx, 0);
  // sort M/Z and intensity by M/Z
  for (i = 0; i < nvals; i++) {
    pOmz[i] = pMz[idx[i]];
    pOint[i] = pIntensity[idx[i]];
    if (calcTic)
      theTic += pOint[i];
  }
  PROTECT(ans = _new_versioned_Spectrum1(msLevel, peaksCount, rt,
					 acquisitionNum, scanIndex,
					 PROTECT(ScalarReal(theTic)),
					 oMz, oInt, fromFile, centroided,
					 smoothed, polarity, versions));
  UNPROTECT(4);
  return(ans);
  return R_NilValue;
}



/*************************************************************
 * Multi_Spectrum1_constructor_mz_sorted
 *
 * Simple C-function to create a list of Spectrum1 object based on the
 * submitted values.
 * All arguments (except mz, intensity) are supposed to be vectors of
 * length equal to the number of spectra. Argument nvalues defining the
 * number of M/Z and intensity values per spectrum.
 * This function is called by the R-function: Spectra1, argument checking is supposed
 * to take place there.
 * This constructor ensures that M/Z-intensity pairs are ordered by M/Z.
 * Adds in addition the class version(s), see issue #163.
 */
SEXP Multi_Spectrum1_constructor_mz_sorted(SEXP msLevel, SEXP peaksCount,
					   SEXP rt, SEXP acquisitionNum,
					   SEXP scanIndex, SEXP tic,
					   SEXP mz, SEXP intensity,
					   SEXP fromFile, SEXP centroided,
					   SEXP smoothed, SEXP polarity,
					   SEXP nvalues, SEXP check,
					   SEXP versions)
{
  int n = LENGTH(nvalues);
  int currentN = 0;
  int startN = 0;
  SEXP cMz, cIntensity, orderMz;
  double *p_rt, *p_tic;
  int *p_peaksCount, *p_acquisitionNum, *p_scanIndex, *p_polarity,
    *p_fromFile, *p_centroided, *p_smoothed, *p_nvalues;
  double *p_orderMz, *p_cMz, *p_cIntensity, *p_mz, *p_intensity;

  double currentTic;
  int calcTic;

  p_rt = REAL(rt);
  p_tic = REAL(tic);
  p_peaksCount = INTEGER(peaksCount);
  p_acquisitionNum = INTEGER(acquisitionNum);
  p_scanIndex = INTEGER(scanIndex);
  p_fromFile = INTEGER(fromFile);
  p_polarity = INTEGER(polarity);
  p_nvalues = INTEGER(nvalues);
  p_centroided = LOGICAL(centroided);
  p_smoothed = LOGICAL(smoothed);
  p_mz = REAL(mz);
  p_intensity = REAL(intensity);

  SEXP out = PROTECT(allocVector(VECSXP, n));
  for (int i = 0; i < n; i++) {
    currentN = p_nvalues[i];
    // Check if TIC is 0
    currentTic = p_tic[i];
    calcTic = 0;
    if (currentTic == 0)
      calcTic = 1;
    // Creating the mz and intensity vectors.
    PROTECT(orderMz = allocVector(REALSXP, currentN));
    PROTECT(cMz = allocVector(REALSXP, currentN));
    PROTECT(cIntensity = allocVector(REALSXP, currentN));
    // Create pointers to enable faster access (http://adv-r.had.co.nz/C-interface.html#c-vectors)
    p_orderMz = REAL(orderMz);
    p_cMz = REAL(cMz);
    p_cIntensity = REAL(cIntensity);
    // Fill the vector with the M/Z values of the current spectrum
    for (int j = 0; j < currentN; j++) {
      p_orderMz[j] = p_mz[startN + j];
    }
    // Get the order of the M/Z values.
    int idx[currentN];
    _get_order_of_double_array(p_orderMz, currentN, 0, idx, 0);
    // Sort the M/Z and intensity values for the current spectrum by M/Z
    for (int j = 0; j < currentN; j++) {
      p_cMz[j] = p_orderMz[idx[j]];
      p_cIntensity[j] = p_intensity[startN + (idx[j])];
      // Calculate TIC if provided one was 0
      if (calcTic > 0)
	currentTic += p_cIntensity[j];
    }
    // Create the new Spectrum1
    SET_VECTOR_ELT(out,
		   i,
		   _new_versioned_Spectrum1(msLevel,
					    PROTECT(ScalarInteger(p_peaksCount[i])),
					    PROTECT(ScalarReal(p_rt[i])),
					    PROTECT(ScalarInteger(p_acquisitionNum[i])),
					    PROTECT(ScalarInteger(p_scanIndex[i])),
					    /* ScalarReal(pTic[i]), */
					    PROTECT(ScalarReal(currentTic)),
					    cMz,
					    cIntensity,
					    PROTECT(ScalarInteger(p_fromFile[i])),
					    PROTECT(ScalarLogical(p_centroided[i])),
					    PROTECT(ScalarLogical(p_smoothed[i])),
					    PROTECT(ScalarInteger(p_polarity[i])),
					    versions));
    UNPROTECT(12);
    startN = startN + currentN;
  }

  UNPROTECT(1);
  return(out);
}

/* /\************************************************************* */
/*  * */
/*  * DONT USE THIS ONE!!! */
/*  * */
/*  * Multi_Spectrum1_constructor */
/*  * */
/*  * Simple C-function to create a list of Spectrum1 object based on the */
/*  * submitted values. */
/*  * All arguments (except mz, intensity) are supposed to be vectors of */
/*  * length equal to the number of spectra. Argument nvalues defining the */
/*  * number of M/Z and intensity values per spectrum. */
/*  * This function is called by the R-function: Spectra1, argument checking is supposed */
/*  * to take place there. */
/*  *\/ */
/* SEXP Multi_Spectrum1_constructor(SEXP msLevel, SEXP peaksCount, SEXP rt, */
/* 				 SEXP acquisitionNum, SEXP scanIndex, SEXP tic, */
/* 				 SEXP mz, SEXP intensity, SEXP fromFile, */
/* 				 SEXP centroided, SEXP smoothed, SEXP polarity, */
/* 				 SEXP nvalues, SEXP check) */
/* { */
/*   int n = length(nvalues); */
/*   int currentN = 0; */
/*   int startN = 0; */
/*   SEXP cMz, cIntensity; */
/*   SEXP out = PROTECT(allocVector(VECSXP, n)); */
/*   double *pRt, *pTic; */
/*   int *pPeaksCount, *pAcquisitionNum, *pScanIndex, *pPolarity, */
/*     *pFromFile, *pCentroided, *pSmoothed, *pNvalues; */

/*   pRt = REAL(rt); */
/*   pTic = REAL(tic); */
/*   pPeaksCount = INTEGER(peaksCount); */
/*   pAcquisitionNum = INTEGER(acquisitionNum); */
/*   pScanIndex = INTEGER(scanIndex); */
/*   pFromFile = INTEGER(fromFile); */
/*   pPolarity = INTEGER(polarity); */
/*   pNvalues = INTEGER(nvalues); */
/*   pCentroided = LOGICAL(centroided); */
/*   pSmoothed = LOGICAL(smoothed); */

/*   for (int i = 0; i < n; i++) { */
/*     // Creating the mz and intensity vectors. */
/*     currentN = pNvalues[i]; */
/*     PROTECT(cMz = allocVector(REALSXP, currentN)); */
/*     PROTECT(cIntensity = allocVector(REALSXP, currentN)); */
/*     for (int j = 0; j < currentN; j++) { */
/*       REAL(cMz)[j] = REAL(mz)[startN + j]; */
/*       REAL(cIntensity)[j] = REAL(intensity)[startN + j]; */
/*     } */

/*     SET_VECTOR_ELT(out, i, _new_Spectrum1(msLevel, */
/* 					  ScalarInteger(pPeaksCount[i]), */
/* 					  ScalarReal(pRt[i]), */
/* 					  ScalarInteger(pAcquisitionNum[i]), */
/* 					  ScalarInteger(pScanIndex[i]), */
/* 					  ScalarReal(pTic[i]), */
/* 					  cMz, */
/* 					  cIntensity, */
/* 					  ScalarInteger(pFromFile[i]), */
/* 					  ScalarLogical(pCentroided[i]), */
/* 					  ScalarLogical(pSmoothed[i]), */
/* 					  ScalarInteger(pPolarity[i]))); */
/*     UNPROTECT(2); */
/*     startN = startN + currentN; */
/*   } */

/*   UNPROTECT(1); */
/*   return(out); */
/* } */
