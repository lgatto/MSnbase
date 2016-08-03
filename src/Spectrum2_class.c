#include "MSnbase.h"

/************************************************************
 * The Spectrum2 constructor in C.
 *
 */

/* This constructor uses functions defined in Rinternals.h */
static SEXP _new_Spectrum2(SEXP msLevel, SEXP peaksCount, SEXP rt,
			   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
			   SEXP mz, SEXP intensity, SEXP fromFile,
			   SEXP centroided, SEXP smoothed, SEXP polarity,
			   SEXP merged, SEXP precScanNum, SEXP precursorMz,
			   SEXP precursorIntensity, SEXP precursorCharge,
			   SEXP collisionEnergy)
{
  SEXP classdef, ans;
  if (asInteger(msLevel) < 2)
    error("_new_Spectrum2: msLevel should be >= 2, but I got: %d!\n", asInteger(msLevel));

  PROTECT(classdef = R_getClassDef("Spectrum2"));
  PROTECT(ans = R_do_new_object(classdef));
  /* Fill with data; should I do PROTECT here for all slots too???*/
  R_do_slot_assign(ans, install("msLevel"), msLevel);
  R_do_slot_assign(ans, install("msLevel"), msLevel);
  R_do_slot_assign(ans, install("peaksCount"), peaksCount);
  R_do_slot_assign(ans, install("rt"), rt);
  R_do_slot_assign(ans, install("acquisitionNum"), acquisitionNum);
  R_do_slot_assign(ans, install("scanIndex"), scanIndex);
  R_do_slot_assign(ans, install("tic"), tic);
  R_do_slot_assign(ans, install("mz"), mz);
  R_do_slot_assign(ans, install("intensity"), intensity);
  R_do_slot_assign(ans, install("fromFile"), fromFile);
  R_do_slot_assign(ans, install("centroided"), centroided);
  R_do_slot_assign(ans, install("smoothed"), smoothed);
  R_do_slot_assign(ans, install("polarity"), polarity);
  R_do_slot_assign(ans, install("merged"), merged);
  R_do_slot_assign(ans, install("precScanNum"), precScanNum);
  R_do_slot_assign(ans, install("precursorMz"), precursorMz);
  R_do_slot_assign(ans, install("precursorIntensity"), precursorIntensity);
  R_do_slot_assign(ans, install("precursorCharge"), precursorCharge);
  R_do_slot_assign(ans, install("collisionEnergy"), collisionEnergy);
  UNPROTECT(2);
  return(ans);
}

/*
 * Uses PROTECT for each newly generated SEXP
 */
static SEXP _new_Spectrum2_memsafe(SEXP msLevel, SEXP peaksCount, SEXP rt,
				   SEXP acquisitionNum, SEXP scanIndex, SEXP tic,
				   SEXP mz, SEXP intensity, SEXP fromFile,
				   SEXP centroided, SEXP smoothed, SEXP polarity,
				   SEXP merged, SEXP precScanNum, SEXP precursorMz,
				   SEXP precursorIntensity, SEXP precursorCharge,
				   SEXP collisionEnergy)
{
  SEXP classdef, ans;
  if (asInteger(msLevel) < 2)
    error("_new_Spectrum2_memsafe: msLevel should be >= 2, but I got: %d!\n",
	  asInteger(msLevel));

  PROTECT(classdef = R_getClassDef("Spectrum2"));
  PROTECT(ans = R_do_new_object(classdef));
  /* Fill with data; should I do PROTECT here for all slots too???*/
  R_do_slot_assign(ans, PROTECT(install("msLevel")), msLevel);
  R_do_slot_assign(ans, PROTECT(install("msLevel")), msLevel);
  R_do_slot_assign(ans, PROTECT(install("peaksCount")), peaksCount);
  R_do_slot_assign(ans, PROTECT(install("rt")), rt);
  R_do_slot_assign(ans, PROTECT(install("acquisitionNum")), acquisitionNum);
  R_do_slot_assign(ans, PROTECT(install("scanIndex")), scanIndex);
  R_do_slot_assign(ans, PROTECT(install("tic")), tic);
  R_do_slot_assign(ans, PROTECT(install("mz")), mz);
  R_do_slot_assign(ans, PROTECT(install("intensity")), intensity);
  R_do_slot_assign(ans, PROTECT(install("fromFile")), fromFile);
  R_do_slot_assign(ans, PROTECT(install("centroided")), centroided);
  R_do_slot_assign(ans, PROTECT(install("smoothed")), smoothed);
  R_do_slot_assign(ans, PROTECT(install("polarity")), polarity);
  R_do_slot_assign(ans, PROTECT(install("merged")), merged);
  R_do_slot_assign(ans, PROTECT(install("precScanNum")), precScanNum);
  R_do_slot_assign(ans, PROTECT(install("precursorMz")), precursorMz);
  R_do_slot_assign(ans, PROTECT(install("precursorIntensity")), precursorIntensity);
  R_do_slot_assign(ans, PROTECT(install("precursorCharge")), precursorCharge);
  R_do_slot_assign(ans, PROTECT(install("collisionEnergy")), collisionEnergy);
  UNPROTECT(21);
  return(ans);
}


/* SEXP _new_Spectrum2_from_C(int msLevel, int peaksCount, double rt, */
/* 			   int acqui) { */
/* } */

/*NOTE: this one does not ensure M/Z ordering!*/
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
 * This constructor ensures that M/Z-intensity pairs are ordered by M/Z. Also, it
 * calculates the TIC if not provided (0).
 */
SEXP Multi_Spectrum2_constructor_mz_sorted(SEXP msLevel,
					   SEXP peaksCount,
					   SEXP rt,
					   SEXP acquisitionNum,
					   SEXP scanIndex,
					   SEXP tic,
					   SEXP mz,
					   SEXP intensity,
					   SEXP fromFile,
					   SEXP centroided,
					   SEXP smoothed,
					   SEXP polarity,
					   SEXP merged,
					   SEXP precScanNum,
					   SEXP precursorMz,
					   SEXP precursorIntensity,
					   SEXP precursorCharge,
					   SEXP collisionEnergy,
					   SEXP nvalues,
					   SEXP check)
{
  int n = LENGTH(nvalues);
  int currentN = 0;
  int startN = 0;
  int calcTic;
  double currentTic;
  // Pointers
  double *p_rt, *p_tic, *p_merged, *p_precursorMz, *p_precursorIntensity,
    *p_collisionEnergy, *p_mz, *p_intensity;

  int *p_peaksCount, *p_acquisitionNum, *p_scanIndex, *p_polarity,
    *p_fromFile, *p_centroided, *p_smoothed, *p_nvalues,
    *p_precScanNum, *p_precursorCharge, *p_msLevel;

  // Creating pointers to the C-arrays.
  p_msLevel = INTEGER(msLevel);
  p_peaksCount = INTEGER(peaksCount);
  p_rt = REAL(rt);
  p_acquisitionNum = INTEGER(acquisitionNum);
  p_scanIndex = INTEGER(scanIndex);
  p_tic = REAL(tic);
  p_mz = REAL(mz);
  p_intensity = REAL(intensity);
  p_fromFile = INTEGER(fromFile);
  p_centroided = LOGICAL(centroided);
  p_smoothed = LOGICAL(smoothed);
  p_polarity = INTEGER(polarity);
  p_merged = REAL(merged);
  p_precScanNum = INTEGER(precScanNum);
  p_precursorMz = REAL(precursorMz);
  p_precursorIntensity = REAL(precursorIntensity);
  p_precursorCharge = INTEGER(precursorCharge);
  p_collisionEnergy = REAL(collisionEnergy);
  p_nvalues = INTEGER(nvalues);

  SEXP out = PROTECT(allocVector(VECSXP, n));
  for (int i = 0; i < n; i++) {
    double *p_cIntensity, *p_cMz, *p_orderMz;
    SEXP cMz, cIntensity, orderMz;
    currentN = p_nvalues[i];
    // Check if TIC is 0
    currentTic = p_tic[i];
    calcTic = 0;
    if (currentTic == 0)
      calcTic = 1;
    // Creating the mz and intensity vectors.
    cMz = PROTECT(allocVector(REALSXP, currentN));
    cIntensity = PROTECT(allocVector(REALSXP, currentN));
    // And the one we need for sorting
    orderMz = PROTECT(allocVector(REALSXP, currentN));
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
    // Sort the M/Z and intensity values.
    for (int j = 0; j < currentN; j++) {
      p_cMz[j] = p_orderMz[idx[j]];
      p_cIntensity[j] = p_intensity[startN + idx[j]];
      // Calculate TIC if was only 0
      if (calcTic > 0)
    	currentTic += p_cIntensity[j];
    }
    // Create the new Spectrum2
    SET_VECTOR_ELT(out, i, _new_Spectrum2(PROTECT(ScalarInteger(p_msLevel[i])),
					  PROTECT(ScalarInteger(p_peaksCount[i])),
					  PROTECT(ScalarReal(p_rt[i])),
					  PROTECT(ScalarInteger(p_acquisitionNum[i])),
					  PROTECT(ScalarInteger(p_scanIndex[i])),
					  PROTECT(ScalarReal(currentTic)),
					  cMz,
					  cIntensity,
					  PROTECT(ScalarInteger(p_fromFile[i])),
					  PROTECT(ScalarLogical(p_centroided[i])),
					  PROTECT(ScalarLogical(p_smoothed[i])),
					  PROTECT(ScalarInteger(p_polarity[i])),
					  PROTECT(ScalarReal(p_merged[i])),
					  PROTECT(ScalarInteger(p_precScanNum[i])),
					  PROTECT(ScalarReal(p_precursorMz[i])),
					  PROTECT(ScalarReal(p_precursorIntensity[i])),
					  PROTECT(ScalarInteger(p_precursorCharge[i])),
					  PROTECT(ScalarReal(p_collisionEnergy[i]))));
    UNPROTECT(19);
    startN = startN + currentN;
  }

  UNPROTECT(1);
  return(out);
}

/* This is definitely memory safe, i.e. uses tons of PROTECT calls */
SEXP Multi_Spectrum2_constructor_mz_sorted_memsafe(SEXP msLevel,
						   SEXP peaksCount,
						   SEXP rt,
						   SEXP acquisitionNum,
						   SEXP scanIndex,
						   SEXP tic,
						   SEXP mz,
						   SEXP intensity,
						   SEXP fromFile,
						   SEXP centroided,
						   SEXP smoothed,
						   SEXP polarity,
						   SEXP merged,
						   SEXP precScanNum,
						   SEXP precursorMz,
						   SEXP precursorIntensity,
						   SEXP precursorCharge,
						   SEXP collisionEnergy,
						   SEXP nvalues,
						   SEXP check)
{
  int n = LENGTH(nvalues);
  int currentN = 0;
  int startN = 0;
  int calcTic;
  double currentTic;
  // Pointers
  double *p_rt, *p_tic, *p_merged, *p_precursorMz, *p_precursorIntensity,
    *p_collisionEnergy, *p_mz, *p_intensity;
  double *p_cIntensity, *p_cMz, *p_orderMz;
  SEXP out, cMz, cIntensity, orderMz, spectrum;

  int *p_peaksCount, *p_acquisitionNum, *p_scanIndex, *p_polarity,
    *p_fromFile, *p_centroided, *p_smoothed, *p_nvalues,
    *p_precScanNum, *p_precursorCharge, *p_msLevel;

  // Creating pointers to the C-arrays.
  p_msLevel = INTEGER(msLevel);
  p_peaksCount = INTEGER(peaksCount);
  p_rt = REAL(rt);
  p_acquisitionNum = INTEGER(acquisitionNum);
  p_scanIndex = INTEGER(scanIndex);
  p_tic = REAL(tic);
  p_mz = REAL(mz);
  p_intensity = REAL(intensity);
  p_fromFile = INTEGER(fromFile);
  p_centroided = LOGICAL(centroided);
  p_smoothed = LOGICAL(smoothed);
  p_polarity = INTEGER(polarity);
  p_merged = REAL(merged);
  p_precScanNum = INTEGER(precScanNum);
  p_precursorMz = REAL(precursorMz);
  p_precursorIntensity = REAL(precursorIntensity);
  p_precursorCharge = INTEGER(precursorCharge);
  p_collisionEnergy = REAL(collisionEnergy);
  p_nvalues = INTEGER(nvalues);

  out = PROTECT(allocVector(VECSXP, n));
  for (int i = 0; i < n; i++) {
    currentN = p_nvalues[i];
    // Check if TIC is 0
    currentTic = p_tic[i];
    calcTic = 0;
    if (currentTic == 0)
      calcTic = 1;
    // Creating the mz and intensity vectors.
    cMz = PROTECT(allocVector(REALSXP, currentN));
    cIntensity = PROTECT(allocVector(REALSXP, currentN));
    // And the one we need for sorting
    orderMz = PROTECT(allocVector(REALSXP, currentN));
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
    // Sort the M/Z and intensity values.
    for (int j = 0; j < currentN; j++) {
      p_cMz[j] = p_orderMz[idx[j]];
      p_cIntensity[j] = p_intensity[startN + idx[j]];
      // Calculate TIC if was only 0
      if (calcTic > 0)
    	currentTic += p_cIntensity[j];
    }
    /* Protect each newly generated SEXP */
    PROTECT(spectrum = _new_Spectrum2_memsafe(PROTECT(ScalarInteger(p_msLevel[i])),
					      PROTECT(ScalarInteger(p_peaksCount[i])),
					      PROTECT(ScalarReal(p_rt[i])),
					      PROTECT(ScalarInteger(p_acquisitionNum[i])),
					      PROTECT(ScalarInteger(p_scanIndex[i])),
					      PROTECT(ScalarReal(currentTic)),
					      cMz,
					      cIntensity,
					      PROTECT(ScalarInteger(p_fromFile[i])),
					      PROTECT(ScalarLogical(p_centroided[i])),
					      PROTECT(ScalarLogical(p_smoothed[i])),
					      PROTECT(ScalarInteger(p_polarity[i])),
					      PROTECT(ScalarReal(p_merged[i])),
					      PROTECT(ScalarInteger(p_precScanNum[i])),
					      PROTECT(ScalarReal(p_precursorMz[i])),
					      PROTECT(ScalarReal(p_precursorIntensity[i])),
					      PROTECT(ScalarInteger(p_precursorCharge[i])),
					      PROTECT(ScalarReal(p_collisionEnergy[i]))));
    SET_VECTOR_ELT(out, i, spectrum);
    UNPROTECT(20); /*cMz, cIntensity, orderMz, spectrum*/
    startN = startN + currentN;
  }

  UNPROTECT(1);
  return(out);
}
