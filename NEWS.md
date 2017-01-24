# MSnbase 2.1

## Changes in version 2.1.10
- (nothing yet)

## Changes in version 2.1.9
- When fully trimmed, an (empty) spectrum has peaksCount of 0L - see
  https://github.com/lgatto/MSnbase/issues/184 <2017-01-20 Fri>
- Add filterEmptySpectra,MSnExp method (see issue #181)
  <2017-01-20 Fri>
- Add a section about notable on-disk and in-memory differences (was
  issue #165) <2017-01-20 Fri>  

## Changes in version 2.1.8
- Remove order option altogether <2017-01-19 Thu> (superseeds setting
  default sorting using "auto" on R < 3.3 and "radix" otherwise
  <2017-01-03 Tue>)

## Changes in version 2.1.7
- Setting default sorting using "auto" on R < 3.3 and "radix"
  otherwise <2017-01-03 Tue>  
- filterMz returns an empty spectrum when no data is within the mz
  range (see issue #181) <2017-01-16 Mon>
- Performance improvement: a new private .firstMsLevel will
  efficiently return the first MS level in an MSnExp and
  OnDiskMSnExp. See issue #183 for details/background <2017-01-18 Wed>

## Changes in version 2.1.6

- Migrate io and dev vignettes to BiocStyle's html_document2 style
  <2016-12-23 Fri>
- Update show method to display class.
- Migrated to NEWS.md <2016-12-23 Fri>
- Update DESCRIPTION (and README) to reflect wider usage of MSnbase
  (replaced MS-based proteomics by mass spectrometry and proteomics)
  <2016-12-23 Fri>
 
## Changes in version 2.1.5

 - Fix (unexported) navMS example code <2016-12-14 Wed>
 
## Changes in version 2.1.4

 - Jo added netCDF support <2016-11-30 Wed>
 
## Changes in version 2.1.3

 - FeaturesOfInterest collections can now be assigned names -
   addresses issue #172 <2016-11-25 Fri>
 
## Changes in version 2.1.2

 - Update readMSnSet2 to save filename <2016-11-09 Wed>
 - Ensure that header information is read too if spectra data is
   loaded for onDiskMSnExp objects (see issue #170) <2016-11-24 Thu>
 
## Changes in version 2.1.1

 - Fix typo in impute man page <2016-10-19 Wed>
 - Cite Lazar 2016 in vignette imputation section <2016-10-28 Fri>
 
## Changes in version 2.1.0

 - New version for Bioconductor devel
 
# MSnbase 2.0.0

 - New version for Bioconductor release version 3.4
 
# MSnbase 1.99

## Changes in version 1.99.6
 - Reverting to old initialize,Spectrum (see issue #163)
   <2016-10-07 Fri>
 - Setting Spectrum class versions outside of prototype (see issue
   #163). For this, there is now a vector of class version in
   .MSnbaseEnv <2016-10-10 Mon>

## Changes in version 1.99.5

 - Add removeReporters,OnDiskMSnExp (see issue #161 for details)
   <2016-10-07 Fri>
 
## Changes in version 1.99.4

 - Fix bug in readMzTabData_v0.9 <2016-10-07 Fri>
 
## Changes in version 1.99.3

 - Injection time is now added to the header when reading mzML files
   using readMSData2 (see issue #159) <2016-10-04 Tue>
 
## Changes in version 1.99.2

 - Added isolationWindow,MSnExp method <2016-09-23 Fri>
 - Updated readMSData Spectrum2 class to support MS levels > 2
   <2016-09-30 Fri>
 
## Changes in version 1.99.1

 - Fix MS level test in quantify,OnDiskMSnExp to support MS2+ level
   quantitation <2016-09-14 Wed>
 
## Changes in version 1.99.0

 - Add normalize method for onDiskMSnExp <2016-06-07 Tue>.
 - Fix bug in readMSData <2016-06-07 Tue>.
 - Added documentation and unit test for trimMz <2016-06-01 Wed>.
 - Added trimMz method for onDiskMSnExp objects <2016-05-31 Tue>.
 - Added (internal) method spectrapply to apply a function to all
   spectra of an onDiskMSnExp object; does the data import,
   subsetting, application of lazy processing steps
   etc. <2016-05-27 Fri>.
 - Added a section to the MSnbase-development vignette
   <2016-05-27 Fri>.
 - Finished with all pSet inherited methods and docs <2016-05-26 Thu>.
 - Added rtlim argument to spectra method for onDiskMSnExp
   <2016-05-26 Thu>.
 - Methods rtime, tic, ionCount, polarity and acquisitionNum
   implemented <2016-05-25 Wed>.
 - Documentation added for most methods <2016-05-25 Wed>.
 - Methods peaksCount and spectra for onDiskMSnExp objects implemented
  <2016-05-24 Tue>.
 - Performance and validation tests for the Spectrum1 C-constructor
   <2016-05-23 Mon>.
 - Spectrum1 constructor implemented in C (presently not exported)
   <2016-05-20>
 - length, scanIndex, acquisitionNum and centroided implemented
   <2016-05-20>
 - Implemented fromFile and msLevel for onDiskMSnExp <2016-05-19 Thu>
 - Implemented an onDiskMSnExp class for on-the-fly data import
   <2016-05-19 Thu>
 - Fixed the validate method for pSet to play nicely with onDiskMSnExp
   objects <2016-05-19 Thu>
 - Added slot onDisk to pSet objects (TRUE for onDiskMSnExp objects,
   FALSE otherwise).  The getter method isOnDisk checks for the
   presence of the slot in the object (backward compatibility)
   <2016-05-19 Thu>
 - Implemented a ProcessingStep class that helps to keep track how
   (spectra) data should be processed on the fly <2016-05-18 Wed>
 - In onDiskMSnExp validity, check that the assaydata is empty
   <2016-06-28 Tue>
 - Pass neutralLoss in plot,Spectrum,Spectrum-method to
   .calculateFragments; fixes #146 <2016-08-12 Fri>
 - Allow the user to specify the `cex`, `lwd`, `pch` for peaks and
   fragments in plot,Spectrum,Spectrum-method; closes #148
   <2016-08-12 Fri>
 - Update centroided with an na.fail argument (see issue #150 for
   details) <2016-08-12 Fri>
 - Fix warning in readMgfData if TITLE contains multiple "=" <2016-08-24 Wed>

# MSnbase 1.21

## Changes in version 1.21.8

 - Update MSnSet validity method to guard agains empty string feature
   names <2016-06-20 Mon>
 - Simplify show,MSnExp method to work for various MS level cases and
   onDiskMSnExp - addressed issue #98 <2016-06-28 Tue>
 - removeMultipleAssignment also removes features that were not
   assigned (i.e. that have fcol (nprots) NA) <2016-07-02 Sat>
 - New smoothed slot/accessor/replacement methods <2016-07-08 Fri>
 - Update reporter masses and add TMT10 ETD/HCD <2016-07-20 Wed>
 - returning empty spectrum when fliterMz has empty range - see issue
   #134 <2016-07-22 Fri>
 - (mz, intensity) values are reordered based on order(mz) - see issue
   #135 <2016-07-26 Tue>
 - fix bug in bin,Spectrum - see issue #137 <2016-08-05 Fri>

## Changes in version 1.21.7

 - Update iPQF reference <2016-06-01 Wed>
 - Fix a bug in normalize method for MSnExp objects: assigning
   normalized spectra directly to assayData is not possible, as the
   environment is locked. See PR #91.
 - readMSData: if no phenodata is provided it creates an empty one
   with rownames corresponding to the file names. See PR #91.
 - Lock itraqdata's assaydata bindings <2016-06-08 Wed>

## Changes in version 1.21.6

 - MSmap unit test <2016-05-23 Mon>
 - Fix bug in as(MSmap, "data.frame") <2016-05-23 Mon>

## Changes in version 1.21.5

 - Added Johannes as contributor <2016-05-12 Thu>
 - Deprecate MzTab v0.9 <2016-05-19 Thu>
 - Fix old googlecode URLs to old MzTab <2016-05-19 Thu>

## Changes in version 1.21.4

 - More MzTab and Spectrum1 unit testing <2016-05-08 Sun>
 - Speed up readMSData (PR #86 by jotsetung) <2016-05-12 Thu>
 - Replace example file URL to use github instead of googlecode
    <2016-05-12 Thu>

## Changes in version 1.21.3

 - No fileNames replacement method <2016-05-07 Sat>
 - fileNames unit tests <2016-05-07 Sat>
 - add fileNames to class that had fileName accessor (MSmap, MzTab)
   <2016-05-07 Sat>

## Changes in version 1.21.2

 - Check for rownames/fnames in readMSnSet2 and unit test
   <2016-05-05 Thu>

## Changes in version 1.21.1

 - Fix wrong indexing in readMSdata, msLevel==1 (PR #85 by jotsetung)
   <2016-05-04 Wed>
 - grep/getEcols have a 'n' param specifying which line to grep/get
   <2016-05-04 Wed>

## Changes in version 1.21.0

 - Version bump for new Bioc devel

# MSnbase 1.20

## Changes in version 1.20.0

 - Version bump for Bioc release version 3.3

# MSnbase 1.19

## Changes in version 1.19.24
 - more unit tests and bug fixes <2016-04-30 Sat> <2016-05-02 Mon>

## Changes in version 1.19.23

 - more unit tests and bug fixes <2016-04-27 Wed> <2016-04-28 Thu>

## Changes in version 1.19.22

 - more unit tests <2016-04-23 Sat> <2016-04-24 Sun> <2016-04-26 Tue>
 - remove readIspyData functions <2016-04-23 Sat>
 - Fixed bug in error catching in
   utils.mergeSpectraAndIdentificationData <2016-04-23 Sat>

## Changes in version 1.19.21

 - nadata.R unit tests and bugs fixed <2016-04-22 Fri>

## Changes in version 1.19.20

 - Moved makeNaData[2] and whichNA (was in pRoloc) <2016-04-21 Thu>

## Changes in version 1.19.19

 - new naplot function to visualise missing values as a heatmap and
   barplots along the samples and features. <2016-04-04 Mon>

## Changes in version 1.19.18

 - fix typo in man <2016-04-02 Sat>

## Changes in version 1.19.17

 - Write support of mzTab has been dropped (writeMzTabData, makeMTD,
   makePEP and makePRT are now defunct) <2016-03-18 Fri>

## Changes in version 1.19.16

 - add estimateNoise,[Spectrum|MSnExp]-method; closes #78 <2016-03-10>
 - import a lot of functions from recommended packages, namely graphics, stats
   and utils to avoid many "Undefined global functions or variables" NOTEs in
   R CMD check <2016-03-10>

## Changes in version 1.19.15

 - limit readMSData unit test due to Windows-only error <2016-03-09 Wed>
 - fix unit test for utils.colSd(x, na.rm=TRUE)

## Changes in version 1.19.14

 - Document change in nQuants param fcol to groupBy <2016-03-02 Wed>

## Changes in version 1.19.13

 - Fixed bug in bin_Spectrum, reported by Weibo Xie <2016-02-16 Tue>
 - Added unit test for bug above <2016-02-16 Tue>
 - Merged 'apply functions columnwise by group' <2016-02-16 Tue>
 - In nQuants, the fcol argument has been replaced with groupBy to
   make the signature consistent with featureCV <2016-02-16 Tue>

## Changes in version 1.19.12

 - Moved polarity slot from Spectrum1 (0.1.0 -> 0.2.0) to Spectrum
   (0.2.0 -> 0.3.0) superclass; also bumped Spectrum2 class (0.1.0 ->
   0.2.0) and MSnExp (0.3.0 -> 0.3.1, to trace changes in
   Spectrum2). Wrote MSnExp and Spectrum2 updateObject
   methods. <2016-02-11 Thu>

## Changes in version 1.19.11

 - Fix trimws generic/methods with useAsDefault (see issue #75)
   <2016-02-02 Tue>
 - add exprs,MSnSet-method alias (since exprs is now exported)
   <2016-02-02 Tue>

## Changes in version 1.19.10

 - export exprs, fvarLabels, and sampleNames[<-] <2016-01-30 Sat>

## Changes in version 1.19.9

 - new trimws method for data.frames and MSnSets <2016-01-29 Fri>

## Changes in version 1.19.8

 - readMSnSet2 now also accepts a data.frame as input <2016-01-29 Fri>
 - selective ggplot2 import <2016-01-29 Fri>

## Changes in version 1.19.7

 - new sampleNames<- for pSet and MSnExp objects <2015-12-15 Tue>
 - Fix bug preventing to write MS1 to mgf (fixes issue #73 reported by
    meowcat) <2015-12-18 Fri>

## Changes in version 1.19.6

 - MSnExp feautreNames are now X01, X02 (0 after X) to maintain
   numerical sorting of ASCII characters; contributed by sgibb
   <2015-12-14 Mon>
 - Update MSnbase:::subsetBy to use split instead of lapply, which
   makes topN faster. This nevertheless changes the order of the
   resulting MSnSet (see issue #63 for details and background);
   contributed by sgibb <2015-12-14 Mon>

## Changes in version 1.19.5

 - Merged pull request #67 from lgatto/featureCV by sgibb: featureCV
 ignores its na.rm argument <2015-12-12 Sat>

## Changes in version 1.19.4

 - Replacement method for MSnSetList names <2015-11-24 Tue>

## Changes in version 1.19.3

 - New argument fcol to selectFeatureData to select feature variables
   using a vector <2015-10-28 Wed>

## Changes in version 1.19.2

 - new selectFeatureData function to subset the feature data
   <2015-10-24 Sat>

## Changes in version 1.19.1

 - Remove generics that are defined in ProtGenerics <2015-10-15 Thu>

## Changes in version 1.19.0

 - Bioc devel 3.3

# MSnbase 1.18

## Changes in version 1.18.0

 - Bioc release 3.2

# MSnbase 1.17

## Changes in version 1.17.16
 - new default parameter 'feature weight' in iPQF by Martina Fisher
   (see PR#65) <2015-10-07 Wed>

## Changes in version 1.17.15

 - Fix typo in plot man <2015-08-23 Sun>

## Changes in version 1.17.14

 - Fix typo in impute man <2015-08-17 Mon>

## Changes in version 1.17.13

 - partly rewrite writeMgfData <2015-05-16 Thu>
 - initial hmap function <2015-07-16 Thu>
 - fix bug in plotting MS1 spectra (closes issue #59) <2015-07-16 Thu>
 - new image implementation, based on @vladpetyuk's
   vp.misc::image_msnset <2015-07-25 Sat>
 - Changed the deprecated warning to a message when reading MzTab data
   version 0.9, as using the old reader can not only be achived by
   accident and will be kept for backwards file format compatibility
   <2015-07-30 Thu>

## Changes in version 1.17.12

 - fix show MIAPE when some otherInfo at NA <2015-07-15 Wed>

## Changes in version 1.17.11

 - adding unit tests <2015-07-01 Wed>
 - fix abundance column selection when creating MSnSet form MzTab
   <2015-07-01 Wed>
 - new mzTabMode and mzTabType shortcut accessors for mode and type of
   an mzTab data <2015-07-01 Wed>

## Changes in version 1.17.10

 - Fix URL <2015-06-30 Tue>

## Changes in version 1.17.9

 - calculateFragments' "neutralLoss" argument is now a list (was a logical
   before), see #47. <2015-06-24 Wed>
 - add defaultNeutralLoss() function to fill calculateFragments' modified
  "neutralLoss" argument, see #47. <2015-06-24 Wed>

## Changes in version 1.17.8

 - coercion from IBSpectra to MSnSet, as per user request
   <2015-06-23 Tue>
 - new iPQF combineFeature method <2015-06-24 Wed>

## Changes in version 1.17.7

 - Fix support of identification-only mzTab files <2015-06-22 Mon>

## Changes in version 1.17.6

 - Export metadata,MzTab-method <2015-06-19 Fri>
 - Replace spectra,MzTab-method by psms,MzTab-method <2015-06-20 Sat>
 - Change the meaning of calculateFragments' "modifications" argument. Now the
   modification is added to the mass of the amino acid/peptide. Before it was
   replaced. <2015-06-21 Sun>
 - calculateFragments gains the feature to handle N-/C-terminal modifications,
   see #47. <2015-06-21 Sun>
 - update readMzTabData example <2015-06-22 Mon>

## Changes in version 1.17.5

 - new MzTab class to store a simple parsing of an mzTab version 1
   data file. See ?MzTab for details. <2015-06-16 Tue>

## Changes in version 1.17.4

 - New lengths method for FoICollection instances <2015-06-06 Sat>
 - New image2 function for matrix object, that behaves like the method
   for MSnSets <2015-06-10 Wed>
 - image,MSnSet labels x and y axis as Samples and Features
   <2015-06-10 Wed>
 - fixed bug in purityCorrect, reported by Dario Strbenac
   <2015-06-11 Thu>

## Changes in version 1.17.3

 - iTRAQ 8-plex correction factors and impurity matrix
   <2015-05-22 Fri>

## Changes in version 1.17.2

 - new filterZero function <2015-05-01 Fri>

## Changes in version 1.17.1

 - new MSnSetList class <2015-04-19 Sun>
 - new commonFeatureNames function <2015-04-14 Tue>
 - new compareMSnSets function <2015-04-19 Sun>
 - splitting and unsplitting MSnSets/MSnSetLists <2015-04-19 Sun>

## Changes in version 1.17.0

 - new devel version for Bioc 3.2

# MSnbase 1.15

## Changes in version 1.15.18
 - fix failing test_MSnExp::readMSData unit test on Windows i386
   (@.cache$size being different on that arch) [2015-04-14 Tue]
 - merge @vladpetyuk PR #50 fix of combine features bug [2015-04-14 Tue]

## Changes in version 1.15.17
 - add TMT10 paragraph and fig to demo vignette [2015-04-09 Thu]

## Changes in version 1.15.16
 - support for TMT10 plex [2015-04-08 Wed]

## Changes in version 1.15.15
 - using serial parallelisation during quantitation in unit tests
   [2015-04-02 Thu]

## Changes in version 1.15.14
  - new msnset data, used in various examples instead of quantifying
    the itraqdata experiment over and over again [2015-04-01 Wed]

## Changes in version 1.15.13
 - improve nbavg imputation description and add example [2015-03-22 Sun]
 - reduce compareSpectra example timing [2015-03-30 Mon]

## Changes in version 1.15.12
 - average neighbour imputation for ordered fractions along a gradient
   [2015-03-21 Sat]

## Changes in version 1.15.11
 - update malformed Description [2015-03-20 Fri]

## Changes in version 1.15.10
 - update compareSpectra man [2015-03-19 Thu]
 - ExpressionSet <-> MSnSet coercion [2015-03-19 Thu]

## Changes in version 1.15.9

 - use S4Vectors' isEmpty generic [2015-03-10 Tue]

## Changes in version 1.15.8

 - Improved imputation section in demo vignette [2015-03-05 Thu]

## Changes in version 1.15.7
 - Importing ProtGenerics [2015-02-28 Sat]
 - selective import of MALDIquant [2015-03-02 Mon]

## Changes in version 1.15.6

 - add intensity column to the calculateFragments output; closes #47;1
   [2015-02-02 Mon]
 - add method argument to calculateFragments to allow the user choosing the
   highest/closest peak or all peaks in the tolerance range; closes #47;2
   [2015-02-09 Mon]
 - add neutralLoss argument to calculateFragments and calculate loss of water
   and ammonia; closes #47:3 [2015-02-19 Thu]
 - new imputation methods via imputeLCMD and norm [2015-02-09 Mon]
 - vignette updates [2015-02-09 Mon]
 - imputation unit test [2015-02-24 Tue]

## Changes in version 1.15.5

 - update dependency to mzID >= 1.5.2 [2015-01-28 Wed]
 - rewrite addIdentificationData and change its signature [2015-01-28 Wed]
 - add methods addIdentificationData which work on MSnExp and MSnSets using
   filenames (characters), mzID objects and data.frames; see #42; closes #45;
   [2015-01-28 Wed]
 - add a section about MSmaps in the vignette [2015-02-02 Mon]
 - MSmap has a new zeroIsNA argument to set all 0 values of the map to
   NA. This simplifies the resulting plot3D figure. [2015-02-02 Mon]

## Changes in version 1.15.4

 - partly rewrite readMgfData [2015-01-19 Mon]
 - Typo in addIdentification data man [2015-01-20 Tue]
 - use Biocstyle [2015-01-23 Fri]
 - replace require with requireNamespace [2015-01-23 Fri]

## Changes in version 1.15.3

- plot,Spectrum2,character to add fragment ions based on peptide
  sequence [2014-11-20 Thu]
- update vignette with above [2014-11-20 Thu]
- comment parallel code in quantify man [2015-01-09 Fri]

## Changes in version 1.15.2

 - merged sgibb's pull request exporting overloaded methods as well
   [2014-11-02 Sun]
 - updated MSnSet validity, checking that exprs(.) is a matrix
   [2014-11-12 Wed]
 - fix error (missing header extraction) is MSmap example
   [2014-11-18 Tue]

## Changes in version 1.15.1

 - Fixing error when id file has no spectrumFile info (see issue #39)
   and return a warning (instead of an error) when the file used to
   create the MSnExp/MSnSet and mzid file were different
   [2014-10-15 Wed]

## Changes in version 1.15.0

 - Devel version for Bioc 3.1

# MSnbase 1.14

## Changes in version 1.14.0
 - Release version for Bioc 3.0

# MSnbase 1.13

## Changes in version 1.13.16
 - width generic if non-existant [2014-09-03 Wed]
 - fix undocumented S4 methods warnings [2014-09-27 Sat]

## Changes in version 1.13.15
 - msMap(MSmap)<- method [2014-08-12 Tue]
 - Typo in MIAPE man [2014-08-29 Fri]

## Changes in version 1.13.14
 - Fix xic example [2014-08-08 Fri]
 - importing lattice, ggplot2 (was depends) and suggesting rgl
   [2014-08-09 Sat]
 - MSmap infrastructure [2014-08-09 Sat]

## Changes in version 1.13.13
 - fix issue with readMSnSet2 without fdata [2014-07-28 Mon]

## Changes in version 1.13.12
 - Don't import width from IRanges [2014-07-23 Wed]
 - qual slot is populated again [2014-07-23 Wed]

## Changes in version 1.13.11
 - remove Vennerable::Venn from example and DESCRIPTION [2014-06-19 Thu]

## Changes in version 1.13.10
 - new averageMSnSet function to generate an average over a list of
   MSnSets. [2014-06-17 Tue]
 - new non-parametric coefficent of variation function [2014-06-17 Tue]
 - Using/importing IRanges::width [2014-06-18 Wed]

## Changes in version 1.13.9
 - new listOf helper function [2014-06-16 Mon]
 - compfnames to compare and document differences in MSnSet feature
   names [2014-06-16 Mon]
 - suggesting Vennerable (compfnames example) and roxygen2 (generate
   rd) [2014-06-16 Mon]

## Changes in version 1.13.8
 - Add recommended biocView [2014-06-05 Thu]

## Changes in version 1.13.7
 - Bug tracking link [2014-05-26 Mon]

## Changes in version 1.13.6
 - new isEmpty method for Spectrum instances [2014-05-14 Wed]
 - removePeaks does not try for empty Spectra [2014-05-14 Wed]
 - isEmpty unit test [2014-05-14 Wed]

## Changes in version 1.13.5
 - export .get.amino.acids function [2014-05-08 Thu]
 - removePeaks for centroided data [2014-05-08 Thu]
 - add new exported function get.atomic.mass [2014-05-08 Thu]
 - Spectrum[2] prototype sets centroided=FALSE by default (instead of
   logical()) [2014-05-08 Thu]
 - fnamesIn also supports y = "data.frame" [2014-05-14 Wed]
 - Using BiocParallel for parallel support; replaced parallel argument
   with BPPARAM [2014-05-14 Wed]

## Changes in version 1.13.4

 - Document difference between traceable and non-traceable
   FeaturesOfInterest instances [2014-04-30 Wed]
 - ignoring desciptions with length > 1 [2014-04-30 Wed]

## Changes in version 1.13.3

 - FeaturesOfInterest [2014-04-29 Tue]

## Changes in version 1.13.2

 - add compareSpectra to compare Spectrum objects [2014-04-09 Wed]
 - add bin method for Spectrum objects [2014-04-09 Wed]
 - recreate inst/extdata/msx.rda for R 3.1 and Biobase 2.24 [2014-04-13 Sun]
 - add plot,Spectrum,Spectrum method [2014-04-14 Mon]
 - add calculateFragments,character,Spectrum and
   calculateFragments,character,missing methods [2014-04-15 Tue]
 - updated readMSData test unit to ignore object@.__classVersion__
   that fails with latest R/Bioc versions [2014-04-15 Tue]
 - formatRt conversion from 'mm:sec' to sec and unit test [2014-04-17 Thu]
 - Update nprot/npsm.prot/npsm.pep/npep.prot feature variables and
   relevant man/tests/argument defaults [2014-04-17 Thu]

## Changes in version 1.13.1

 - add precursor method to Spectrum2 normalisation method [2014-04-08 Tue]
 - add smooth for MSnExp and Spectrum classes [2014-04-10 Thu]
 - add pickPeaks for MSnExp and Spectrum classes [2014-04-10 Thu]
 - updated show,MSnExp to display only first/last files when > 2 [2014-04-11 Fri]

## Changes in version 1.13.0

- New devel version for Bioc 3.0

# MSnbase 1.11

## Changes in version 1.11.14
 - update dependency to R >= 3.1 [2014-04-05 Sat]

## Changes in version 1.11.13
 - Document [get|grep]Ecols in io vignette [2014-03-31 Mon]
 - typo in readMSnSet man [2014-03-31 Mon]
 - updated affiliation in vignettes [2014-03-31 Mon]

## Changes in version 1.11.12
 - Fixed a bug in readMzTabData reported by Hendrik Weisser [2014-03-26 Wed]

## Changes in version 1.11.11
 - precomputed msx test data now has id data [2014-03-25 Tue]
 - quantitation unit tests [2014-03-25 Tue]
 - quantify method now accepts label-free methods [2014-03-25 Tue]

## Changes in version 1.11.10

 - import mzID [2014-03-20 Thu]
 - dummyiTRAQ id extdata [2014-03-21 Fri]
 - add addIdentificationData method for MSnExp and MSnSet [2014-03-21 Fri]
 - using data from pRolocdata (was pRoloc) [2014-03-23 Sun]
 - added removeNoId,MSnExp,MSnSet methods [2014-03-23 Sun]
 - added idSummary,MSnExp,MSnSet methods [2014-03-24 Mon]
 - Not generating different sample per file when reading raw
   data. pData(.) now has systematically 1 row if not specifically
   provided by the user. A message is also reported by validity,pSet
   if row(pData(.)) > 1. [2014-03-24 Mon]
 - NAnnotatedDataFrame now has default multiplex 1 [2014-03-25 Tue]
 - NAnnotatedDataFrame unit tests [2014-03-25 Tue]

## Changes in version 1.11.9

 - write.exprs can have fDataCol or fcol (for consistence) [2014-03-17 Mon]
 - Fixing bug in combineFeatures(..., is.character(groupBy)) [2014-03-19 Wed]
 - fixed combineFeatures [2014-03-20 Thu]
 - added example, test and doc for combineFeatures with list [2014-03-20 Thu]

## Changes in version 1.11.8

 - adding redundancy handling to combineFeatures (by vladpetyuk, pull
   request #18) [2014-03-14 Fri]
 - updated combineFeatures signature to accomodate above changes
   [2014-03-14 Fri]
 - updated unit tests for new testhat 0.8 [2014-03-14 Fri]

## Changes in version 1.11.7

 - NA

## Changes in version 1.11.6

 - add corresponding xcms functions to the chromatogram and xic
   manual page [2014-02-21 Fri]
 - new bpca imputation methods [2014-02-27 Thu]
 - replacing stop_on_error with option in vignette [2014-02-27 Thu]

## Changes in version 1.11.5

 - typo in MSnSet droplevels man [2014-01-27 Mon]
 - typo in MSnbase-demo vignette [2014-02-20 Thu]
 - fix BPI legend in chromatogram [2014-02-20 Thu]

## Changes in version 1.11.4

 - passing ... to sweep when normalising [2013-12-08 Sun]
 - updated makeMTD to accomodate new MS ontology [2013-12-23 Mon]

## Changes in version 1.11.3

 - updated mzTab example files to new url [2013-11-15 Fri]
 - warning about mzTab versions [2013-11-15 Fri]

## Changes in version 1.11.2

 - move inst/doc to vignettes [2013-10-19 Sat]

## Changes in version 1.11.1

 - document na.rm in combineFeatures Rd [2013-10-18 Fri]

## Changes in version 1.11.0

 - New devel version for Bioc 2.14

# MSnbase 1.9

## Changes in version 1.9.12

 - fix MSnSet -] ExpressionSet [2013-10-13 Sun]
 - MSnSet -] ExpressionSet unit test [2013-10-13 Sun]

## Changes in version 1.9.11

 - MIAPE to MIAME conversion [2013-10-11 Fri]
 - proper MIAME when MSnSet -] ExpressionSet [2013-10-11 Fri]

## Changes in version 1.9.10

 - faster plotMzDelta [2013-09-28 Sat]
 - faster plotMzDelta for mzRramp instances [2013-09-29 Sun]
 - chromatogram method [2013-10-04 Fri]
 - plotMzDelta has subset arg [2013-10-04 Fri]
 - xic method [2013-10-04 Fri]
 - suggesting msdata for chromatogram example [2013-10-04 Fri]
 - renamed plotting arg 'centroided.' [2013-10-04 Fri]

## Changes in version 1.9.9

 - typo in filterNA Rd [2013-09-18 Wed]
 - writeMgfData now has a progress bar [2013-09-24 Tue]
 - centroided(MSnExp) <- TRUE now allowed [2013-09-24 Tue]

## Changes in version 1.9.8

 - using new.env(parent=emptyenv()) to get rid of enclosing env
   when creating new MSnExps [2013-09-17 Tue]
 - new (private) MSnExp.size function [2013-09-17 Tue]

## Changes in version 1.9.7

 - Passing ... to read.table in MSnbase:::readIspy[Silac|15N]Data [2013-09-16 Mon]
 - QualityControl biocView [2013-09-16 Mon]

## Changes in version 1.9.6

 - new as.data.frame.MSnSet method [2013-08-16 Fri]
 - new ms2df function [2013-08-16 Fri]
 - new getEcols and grepEcols helpers for readMSnSet2 [2013-08-16 Fri]

## Changes in version 1.9.5

 - typo in Author[s]@R [2013-05-15 Wed]

## Changes in version 1.9.4

 - new simple MSnSet constructor [2013-05-07 Tue]

## Changes in version 1.9.3

 - Using knitr as VignetteEngine [2013-04-29 Mon]
 - Remove LazyLoad from DESCRIPTION,
   which is default nowadays [2013-04-29 Mon]
 - knitr dependency ] 1.1.0 for VignetteBuilder [2013-04-29 Mon]
 - Adding MSnSet creating sections in io vignette [2013-04-29 Mon]
 - new readMSnSet2 function [2013-04-30 Tue]

## Changes in version 1.9.2

 - clean has now a all param (default FALSE is retain original behavious)
   to remove all 0 intensity values [2013-04-17 Wed]
 - using BiocGenerics::normalize [2013-04-25 Thu]

## Changes in version 1.9.1

 - new logging utility function to update an MSnSet's
   processingData(object)$processing [2013-03-29 Fri]
 - Proper logging in t.MSnSet [2013-03-29 Fri]

## Changes in version 1.9.0

 - new Bioc 2.13 devel version

# MSnbase 1.8

## Changes in version 1.8.0

 - new Bioc 2.12 stable version

# MSnbase 1.7

## Changes in version 1.7.25

 - updated itraqdata to fix issue in vignette when
   combine(exp1, exp2) and different MIAPE versions
   [2013-03-22 Fri]

## Changes in version 1.7.24

 - Mention scale in vignette [2013-03-02 Sat]
 - exprsToRatio matrix method [2013-03-20 Wed]

## Changes in version 1.7.23

 - new private nologging function [2013-02-21 Thu]
 - adding total number of features on plotNA [2013-02-22 Fri]
 - updated msnbase.r [2013-02-26 Tue]

## Changes in version 1.7.22

 - msnbase.r na.rm arg [2013-02-20 Wed]

## Changes in version 1.7.21

 - Added impute method  [2013-02-19 Tue]

## Changes in version 1.7.20

 - Explicitating that normalise and normalize are the same
   methods in the man. [2013-02-13 Wed]

## Changes in version 1.7.19

 - adding MIAPE and pSet accessors: analyserDetails, analyzerDetails,
   ionSourceDetails, instrumentModel, instrumentManufacturer,
   instrumentCustomisations [2013-02-12 Tue]
 - switching back to analyserDetails slot [2013-02-12 Tue]

## Changes in version 1.7.18

 - readMgfData now supports comments and
   PEPMASS with precursor mz and intensity,
   requested by Thomas Taus [2013-02-11 Mon]
 - improved and running read/writeMgfData
   example [2013-02-11 Mon]
 - new scanIndex accessor method [2013-02-11 Mon]

## Changes in version 1.7.17

 - added a analyzer accessors/slot to accomodate new mzTab
   files with more meta-data  [2013-01-30 Wed]

## Changes in version 1.7.16

 - fixing knitr 1.0 compatibility [2013-01-15 Tue]

## Changes in version 1.7.15

 - new scale method [2013-01-11 Fri]
 - renaming scale.mean and scale.median normalisation
   methods to center.mean and centre.median [2013-01-11 Fri]

## Changes in version 1.7.14

 - new unexported readIspy15NData [2013-01-09 Wed]
 - min.int readIspy[Silac|15N]Data arg [2013-01-11 Fri]

## Changes in version 1.7.13

 - msnbase.r v0.1.1 with oh (help) arg [2013-01-08 Tue]
 - msnbase.r coerce ob arg to numeric [2013-01-08 Tue]
 - testing if any features left in readIspyData [2013-01-08 Tue]

## Changes in version 1.7.12

 - updated makeImpuritiesMatrix to create matrix from
   csv file with correction factors [2012-12-23 Sun]
 - makeImpuritiesMatrix test [2012-12-23 Sun]
 - readIspyData: message instead of warning if NA in
   featureData [2012-12-24 Mon]
 - Added msnbase.r script [2012-12-24 Mon]

## Changes in version 1.7.11

 - new 'pattern' argument to filterNA [2012-12-15 Sat]
 - vignette and man updates [2012-12-15 Sat]
 - filterNA(, pattern) tests [2012-12-15 Sat]

## Changes in version 1.7.10

 - new droplevels.MSnSet S3 method [2012-12-14 Fri]
 - fixed errors in vignette and udpates [2012-12-14 Fri]
 - vignette build stops in case of error [2012-12-14 Fri]

## Changes in version 1.7.9

 - Updating processing data on readIspyData [2012-12-05 Wed]
 - filterNA has a droplevels arg [2012-12-05 Wed]
 - featureCV's default cv.norm is 'sum' now [2012-12-11 Tue]
 - fixed featureCV for 1 sample [2012-12-11 Tue]

## Changes in version 1.7.8

 - new featureCV function [2012-12-04 Tue]
 - more MSnSet combineFeatures tests [2012-12-04 Tue]
 - new TMT6 impurity matrix and fixed
   purityCorrect [2012-12-05 Wed]
 - combineFeatures now automatically computes
   feature CVs (using featureCV) and collates this
   in featureData [2012-12-05 Wed]
 - new exprsToRatios method (moved from pRoloc) [2012-12-05 Wed]
 - initial implementation of impurity correction using
   Cramer's rule (see MSnbase:::cramer4) [2012-12-05 Wed]

## Changes in version 1.7.7

 - added scale.mean and scale.median MSnSet normalisation
   method [2012-11-30 Fri]
 - improvements to readMSData [2012-11-30 Fri]
 - small updates to caching code, max level 2 [2012-11-30 Fri]
 - readMSdata test [2012-12-01 Sat]

## Changes in version 1.7.6

 - Fixed bug in readIspyData, reported by
   Claire Mulvey [2012-11-27 Tue]

## Changes in version 1.7.5

 - dropping levels in readIspySilacData [2012-11-06 Tue]
 - fixed plotNA [2012-11-09 Fri]

## Changes in version 1.7.4

 - exporting log method [2012-11-02 Fri]
 - private readIspySilacData function [2012-11-05 Mon]
 - updating '['-MSnSet to log intial/final dims [2012-11-05 Mon]

## Changes in version 1.7.3

 - updating readMzTabData to properly capture
   experiment description [2012-10-12 Fri]

## Changes in version 1.7.2

 - fixed readMSData for MS1 [2012-10-08 Mon]

## Changes in version 1.7.1

 - fixed vignettes [2012-10-02 Tue]

## Changes in version 1.7.0

 - Version bump for next devel release [2012-10-01 Mon]

# MSnbase 1.5

## Changes in version 1.5.25

 - fixed bug when quantifying exp of length 1 (reported
   by Colin Archer), added test [2012-09-26 Wed]
 - fixed parallel default to FALSE [2012-09-26 Wed]

## Changes in version 1.5.24

 - updated clean, removePeaks, combineFeatures, purityCorrect,
   trimMz, plot, MSnSet-class examples to not use readMSData [2012-09-25 Tue]
 - fixed clean,MSnExp-method [2012-09-25 Tue]
 - space in log message in extractPrecSpectra [2012-09-25 Tue]

## Changes in version 1.5.23

 - updated quantify documentation [2012-09-24 Mon]
 - changed all foo.class functions to foo_class [2012-09-24 Mon]

## Changes in version 1.5.22

 - setting parallel default to FALSE [2012-09-22 Sat]
 - enhanced parallel in DESCRIPTION and detectCores()
   passed to registerDoMC [2012-09-23 Sun]

## Changes in version 1.5.21

 - removed registerDoMC no visible global function NOTE [2012-09-21 Fri]
 - fixed NOTE about xvarname which was a bug [2012-09-21 Fri]

## Changes in version 1.5.20

 - an immediate warning is thrown is any(centroided(object))
   in quantify.MSnExp [2012-09-15 Sat]
 - mzTab file and loading time is now recorded in
   processingData [2012-09-15 Sat]
 - temporarily dontrun'ing' the mzTab read/write examples
   as -LS is down (and breaks rols) (note: modifed Rd files,
   not roxygen template) [2012-09-15 Sat]

## Changes in version 1.5.19

 - updated all ReporterIons rda data [2012-09-13 Thu]
 - Using && in testing parallel, require(foreach and doMC) [2012-09-14 Fri]

## Changes in version 1.5.18

 - new log method for MSnSet instances [2012-09-13 Thu]
 - new MAplot methods using generic and ma.plot/mva.pairs
   from affy [2012-09-13 Thu]

## Changes in version 1.5.17

 - changed TMT7[7] mass from 229.26 to 230.17 and
   ReporterIons descriptions [2012-09-12 Wed]

## Changes in version 1.5.16

 - parallel quantify is now always set to
   FALSE on Windows, fixing example checking
   issues [2012-09-11 Tue]
 - fixed types in plotMzDelta man [2012-09-11 Tue]

## Changes in version 1.5.15

 - updating code to ggplot2 v0.9.2 [2012-09-07 Fri]

## Changes in version 1.5.14

 - added platform test to use doMC (d.tenenbaum) [2012-08-31 Fri]

## Changes in version 1.5.13

 - updated fillUp function [2012-08-15 Wed]
 - added tikzDevice to suggests [2012-08-16 Thu]
 - tikzDevice no longer on CRAN - removing from Suggests
   and using pdf as device in vignette [2012-08-16 Thu]

## Changes in version 1.5.12

 - using knitr instead of pgfSweave and misc
   vignette updated[2012-08-13 Mon]
 - spectrum2 reporter plotting params updates [2012-08-14 Tue]
 - added reporterNames to NAMESPACE [2012-08-14 Tue]

## Changes in version 1.5.11

 - type in filterNA log messaging and also
   rounding pNA [2012-06-07 Thu]
 - typo in demo vignette - Gb instead of Mb [2012-07-15 Sun]

## Changes in version 1.5.10

 - Fixed readMSData instrumentInfo handling
   (reported by Gopuraja Dharmalingam) [2012-06-05 Tue]
 - new multiple file loading test in test_io [2012-06-05 Tue]

## Changes in version 1.5.9

 - topN now properly updates processingData [2012-06-01 Fri]
 - combineFeatures updates featureNames based on the
   groupBy argument - updated demo vignette and man
   accordinlgy [2012-06-01 Fri]
 - additional parameters were not passed when normalise
   using vsn [2012-06-01 Fri]

## Changes in version 1.5.8

 - added aa data in environment data.frame [2012-05-22 Tue]
 - fixed MSnbase:::subsetBy (used by topN) when
   ncol(object) == 1 [2012-05-25 Fri]
 - new nQuants function [2012-05-31 Thu]
 - minor vignettes updates [2012-05-31 Thu]
 - nQuants now return a matrix with col names,
   taken form sampleNames(object) [2012-05-31 Thu]

## Changes in version 1.5.7

 - changed title method to exptitle to avoid conflict/confusion
   with graphics::title and consitency with expinfo method [2012-05-15 Tue]
 - changed email to expemail [2012-05-15 Tue]
 - Added rols to Suggests [2012-05-15 Tue]

## Changes in version 1.5.6

 - new ionSource, analyser, detectorType, title accessor methods
   for MIAPE, pSet and MSnSet classes [2012-05-06 Sun]
 - updated quantify example to use data(itraqdata) instead of
   reading dummyiTRAQ.mzXML [2012-05-09 Wed]
 - initial mzTab write support [2012-05-10 Thu]
 - mzTab read support [2012-05-10 Thu]
 - updated demo and io vignettes with mzTab info [2012-05-10 Thu]

## Changes in version 1.5.5

 - added email slot in MIAPE and accessor [2012-05-04 Fri]
 - added expinfo methods to pSet and MSnSet [2012-05-04 Fri]
 - updated itraqdata [2012-05-04 Fri]
 - misc man typos fixed [2012-05-04 Fri]
 - quantify now properly propagates processingData [2012-05-04 Fri]

## Changes in version 1.5.4

 - automatically populating experiment data instrument info
   while reading data [2012-04-30 Mon]
 - msInfo fixed and exported [2012-04-30 Mon]
 - update demo vignette with 14 fractions analysis paragraph [2012-05-01 Tue]

## Changes in version 1.5.3

 - extractSpectra is now defunct [2012-04-20 Fri]
 - caching full header in level 1; this is required when and MSnExp
   instance with *many* spectra (created from many raw files) is
   quantified - calling header(object) is a too big overhead compared
   to actual reporter quantification. [2012-04-20 Fri]
 - The header() method now uses the cached dataframe if level ]= 1;
   the (unexported) .header function can be used to generate the
   dataframe using the assayData slot data. [2012-04-20 Fri]
 - Setting processingData in MSnSet initialisation. [2012-04-20 Fri]
 - Dropping index column from header. [2012-04-20 Fri]
 - new Spectrum class v0.2.0 has tic slot. [2012-04-21 Sat]
 - *tic* method (data stored as a Spectrum slot) now returns
   _total ion current_ (as commonly used) and _total ion count_ is
   obtain using *ionCount*. [2012-04-21 Sat]
 - fixed normalisation boxplot titles and other tic/ionCount
   changes in demo vignette. [2012-04-21 Sat]
 - removed qual subetting in MSnSet's "[" method [2012-04-24 Tue]

## Changes in version 1.5.2

 - transposing and MSnSet does _silently_ drop the
   protocolData now [2012-04-03 Tue]
 - fixed MSnExp pData creation for multiple files, feature names
   have a .fileNumber extension now. [2012-04-19 Thu]
 - testing for uniqueness of files (filenames) in
   readMSData [2012-04-19 Thu]
 - updated itraqdata.RData [2012-04-19 Thu]
 - added a paralle argument to quantify and using
   paralle = FALSE in vignette, to avoid duplicated display
   of the command [2012-04-19 Thu]
 - defined "reporterNames<-" generics [2012-04-19 Thu]
 - fixed warning in readMSData where all not used for
   comparison of verctors of length ] 1 [2012-04-20 Fri]

## Changes in version 1.5.1

 - updated NA warning message in readIspyData [2012-03-05 Mon]
 - fixed combineFeatures/combineMatrixFeatures for 1 sample [2012-03-20 Tue]
 - image method for MSnSet instances [2012-03-20 Tue]
 - fixed bug in plotNA (first t was wring) [2012-03-21 Wed]
 - import plot from stats4 [2012-03-21 Wed]
 - using reshape2 [2012-03-30 Fri]

## Changes in version 1.5.0

 - new devel version bump

# MSnbase 1.3

## Changes in version 1.3.15

 - new updateFeatureNames function [2012-02-17 Fri]
 - updated vignettes to illustrate vertical/horizontal combine [2012-02-17 Fri]
 - typo in normalised.Rd [2012-02-19 Sun]
 - TODO combine unit tests

## Changes in version 1.3.14

 - new is.na.MSnSet [2012-02-16 Thu]
 - updated vignette and NA related man pages with cross-links [2012-02-16 Thu]

## Changes in version 1.3.13

 - new plotNA method + doc [2012-02-15 Wed]
 - new filterNA method + doc + tests [2012-02-15 Wed]
 - added a check on 'n' in topN [2012-02-16 Thu]
 - created a .Rinstignore [2012-02-16 Thu]
 - Update package Rd [2012-02-16 Thu]

## Changes in version 1.3.12

 - new topN methods + doc + tests [2012-02-14 Tue]

## Changes in version 1.3.11

 - changed explicit close(file) in writeMgf methods to
   on.exit(close(file)) [2012-02-12 Sun]
 - typo in vignette [2012-02-13 Mon]

## Changes in version 1.3.10

 - type in writeMgfData man [2012-02-07 Tue]
 - updated TITLE in writeMgfData [2012-02-08 Wed]

## Changes in version 1.3.9

 - updated demo vignette [2012-02-03 Fri]
 - sorting numeric subsets in "[" pSet, as unsorted
   numerical indexes fails  [2012-02-03 Fri]
 - Added match.arg in combineFeatures so that
   a unique default value (the first) is used when no
   fun is specified [2012-02-03 Fri]


## Changes in version 1.3.8

 - Modified trimMz warning to report acquisition number [2012-02-01 Wed]
 - add 'experimentData(object, value) <- ' method for signature eSet
   and MIAPE [2012-02-02 Thu]
 - combine methods for MIAPE instances [2012-02-02 Thu]
 - combine methods for MSnProcess instances [2012-02-02 Thu]
 - changed qual drop warning into message in combineFeatures,
   updated test_MSnSet accordingly [2012-02-02 Thu]
 - new updateFvarLabels and updateSampleNames function [2012-02-03 Fri]
 - combine method for MSnSets [2012-02-03 Fri]
 - Updated demo vignette figure 8 [2012-02-03 Fri]

## Changes in version 1.3.7

 - Speeded up writeMgfData [2012-01-28 Sat]
 - fixes for ggplot2 0.9.0
 - added import(grid) and import(reshape) [2012-01-30 Mon]
 - importFrom(plyr, ...) instead of only llply [2012-01-31 Tue]
 - loading reshape and grid in vignette [2012-01-31 Tue]
 - fixed chunk 21 (label = quantitation-plot) [2012-01-31 Tue]

## Changes in version 1.3.6

 - Updated NoteAboutSpeedAndMemory since parallel processing
   has been added. [2011-12-18 Sun]
 - Added CITATION [2012-01-27 Fri]
 - Added information to header output: acquisition number
   and precursor intensity [2012-01-27 Fri]
 - Added a test in plot.Spectrum2 for empty dataframe [2012-01-27 Fri]
 - moved foreach, doMC to enhances [2012-01-27 Fri]

## Changes in version 1.3.5

 - added a gc() before mzR::close(msdata)... seems to help
   with Rcpp and ref classes issue. [2011-12-09 Fri]
 - added a show parameter to getCacheEnv to define .cache
   should be printed out before being returned. [2011-12-09 Fri]
 - added cache unit test [2011-12-09 Fri]
 - readMzXMLData is now defunct and remove xcms from Imports [2011-12-16 Fri]

## Changes in version 1.3.4

 - fixed bug in show MSnExp method for MS1 experiments.
   When loading MS1 spectra, cache is set to 0.
   Bug reported by Jesse Meyer. [2011-12-06 Tue]
 - fixed another bug/typo in readMSData [2011-12-06 Tue]
 - now running extractSprectum example again [2011-12-06 Tue]
 - setting default cache to 0, as cache=1 introduces
   unstabel behavious... will investigate that [2011-12-06 Tue]

## Changes in version 1.3.3

 - added parallel computation for MSnExp quantitation using
   foreach with llply(..., .parallel=TRUE) [2011-12-03 Sat]
 - TODO document above in quantify-methods.Rd
 - added foreach and doMC in Suggests [2011-12-03 Sat]
 - added Spectrum removePeaks and clean'ing in
   readMSData [2011-12-05 Mon]

## Changes in version 1.3.2

 - \dontrun{} extractSpectrum example, as this seems to be
   a major offender producing the intermittent check
   'Error in function (x)  : attempt to apply non-function'
   error [2011-11-07 Mon]
 - typo in Author@R [2011-11-14 Mon]
 - modified utils.removePeaks and utils.clean to call sapply
   instead of IRanges:sapply [2011-12-01 Thu]

## Changes in version 1.3.1

 - Herve added BioGenerics as a dependency and import
   statement in NAMESPACE  [2011-11-29 Tue]

## Changes in version 1.3.0

 - Version bump for Bioc 2.10 devel

# MSnbase 1.2

## Changes in version 1.2.0

 - Version bump for Bioc 2.9 release

# MSnbase 1.1

## Changes in version 1.1.28

 - added pgfSweave to Suggests [2011-10-23 Sun]

## Changes in version 1.1.27

 - fixed bug in readIspyData [2011-10-07 Fri]
 - removed (unexported) ratio code [2011-10-12 Wed]
 - added processing description when MSnSet['ing [2011-10-12 Wed]
 - man typos corrected [2011-10-14 Fri]
 - changed readMzXMLData to readMSData in tests [2011-10-14 Fri]
 - expecting warning for readMzXMLData in test_io [2011-10-14 Fri]

## Changes in version 1.1.26

 - extractSpectra deprecated [2011-10-05 Wed]
 - small changes in cache.R functions [2011-10-05 Wed]
 - changed [, extractSpectra and extractPrecSpectra
   to update the .cache slot [2011-10-05 Wed]
 - MSnExp show method uses .cache when level==1. [2011-10-05 Wed]
 - Finished level 1 cache implementation: leads to an
   average 11.8 times faster MSnExp show method. [2011-10-05 Wed]

## Changes in version 1.1.25

 - added .cache slot to pSet class [2011-10-03 Mon]
 - added pSet initialize method to set .cache and
   .cache$level as '0'. [2011-10-03 Mon]
 - Check that level is defined in .cache and env is locked
   in pSet validity method. [2011-10-03 Mon]
 - updated itraqdata.RData  [2011-10-03 Mon]
 - Deprecated readMzXMLData, added defunct.Rd [2011-10-03 Mon]
 - changed readMzXMLData to readMSData in man [2011-10-03 Mon]
 - updated show MSnExp for speed [2011-10-03 Mon]
 - updated read*Data is support cache = [0|1]. [2011-10-03 Mon]
 - updated precScanNum to use sapply(spectra(...), precScanNum)
   instead of unlist(eapply(assayData(...), precScanNum)) to
   preserve splectra order. [2011-10-03 Mon]
 - new cache.R file with @.cache related code [2011-10-03 Mon]

## Changes in version 1.1.24

 - in readMgfData, SCANS is not used to populate
   acquisitionNum anymore, as several scans might
   have been combined upstreams [2011-09-26 Mon]
 - added fillUp function in utils.R and exported [2011-09-27 Tue]
 - new MSnbase-io vignette [2011-09-28 Wed]

## Changes in version 1.1.23

 - fixed writeMgfData("MSnExp") [2011-09-21 Wed]
 - readMgfData now works when mz and intensty are separated
   by a '\t' (as exported by PRIDE Inspector) [2011-09-21 Wed]
 - show("MSnExp") now works without retention time [2011-09-21 Wed]
 - updated mgf2Spectrum2 to make it faster [2011-09-21 Wed]
 - fixed missing fromFile slot in data created from readMgfData
   that prevented calling header [2011-09-21 Wed]
 - readMgfData now creates a fData based on the peak header [2011-09-21 Wed]
 - modified getCurveWidth to work with centroided data [2011-09-21 Wed]
 - fixed bug getCurveWidth [2011-09-21 Wed]

## Changes in version 1.1.22

 - removed (internal) Mascot query link column in readIspyData
   to work with latest ouput version [2011-09-12 Mon]
 - removed the fillUp part in readIspyData [2011-09-14 Wed]
 - exported readIspyData [2011-09-20 Tue]
 - removed link to proteomics sig list [2011-09-20 Tue]
 - removed url in DESCRIPTION [2011-09-20 Tue]
 - Spectrum2 slot ms1scan renamed to precScanNum and populated using
   mzR::header()$precursorScanNum. Updated affected methods/functions
   and manual pages. New accessor method precScanNum is exported.
   THIS CHANGE IS NOT COMPATIBLE WITH PREVIOUSLY CREATED
   MSnSet -BJECTS! [2011-09-20 Tue]


## Changes in version 1.1.21

 - updated write.exprs to add fData columns [2011-09-08 Thu]
 - added write.exprs and readMSnSet unit test [2011-09-08 Thu]

## Changes in version 1.1.20

 - added a readMSnSet function [2011-09-08 Thu]
 - added a write.MSnSet("MSnSet") method  [2011-09-08 Thu]
 - vignette/man updates - mainly data import section documenting
   readMSData and readMgfData [2011-09-08 Thu]


 - incorporating mzR io frame work [2011-09-05 Mon]
 - use mzR's peaksCount and header generics [2011-09-05 Mon]
 - added test to check that readMzXMLData and readMSData give same output [2011-09-05 Mon]
 - added readMSData.Rd doc file [2011-09-05 Mon]
 - added Author@R field in DESCRIPTION [2011-09-07 Wed]
 - compressed/resaved (using resaveRdaFiles) itraqdata.RData file to  [2011-09-07 Wed]

## Changes in version 1.1.18

 - read support for mgf files contributed by Guangchuang Yu [2011-07-06 Wed]
 - added as.ExpressionSet.MSnSet and setAs methods and updated MSnSet doc [2011-07-09 Sat]
 - fixed read/write mgf compatibility [2011-09-01 Thu]
 - added mgf io test [2011-09-01 Thu]
 - exported and document mfg read/write support [2011-09-01 Thu]
 - added centroided parameter to rawToSpectrum[1|2] to set this directly
   at object creation in readMzXmlData [2011-09-01 Thu]
 - other minor changes in Rd files [2011-09-01 Thu]
 - added warning checks for combineFeatures when verbose=TRUE in test_MSnSet.R [2011-09-02 Fri]

## Changes in version 1.1.17

 - updated itraqdata manual [2011-06-27 Mon]
 - fixed bug in MSnSet validity method - msg was initialised
   to NULL when testing Biobase:::isValidVersion and
   Biobase::assayDataValidMembers [2011-07-02 Sat]
 - added t.MSnSet method [2011-07-02 Sat]
 - document t.MSnSet method [2011-07-05 Tue]
 - test for t-MSnSet [2011-07-05 Tue]
 - new "["-MSnSet method to properly subset qual slot [2011-07-06 Wed]
 - updated MSnSet validity method to check qual dims [2011-07-06 Wed]
 - combineFeatures now drops the spectrum-specific qual slot [2011-07-06 Wed]
 - test for "["-MSnSet [2011-07-06 Wed]

## Changes in version 1.1.16

 - added addVigs2WinMenu call in .onAttach [2011-06-26 Sun]
 - added plotMzDelta paragraph in vignette [2011-06-26 Sun]

## Changes in version 1.1.15

 - new names and description methods for ReporterIons [2011-06-16 Thu]
 - added validObject() checks [2011-06-16 Thu]
 - new removeReporters method [2011-06-16 Thu]
 - plotMzDelta is exported and documented [2011-06-17 Fri]

## Changes in version 1.1.14

 - Adding plotMzDelta QC plot (not yet exported), contributed by Guangchuang Yu [2011-06-14 Tue]
 - created a locked environment to store amino.acids dataframe [2011-06-14 Tue]
 - TODO plotMzDelta documentation and vignette section - man DONE in v 1.1.15 , vignette DONE in v 1.1.16

## Changes in version 1.1.13

 - changed pSet [[ method [2011-06-13 Mon]
 - additional updates pSet [[ method [2011-06-13 Mon]
 - added MSnExp subsetting error tests [2011-06-13 Mon]
 - created new plotting-dataframe.R file with former plotting
   utils functions, now renamed plot*.header [2011-06-13 Mon]

## Changes in version 1.1.12

 - added invisible(NULL) to all show methods [2011-06-08 Wed]
 - added centroided argument to plot.MSnExp [2011-06-09 Thu]

## Changes in version 1.1.11

 - harmonised MSnExp and Spectrum plot axes labels [2011-05-18 Wed]
 - Added plotting customisation section in vignette [2011-05-18 Wed]
 - updated signature of plot2d method to "MSnExp" only [2011-05-18 Wed]
 - added/exported plotDensity methods [2011-05-18 Wed]
 - started QC vignette section [2011-05-18 Wed]
 - added preprocSelection and preprocSelectionTable functions [2011-05-18 Wed]
 - TODO document preprocSelection[Table] in vignette
 - reduces plot2d-figure and plotDensity-figure sizes using png [2011-05-19 Thu]
 - added plotDensity doc [2011-05-19 Thu]
 - added round param to preprocSelection[Table] [2011-05-19 Thu]
 - preprocSelection[Table] documented and exported [2011-05-19 Thu]
 - fixed bug in plot.Spectrum1 [2011-05-24 Tue]
 - changed removePeaks setGeneric explicit signature [2011-05-26 Thu]
 - added MassSpectrometry biocView [2011-05-27 Fri]

## Changes in version 1.1.10

 - minor updates in demo vignette [2011-05-13 Fri]
 - added plot argument to plot methods [2011-05-16 Mon]
 - fix in makeImpuritiesMatrix [2011-05-16 Mon]
 - added meanSdPlot MSnSet method [2011-05-17 Tue]
 - minor cosmetic fix in purityCorrect error message [2011-05-17 Tue]
 - added method="sum" to Spectrum/MSnExp normalisation [2011-05-17 Tue]
 - typo in MSnSet-class.Rd corrected [2011-05-17 Tue]

## Changes in version 1.1.9

 - exporting normali[s|z]e methods for MSnSet, Spectrum and MSnExp [2011-05-12 Thu]
 - added quantile normalisation [2011-05-12 Thu]
 - added quantile.robust normalisation [2011-05-12 Thu]
 - added vsn2 normalisation [2011-05-12 Thu]
 - added normalise manual  [2011-05-12 Thu]
 - included normalisation section in vignette [2011-05-13 Fri]
 - more vignette updates [2011-05-13 Fri]
 - updated plotting methods to round MZ in title [2011-05-13 Fri]

## Changes in version 1.1.8

 - added writeMgfData method from spectra and experiment [2011-05-09 Mon]
 - added writeMgfData manual [2011-05-09 Mon]
 - added itraqdata data set and updated vignette to use it [2011-05-11 Wed]
 - updated man pages to use tiny dummyiTRAQ.mzXML [2011-05-11 Wed]
 - added makeImpuritiesMatrix function [2011-05-11 Wed]

## Changes in version 1.1.7

 - reporter ions purity correction method, man and test [2011-05-09 Mon]
 - updated vignette [2011-05-09 Mon]

## Changes in version 1.1.6

 - added incomplete dissociation and spectral counting sections
   to demo vignette [2011-05-06 Fri]
 - added bioc-sig-proteomics link to foreword [2011-05-06 Fri]
 - type in foreword [2011-05-06 Fri]

## Changes in version 1.1.5

 - added combineFeatures example in demo vignette [2011-05-05 Thu]

## Changes in version 1.1.4

 - readIspyData updated to return updated factors [2011-05-03 Tue]
 - added unexported/undocumented combineFeatures function for MSnSets [2011-05-03 Tue]
 - added basic tests for combineFeatures [2011-05-03 Tue]
 - added combineFeatures manual [2011-05-04 Wed]
 - exportig combineFeatures  [2011-05-04 Wed]

## Changes in version 1.1.3

 - as.data.frame.Spectrum columns now are (first) mz and (second) intensity [2011-04-27 Wed]
 - exporting as.data.frame.Spectrum and coerce [2011-04-27 Wed]

## Changes in version 1.1.2

 - Simplified quantify generic signature - now only object argument [2011-04-19 Tue]
 - Added strict parameter to quantify method, man updated, added relevant test
 - Added illustrative plot for quantitation methods in MSnbase-demo vignette [2011-04-19 Tue]
 - Added illustrative plot for data pre-processing (removePeaks and clean)
   in MSnbase-demo vignette [2011-04-20 Wed]
 - No warnings are issued anymore when peaks expands outside of mz(reporters) +/- width(reporters).
   See ?quantify on how to check this manually. [2011-04-19 Tue]
 - No warnings are issued anymore when reporter peaks are missing.
   See ?quantify on how to check this manually. [2011-04-20 Wed]
 - pSet validity warns if length(unique(msLevel(object))) ] 1, rather than != 1.
   The latter triggered a warning for a new("MSnExp"). [2011-04-20 Wed]

## Changes in version 1.1.1

 - added setAs data.frame and as.data.frame methods for
   Spectrum objects [2011-03-29 Tue]
 - support for uncentroided MS2 spectra plots [2011-03-31 Thu] [2011-04-02 Sat]
 - support for uncentroided MS1 spectra plots [2011-04-02 Sat]
 - minor modification to readIspyData [2011-04-04 Mon]
 - removed centroided slot from MSnProcess and added to
   individial Spectrum instances. Relevant for Velos
   HCD (profile)/CID (uncentroided) data [2011-04-04 Mon]
 - modified readMzXmlData accordingly [2011-04-04 Mon]
 - added validObject(new(...)) tests for each class [2011-04-04 Mon]
 - added centroided[<-] methods to Spectrum and pSet [2011-04-04 Mon]
 - added 'keepAll' parameter to readIspyData [2011-04-11 Mon]

# MSnbase 0.99

## Changes in version 0.99.4

 - removed pgfSweave from Suggests field [2011-04-04 Mon]

## Changes in version 0.99.3

 - Using Sweave rather than pgfSweave to build on lamb1. [2011-03-28 Mon]

## Changes in version 0.99.2

 - updated references in package Rd [2011-03-24 Thu]
 - added readIspyData.Rd  [2011-03-24 Thu]
 - removed old readMzXMLData function [2011-03-24 Thu]
 - cleaning up code [2011-03-24 Thu]
 - fixed MSnSet initialize: setting experimentData as MIAPE [2011-03-25 Fri]
 - fixed MSnSet initialize: identical featureNames in assayData and
   featureData for empty MSnExp istances [2011-03-25 Fri]
 - better test organisation and added object validity tests [2011-03-25 Fri]
 - updated readIspyData function [2011-03-25 Fri]

## Changes in version 0.99.1

 - fake vignettes are copied back in inst/doc/. to
   make sur tar tarball can be extracted and build
   again [2011-03-24 Thu]
 - Making clean after make all [2011-03-24 Thu]

## Changes in version 0.99.0

 - MSnProcess@process get's one line when subsetting and MSnExp object [2011-03-21 Mon]
 - Updated vignette Makefile [2011-03-22 Tue]
 - Added fake vignettes [2011-03-22 Tue]
 - Version set to 0.99.0 [2011-03-22 Tue]
 - Minor changes in quantify-method.Rd [2011-03-22 Tue]
 - added .Rbuildignore to ignore 'sweave-cache.*' [2011-03-22 Tue]
 - added R/zzz.R with start-up message [2011-03-23 Wed]
 - Added width generic [2011-03-23 Wed]
 - IRanges in now imported [2011-03-23 Wed]

# MSnbase 0.2

## Changes in version 0.2.0 [2011-03-17 Thu]

 - MSnSet now extends eSet, reimplementing the ExpressionSet class
   with exception of the experimentData slot, that must now be a MIAPE
   instance.
 - updated MIAPE substantially and added accessor methods.
 - updated docs and vignettes to reflect above changes.

# MSnbase 0.1

## Changes in version 0.1.5

 - added zoo to Suggest, as zoo::rollapply is used in vignette.
 - changes rollapply call to zoo:::rollapply.zoo, to make it
   work with zoo_1.6-4, rather that depending on zoo_1.7-0, which
   is not yet on CRAN (only r-forge).
 - added vns to Suggests, as it is used in the demo vignette.

## Changes in version 0.1.4

 - added pgfSweave to Suggests, as is called in vignette Makefile

## Changes in version 0.1.3

 - Updates to ReporterIons man page
 - Moved some Depends to Imports and updated NAMESPACE
 - corrected a IRanges:::sapply ot IRanges::sapply in utils.R
 - added inst/tests/test_MSnProcess.R

## Changes in version 0.1.2

 - MSnProcess@MSnbaseVersion is not set in the initialize method instead to the prototype
 - changed R CMD to $(R_HOME)/bin/R CMD in inst/doc/Makefile
 - removed first and last (was 5th) MS1 spectra and associated MS2 spectra
   from dummy file to reduce size below 2MB.
 - Updated tests to reflect new dummy data set
 - Edition of dummy brakes readMzXMLData with msLevel=1 -- comment related tests
 - Updated MSnbase-demo vignette and man to use good spectra un plots (i.e "X43" instead of "X64")

## Changes in version 0.1.1

 - Added LazyLoad: yes in DESCRIPTION and modified docs accordingly
 - Added NEWS file

## Changes in version 0.1.0

 - Added MSnbase-demo vignette
 - Added MSnbase-development vignette
