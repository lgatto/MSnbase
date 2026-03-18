# The "MIAPE" Class for Storing Proteomics Experiment Information

The Minimum Information About a Proteomics Experiment. The current
implementation is based on the MIAPE-MS 2.4 document.

## Slots

- `title`::

  Object of class `character` containing a single-sentence experiment
  title.

- `abstract`::

  Object of class `character` containing an abstract describing the
  experiment.

- `url`::

  Object of class `character` containing a URL for the experiment.

- `pubMedIds`::

  Object of class `character` listing strings of PubMed identifiers of
  papers relevant to the dataset.

- `samples`::

  Object of class `list` containing information about the samples.

- `preprocessing`::

  Object of class `list` containing information about the pre-processing
  steps used on the raw data from this experiment.

- `other`::

  Object of class `list` containing other information for which none of
  the above slots applies.

- `dateStamp`::

  Object of class `character`, giving the date on which the work
  described was initiated; given in the standard 'YYYY-MM-DD' format
  (with hyphens).

- `name`::

  Object of class `character` containing the name of the (stable)
  primary contact person for this data set; this could be the
  experimenter, lab head, line manager, ...

- `lab`::

  Object of class `character` containing the laboratory where the
  experiment was conducted.

- `contact`::

  Object of class `character` containing contact information for lab
  and/or experimenter.

- `email`::

  Object of class `character` containing tmail contact information for
  the primary contact person (see `name` above).

- `instrumentModel`::

  Object of class `character` indicating the model of the mass
  spectrometer used to generate the data.

- `instrumentManufacturer`::

  Object of class `character` indicating the manufacturing company of
  the mass spectrometer.

- `instrumentCustomisations`::

  Object of class `character` describing any significant (i.e. affecting
  behaviour) deviations from the manufacturer's specification for the
  mass spectrometer.

- `softwareName`::

  Object of class `character` with the instrument management and data
  analysis package(s) name(s).

- `softwareVersion`::

  Object of class `character` with the instrument management and data
  analysis package(s) version(s).

- `switchingCriteria`::

  Object of class `character` describing the list of conditions that
  cause the switch from survey or zoom mode (MS1) to or tandem mode (MSn
  where n \> 1); e.g. 'parent ion” mass lists, neutral loss criteria and
  so on \[applied for tandem MS only\].

- `isolationWidth`::

  Object of class `numeric` describing, for tandem instruments, the
  total width (i.e. not half for plus-or-minus) of the gate applied
  around a selected precursor ion m/z, provided for all levels or by MS
  level.

- `parameterFile`::

  Object of class `character` giving the location and name under which
  the mass spectrometer's parameter settings file for the run is stored,
  if available. Ideally this should be a URI+filename, or most
  preferably an LSID, where feasible.

- `ionSource`::

  Object of class `character` describing the ion source (ESI, MALDI,
  ...).

- `ionSourceDetails`::

  Object of class `character` describing the relevant details about the
  ion source. See MIAPE-MI docuement for more details.

- `analyser`::

  Object of class `character` describing the analyzer type (Quadrupole,
  time-of-flight, ion trap, ...).

- `analyserDetails`::

  Object of class `character` describing the relevant details about the
  analyzer. See MIAPE-MI document for more details.

- `collisionGas`::

  Object of class `character` describing the composition of the gas used
  to fragment ions in the collision cell.

- `collisionPressure`::

  Object of class `numeric` providing the pressure (in bars) of the
  collision gas.

- `collisionEnergy`::

  Object of class `character` specifying for the process of imparting a
  particular impetus to ions with a given m/z value, as they travel into
  the collision cell for fragmentation. This could be a global figure
  (e.g. for tandem TOF's), or a complex function; for example a gradient
  (stepped or continuous) of m/z values (for quads) or activation
  frequencies (for traps) with associated collision energies (given in
  eV). Note that collision energies are also provided for individual
  `"`[`Spectrum2`](https://lgatto.github.io/MSnbase/reference/Spectrum2-class.md)`"`
  instances, and is the preferred way of accessing this data.

- `detectorType`::

  Object of class `character` describing the type of detector used in
  the machine (microchannel plate, channeltron, ...).

- `detectorSensitivity`::

  Object of class `character` giving and appropriate measure of the
  sensitivity of the described detector (e.g. applied voltage).

## Methods

The following methods as in `"MIAME"`:

- `abstract(MIAPE)`::

  An accessor function for `abstract`.

- `expinfo(MIAPE)`::

  An accessor function for `name`, `lab`, `contact`, `title`, and `url`.

- `notes(MIAPE), notes(MIAPE) <- value`::

  Accessor functions for `other`. `notes(MIAME) <- character` *appends*
  character to notes; use `notes(MIAPE) <- list` to replace the notes
  entirely.

- `otherInfo(MIAPE)`::

  An accessor function for `other`.

- `preproc(MIAPE)`::

  An accessor function for `preprocessing`.

- `pubMedIds(MIAPE), pubMedIds(MIAME) <- value`::

  Accessor function for `pubMedIds`.

- `expemail(MIAPE)`::

  An accessor function for `email` slot.

- `exptitle(MIAPE)`::

  An accessor function for `title` slot.

- `analyzer(MIAPE)`::

  An accessor function for `analyser` slot. `analyser(MIAPE)` is also
  available.

- `analyzerDetails(MIAPE)`::

  An accessor function for `analyserDetails` slot. `analyserDetails` is
  also available.

- `detectorType(MIAPE)`::

  An accessor function for `detectorType` slot.

- `ionSource(MIAPE)`::

  An accessor function for `ionSource` slot.

- `ionSourceDetails(MIAPE)`::

  An accessor function for `ionSourceDetails` slot.

- `instrumentModel(MIAPE)`::

  An accessor function for `instrumentModel` slot.

- `instrumentManufacturer(MIAPE)`::

  An accessor function for `instrumentManufacturer` slot.

- `instrumentCustomisations(MIAPE)`::

  An accessor function for `instrumentCustomisations` slot.

- `as(,"MIAME")`::

  Coerce the object from `MIAPE` to `MIAME` class. Used when converting
  an `MSnSet` into an `ExpressionSet`.

MIAPE-specific methods, including MIAPE-MS meta-data:

- `show(MIAPE)`::

  Displays the experiment data.

- `msInfo(MIAPE)`::

  Displays 'MIAPE-MS' information.

## Extends

Class `"MIAxE"`, directly. Class `"Versioned"`, by class "MIAxE",
distance 2.

## References

About MIAPE: <http://www.psidev.info/index.php?q=node/91>, and
references therein, especially 'Guidelines for reporting the use of mass
spectrometry in proteomics', Nature Biotechnology 26, 860-861 (2008).

## Author

Laurent Gatto
