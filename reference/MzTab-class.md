# Parse `MzTab` files

The `MzTab` class stores the output of a basic parsing of a `mzTab`
file. It contain the metadata (a `list`), comments (a `character`
vector), and the at least of of the following data types: proteins,
peptides, PSMs and small molecules (as `data.frames`).

At this stage, the metadata and data are only minimally parsed. The
column names are kept as they are defined in the original files and are
thus not all going to be valid colnames. To access them using the dollar
operator, use backticks. More specific data extraction and preparation
are delegated to more specialised functions, such as the
`as(., to = "MSnSetList")` and `readMzTabData` for proteomics data.

Note that no attempts are made to verify the validitiy of the mzTab
file.

## Objects from the Class

Objects can be created by calls the the constructor `MzTab` that takes a
single `mzTab` file as input.

The objects can subsequently be coerced to
[MSnSetList](https://lgatto.github.io/MSnbase/reference/MSnSetList-class.md)
instances with `as(object, "MSnSetList")`. The resulting `MSnSetList`
contains possibly empty `MSnSet` instances for proteins, peptide and
PSMs, respectively named `"Proteins"`, `"Peptides"` and `"PSMs"`.

The assaydata slots of the two former are populated with the
`protein_abundance_assay[1-n]` and `peptide_abundance_assay[1-n]`
columns in the `mzTab` file. No abundance values are defined for the
latter. The respective feature names correspond to protein accessions,
peptide sequences and PSM identifiers, possibly made unique as by
appending sequence numbers to duplicates.

## Slots

- `Metadata`::

  Object of class `"list"` storing the metadata section.

- `Filename`::

  Object of class `"character"` storing the orginal file name.

- `Proteins`::

  Object of class `"data.frame"` storing the protein data.

- `Peptides`::

  Object of class `"data.frame"` storing the peptide data.

- `PSMs`::

  Object of class `"data.frame"` storing the PSM data.

- `SmallMolecules`::

  Object of class `"data.frame"` storing the small molecules data.

- `MoleculeFeatures`::

  Object of class `"data.frame"` storing the molecule features.

- `MoleculeEvidence`::

  Object of class `"data.frame"` storing the molecule evidence.

- `Comments`::

  Object of class `"character"` storing the comments that were present
  in the file.

## Accessors

- metadata:

  `signature(x = "MzTab")`: returns the meta data `list`.

- mzTabMode:

  `signature(x = "MzTab")`: returns the mode (complete or summary) of
  the `mzTab` data. A shortcut for `` metadata(x)$`mzTab-mode` ``.

- mzTabType:

  `signature(x = "MzTab")`: returns the type (quantification or
  identification) of the `mzTab` data. A shortcut for
  `` metadata(x)$`mzTab-type` ``.

- fileName:

  `signature(object = "MzTab")`: returns the file name of the original
  `mzTab` file.

- peptides:

  `signature(object = "MzTab")`: returns the peptide `data.frame`.

- proteins:

  `signature(object = "MzTab")`: returns the proteins `data.frame`.

- psms:

  `signature(object = "MzTab")`: returns the PSMs `data.frame`.

- smallMolecules:

  `signature(object = "MzTab")`: returns the small molecules (SML)
  `data.frame`.

- moleculeFeatures:

  `signature(object = "MzTab")`: returns the small molecules features
  (SMF) `data.frame`.

- moleculeEvidence:

  `signature(object = "MzTab")`: returns the small molecule
  identification evidence (SME) `data.frame`.

- comments:

  `signature(object = "MzTab")`: returns the comments.

## References

The mzTab format is a light-weight, tab-delimited file format for
proteomics data. Version mzTab 1.0 is aimed at proteomics, mzTab-M 2.0
was adapted to metabolomics. See https://github.com/HUPO-PSI/mzTab for
details and specifications.

Griss J, Jones AR, Sachsenberg T, Walzer M, Gatto L, Hartler J,
Thallinger GG, Salek RM, Steinbeck C, Neuhauser N, Cox J, Neumann S, Fan
J, Reisinger F, Xu QW, Del Toro N, Perez-Riverol Y, Ghali F, Bandeira N,
Xenarios I, Kohlbacher O, Vizcaino JA, Hermjakob H. The mzTab data
exchange format: communicating mass-spectrometry-based proteomics and
metabolomics experimental results to a wider audience. Mol Cell
Proteomics. 2014 Oct;13(10):2765-75. doi: 10.1074/mcp.O113.036681. Epub
2014 Jun 30. PubMed PMID: 24980485; PubMed Central PMCID: PMC4189001.

Hoffmann N, Rein J, Sachsenberg T, et al. mzTab-M: A Data Standard for
Sharing Quantitative Results in Mass Spectrometry Metabolomics. Anal
Chem. 2019;91(5):3302‐3310. doi:10.1021/acs.analchem.8b04310 PubMed
PMID: 30688441; PubMed Central PMCID: PMC6660005.

## Author

Laurent Gatto, with contributions from Richard Cotton (see
<https://github.com/lgatto/MSnbase/issues/41>) and Steffen Neuman (see
https://github.com/lgatto/MSnbase/pull/500).

## Examples

``` r
## Test files from the mzTab developement repository
fls <- c("Cytidine.mzTab", "MTBLS2.mztab", 
         "PRIDE_Exp_Complete_Ac_1643.xml-mztab.txt", 
         "PRIDE_Exp_Complete_Ac_16649.xml-mztab.txt", 
         "SILAC_CQI.mzTab", "SILAC_SQ.mzTab", 
         "iTRAQ_CQI.mzTab", "iTRAQ_SQI.mzTab", 
         "labelfree_CQI.mzTab", "labelfree_SQI.mzTab", 
         "lipidomics-HFD-LD-study-PL-DG-SM.mzTab", 
         "lipidomics-HFD-LD-study-TG.mzTab")

baseUrl <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/1_0-Proteomics-Release/"

## a list of mzTab objects
mzt <- sapply(file.path(baseUrl, fls), MzTab)
stopifnot(length(mzt) == length(fls))
mzt[[4]]
#> Object of class "MzTab".
#>  Description: date of export: Mon Jun 16 10:57:56 BST 2014
#>  Mode: Summary 
#>  Type: Quantification 
#>  Available data: Proteins  PSMs  

dim(proteins(mzt[[4]]))
#> [1] 1249   19
dim(psms(mzt[[4]]))
#> [1] 8761   18

prots4 <- proteins(mzt[[4]])
class(prots4)
#> [1] "data.frame"
prots4[1:5, 1:4]
#>   accession
#> 1 223462890
#> 2  19855078
#> 3  21450277
#> 4   6978545
#> 5  51315739
#>                                                                                                                                                                                                  description
#> 1                                                                                                                                                                               Spna2 protein [Mus musculus]
#> 2 RecName: Full=Sodium/potassium-transporting ATPase subunit alpha-3; Short=Na(+)/K(+) ATPase alpha-3 subunit; AltName: Full=Na(+)/K(+) ATPase alpha(III) subunit; AltName: Full=Sodium pump subunit alpha-3
#> 3                                                                                                                              sodium/potassium-transporting ATPase subunit alpha-1 precursor [Mus musculus]
#> 4                                                                                                                         sodium/potassium-transporting ATPase subunit alpha-2 precursor [Rattus norvegicus]
#> 5                                                                                                                                                                              RecName: Full=Protein bassoon
#>   taxid              species
#> 1 10090 Mus musculus (Mouse)
#> 2 10090 Mus musculus (Mouse)
#> 3 10090 Mus musculus (Mouse)
#> 4 10090 Mus musculus (Mouse)
#> 5 10090 Mus musculus (Mouse)
```
