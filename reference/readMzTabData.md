# Read an 'mzTab' file

This function can be used to create an
`"`[`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.md)`"`
by reading and parsing an `mzTab` file. The metadata section is always
used to populate the `MSnSet`'s `experimentData()@other$mzTab` slot.

## Usage

``` r
readMzTabData(
  file,
  what = c("PRT", "PEP", "PSM"),
  version = c("1.0", "0.9"),
  verbose = isMSnbaseVerbose()
)
```

## Arguments

- file:

  A `character` with the `mzTab` file to be read in.

- what:

  One of `"PRT"`, `"PEP"` or `"PSM"`, defining which of protein, peptide
  PSMs section should be returned as an `MSnSet`.

- version:

  A `character` defining the format specification version of the mzTab
  file. Default is `"1.0"`. Version `"0.9"` is available of backwards
  compatibility. See
  [`readMzTabData_v0.9`](https://lgatto.github.io/MSnbase/reference/readMzTabData_v0.9.md)
  for details.

- verbose:

  Produce verbose output.

## Value

An instance of class `MSnSet`.

## See also

See [`MzTab`](https://lgatto.github.io/MSnbase/reference/MzTab-class.md)
and
[`MSnSetList`](https://lgatto.github.io/MSnbase/reference/MSnSetList-class.md)
for details about the inners of `readMzTabData`.

## Author

Laurent Gatto

## Examples

``` r
testfile <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/1_0-Proteomics-Release/PRIDE_Exp_Complete_Ac_16649.xml-mztab.txt"

prot <- readMzTabData(testfile, "PRT")

prot
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 1249 features, 4 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData: none
#> featureData
#>   featureNames: X223462890 X19855078 ... X26329627 (1249 total)
#>   fvarLabels: accession description ... protein_coverage (15 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#>   pubMedIds: pubmed:21398567 
#> Annotation:  
#> - - - Processing information - - -
#>  MSnbase version: 2.37.3 

head(fData(prot))
#>            accession
#> X223462890 223462890
#> X19855078   19855078
#> X21450277   21450277
#> X6978545     6978545
#> X51315739   51315739
#> X117938332 117938332
#>                                                                                                                                                                                                           description
#> X223462890                                                                                                                                                                               Spna2 protein [Mus musculus]
#> X19855078  RecName: Full=Sodium/potassium-transporting ATPase subunit alpha-3; Short=Na(+)/K(+) ATPase alpha-3 subunit; AltName: Full=Na(+)/K(+) ATPase alpha(III) subunit; AltName: Full=Sodium pump subunit alpha-3
#> X21450277                                                                                                                               sodium/potassium-transporting ATPase subunit alpha-1 precursor [Mus musculus]
#> X6978545                                                                                                                           sodium/potassium-transporting ATPase subunit alpha-2 precursor [Rattus norvegicus]
#> X51315739                                                                                                                                                                               RecName: Full=Protein bassoon
#> X117938332                                                                                                                                                      spectrin beta chain, brain 1 isoform 1 [Mus musculus]
#>            taxid              species       database database_version
#> X223462890 10090 Mus musculus (Mouse) NCBInr_2010_10  nr_101020.fasta
#> X19855078  10090 Mus musculus (Mouse) NCBInr_2010_10  nr_101020.fasta
#> X21450277  10090 Mus musculus (Mouse) NCBInr_2010_10  nr_101020.fasta
#> X6978545   10090 Mus musculus (Mouse) NCBInr_2010_10  nr_101020.fasta
#> X51315739  10090 Mus musculus (Mouse) NCBInr_2010_10  nr_101020.fasta
#> X117938332 10090 Mus musculus (Mouse) NCBInr_2010_10  nr_101020.fasta
#>                         search_engine best_search_engine_score[1]
#> X223462890 [MS, MS:1001207, Mascot, ]                     6539.67
#> X19855078  [MS, MS:1001207, Mascot, ]                     6331.91
#> X21450277  [MS, MS:1001207, Mascot, ]                     4577.11
#> X6978545   [MS, MS:1001207, Mascot, ]                     4342.81
#> X51315739  [MS, MS:1001207, Mascot, ]                     4177.55
#> X117938332 [MS, MS:1001207, Mascot, ]                     4001.66
#>            search_engine_score[1]_ms_run[1] num_psms_ms_run[1]
#> X223462890                          6539.67                157
#> X19855078                           6331.91                144
#> X21450277                           4577.11                112
#> X6978545                            4342.81                108
#> X51315739                           4177.55                100
#> X117938332                          4001.66                109
#>            num_peptides_distinct_ms_run[1] num_peptides_unique_ms_run[1]
#> X223462890                              92                            NA
#> X19855078                               49                            NA
#> X21450277                               39                            NA
#> X6978545                                42                            NA
#> X51315739                               59                            NA
#> X117938332                              72                            NA
#>            ambiguity_members
#> X223462890                NA
#> X19855078                 NA
#> X21450277                 NA
#> X6978545                  NA
#> X51315739                 NA
#> X117938332                NA
#>                                                                                 modifications
#> X223462890                                                                               <NA>
#> X19855078  32-MOD:00425,525-MOD:00425,606-MOD:00425,725-MOD:00425,739-MOD:00425,940-MOD:00425
#> X21450277                              42-MOD:00425,616-MOD:00425,749-MOD:00425,950-MOD:00425
#> X6978545                               40-MOD:00425,613-MOD:00425,746-MOD:00425,947-MOD:00425
#> X51315739                                                                                <NA>
#> X117938332                                                                               <NA>
#>            protein_coverage
#> X223462890                0
#> X19855078                 0
#> X21450277                 0
#> X6978545                  0
#> X51315739                 0
#> X117938332                0

head(exprs(prot))
#>            protein_abundance_assay[1] protein_abundance_assay[2]
#> X223462890                          1                      0.853
#> X19855078                          NA                         NA
#> X21450277                           1                      0.776
#> X6978545                            1                      0.784
#> X51315739                          NA                         NA
#> X117938332                          1                      0.865
#>            protein_abundance_assay[3] protein_abundance_assay[4]
#> X223462890                      0.864                      0.791
#> X19855078                          NA                         NA
#> X21450277                       0.819                      0.687
#> X6978545                        0.848                      0.693
#> X51315739                          NA                         NA
#> X117938332                      0.861                      0.795

psms <- readMzTabData(testfile, "PSM")

psms
#> MSnSet (storageMode: lockedEnvironment)
#> assayData: 8761 features, 0 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData: none
#> featureData
#>   featureNames: X1661 X2280 ... X20346 (8761 total)
#>   fvarLabels: sequence PSM_ID ... end (18 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#>   pubMedIds: pubmed:21398567 
#> Annotation:  
#> - - - Processing information - - -
#>  MSnbase version: 2.37.3 

head(fData(psms))
#>       sequence PSM_ID accession unique       database database_version
#> X1661   QQVLDR   1661 223462890     NA NCBInr_2010_10  nr_101020.fasta
#> X2280   LVQYLR   2280 223462890     NA NCBInr_2010_10  nr_101020.fasta
#> X2281   LVQYLR   2281 223462890     NA NCBInr_2010_10  nr_101020.fasta
#> X2537   LQQLFR   2537 223462890     NA NCBInr_2010_10  nr_101020.fasta
#> X2809 EAGSVSLR   2809 223462890     NA NCBInr_2010_10  nr_101020.fasta
#> X5465 LSILSEER   5465 223462890     NA NCBInr_2010_10  nr_101020.fasta
#>                    search_engine search_engine_score[1] modifications
#> X1661 [MS, MS:1001207, Mascot, ]                  37.76   0-MOD:01499
#> X2280 [MS, MS:1001207, Mascot, ]                  44.64   0-MOD:01499
#> X2281 [MS, MS:1001207, Mascot, ]                  44.76   0-MOD:01499
#> X2537 [MS, MS:1001207, Mascot, ]                  45.41   0-MOD:01499
#> X2809 [MS, MS:1001207, Mascot, ]                  55.05   0-MOD:01499
#> X5465 [MS, MS:1001207, Mascot, ]                  39.82   0-MOD:01499
#>       retention_time charge exp_mass_to_charge calc_mass_to_charge
#> X1661             NA      1           902.4821            902.5181
#> X2280             NA      1           935.5775            935.5800
#> X2281             NA      1           935.5833            935.5800
#> X2537             NA      1           948.5956            948.5753
#> X2809             NA      1           962.5098            962.5393
#> X5465             NA      1          1090.6232           1090.6230
#>                   spectra_ref pre post start  end
#> X1661 ms_run[1]:spectrum=1661   R    Y    20   25
#> X2280 ms_run[1]:spectrum=2280   K    E   151  156
#> X2281 ms_run[1]:spectrum=2281   K    E   151  156
#> X2537 ms_run[1]:spectrum=2537   R    D   786  791
#> X2809 ms_run[1]:spectrum=2809   K    M  1058 1065
#> X5465 ms_run[1]:spectrum=5465   K    T   442  449
```
