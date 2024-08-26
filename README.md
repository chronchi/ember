# Molecula EMBeddeR: ember

<!-- badges: start -->
<!-- badges: end -->

A framework for precision medicine in breast cancer. Calculates molecular subtype and provides pathway scores for samples individually.

## Installation

To install the package run R and the following lines:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("chronchi/ember")
```

## Basic usage 

Check the tutorial page at [chronchi.github.io/ember/](https://chronchi.github.io/ember/ "Tutorial for ember")

## Citation

If you used the package (and enjoyed it!), use the following for citation:

> Ronchi, C., Haider, S. & Brisken, C. EMBER creates a unified space for independent breast cancer transcriptomic datasets enabling precision oncology. npj Breast Cancer 10, 56 (2024). https://doi.org/10.1038/s41523-024-00665-z

For the bibtex: 

```
@article{Ronchi2024,
  title = {EMBER creates a unified space for independent breast cancer transcriptomic datasets enabling precision oncology},
  volume = {10},
  ISSN = {2374-4677},
  url = {http://dx.doi.org/10.1038/s41523-024-00665-z},
  DOI = {10.1038/s41523-024-00665-z},
  number = {1},
  journal = {npj Breast Cancer},
  publisher = {Springer Science and Business Media LLC},
  author = {Ronchi, Carlos and Haider, Syed and Brisken, Cathrin},
  year = {2024},
  month = jul 
}
```
