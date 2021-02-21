
# scPseudoBulk.avg

<!-- badges: start -->
<!-- badges: end -->

The goal of scPseudoBulk.avg is to prduce pseudo-bulk samples of the clusters identified from scRNA-Seq data
by first filtering (to remove technical artefacts) and averaging expression values of cells for each cluster. 
This would help improve the signal to noise ratio.

## Installation

You can install the released version of scPseudoBulk.avg from [github](https://github.com/anirudhpatir) with:

``` r
devtools::install_github("scPseudoBulk.avg")
```

## Usage

``` r
library(scPseudoBulk.avg)

scbulk = scPseudoBulk(norm_matrix, clusters)

```

# Examples

See [vignette]() for a demo.

