# WEIN - Web-based Engine for Interactive Next-generation sequencing analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

WEIN is an R/Bioconductor package containing a Shiny application for analyzing RNA-Seq data in the context of differential expression. This enables an interactive and at the same time reproducible analysis, keeping the functionality accessible, and yet providing a comprehensive selection of graphs and tables to mine the dataset at hand.

## Installation

WEIN can be easily installed using `BiocManager::install()`:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
    
BiocManager::install("WEIN")
library(WEIN)
WEIN()
```

Or from GitHub:

```r
install.packages("devtools")
devtools::install_github("baj12/WEIN", dependencies = T)
library(WEIN)
WEIN()
```

## Quick start

This command loads the WEIN package:

```r
library("WEIN")
```

The main parameters for WEIN are:

- `dds_obj` - a `DESeqDataSet` object. If not provided, then a `countmatrix` and a `expdesign` need to be provided. If none of the above is provided, it is possible to upload the data during the execution of the Shiny App
- `res_obj` -  a `DESeqResults` object. If not provided, it can be computed during the execution of the application
- `annotation_obj` - a `data.frame` object, with row.names as gene identifiers (e.g. ENSEMBL ids) and a column, `gene_name`, containing e.g. HGNC-based gene symbols. If not provided, it can be constructed during the execution via the `org.eg.XX.db` packages
- `countmatrix` - a count matrix, with genes as rows and samples as columns. If not provided, it is possible to upload the data during the execution of the Shiny App
- `expdesign` -a `data.frame` containing the info on the experimental covariates of each sample. If not provided, it is possible to upload the data during the execution of the Shiny App

The WEIN app can be launched in different modes:

- `WEIN(dds_obj = dds, res_obj = res, annotation_obj = anno)`, where the objects are precomputed in the current session and provided as parameters
- `WEIN(dds_obj = dds)`, as in the command above, but where the result object is assembled at runtime
- `WEIN(countmatrix = countmatrix, expdesign = expdesign)`, where instead of passing the defined `DESeqDataSet` object, its components are given, namely the count matrix (e.g. generated after a run of featureCounts or HTSeq-count) and a data frame with the experimental covariates. The design formula can be constructed interactively at runtime
- `WEIN()`, where the count matrix and experimental design can simply be uploaded at runtime, where all the derived objects can be extracted and computed live. These files have to be formatted as tabular text files, and a function in the package tries to guess the separator, based on heuristics of occurrencies per line of commonly used characters

## Debugging Options

The application includes debugging output that can be enabled/disabled using the global option `wein.debug`. By default, debugging output is disabled to reduce console clutter.

To enable debugging output, set the option before launching the app:
```r
options(wein.debug = TRUE)
WEIN()
```

When enabled, detailed timing information and debug messages will be printed to the console during application execution.

## Large Files Notice

This repository excludes large data files that exceed GitHub's file size limits. These files include:

- `idealState_*.RData` (large session state files)
- `shiny_bookmarks/*` (bookmark files that can become very large)
- `test.Rdata` and `testData.RData` (large test data files)
- Other large `.RData` and `.rds` files

If you need these files for development or testing purposes, please contact the repository maintainer.

## Development

To set up the development environment:

1. Clone the repository
2. Install dependencies with `BiocManager::install()`
3. Load the package with `library(WEIN)`
4. Run the app with `WEIN()`

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bug reports, feature requests, or general feedback.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For additional details regarding the functions of WEIN, please consult the documentation or write an email to bernd.jagla@pasteur.fr.

### Bug reports/Issues/New features

Please use https://github.com/baj12/WEIN/issues for reporting bugs, issues or for suggesting new features to be implemented.
