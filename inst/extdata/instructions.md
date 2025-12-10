*This information is also contained in the `WEIN` package documentation.*

`WEIN` (Web-based Engine for Interactive Next-generation sequencing analysis) is a Bioconductor package containing a Shiny application for interactively analyzing RNA-seq expression data, by interactive exploration of the results of a Differential Expression analysis.

## Getting started

*WEIN* is an R package developed by Bernd Jagla at the Institut Pasteur, based on the original *ideal* package developed by Federico Marini. To install the package, start R and enter:

```r
# Installation instructions will be added here once the package is ready for distribution
```

Once *WEIN* is installed, it can be loaded by the following command.

```r
library("WEIN")
```

## Introduction

*WEIN* (Web-based Engine for Interactive Next-generation sequencing analysis) is a Bioconductor package containing a Shiny application for analyzing RNA-Seq data in the context of differential expression. This enables an interactive and at the same time analysis, keeping the functionality accessible, and yet providing a comprehensive selection of graphs and tables to mine the dataset at hand.

*WEIN* is an R package which fully leverages the infrastructure of the Bioconductor project in order to deliver an interactive yet reproducible analysis for the detection of differentially expressed genes in RNA-Seq datasets. Graphs, tables, and interactive HTML reports can be readily exported and shared across collaborators. The dynamic user interface displays a broad level of content and information, subdivided by thematic tasks. All in all, it aims to enforce a proper analysis, by reaching out both life scientists and experienced bioinformaticians, and also fosters the communication between the two sides, offering robust statistical methods and high standard of accessible documentation.

It is structured in a similar way to the *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)*, also designed as an interactive companion tool for RNA-seq analysis focused rather on the exploratory data analysis e.g. using principal components analysis as a main tool.

The interactive/reactive design of the app, with a dynamically generated user interface makes it easy and immediate to apply the gold standard methods (in the current implementation, based on *[DESeq2](http://bioconductor.org/packages/DESeq2)*) in a way that is information-rich and accessible also to the bench biologist, while also providing additional insight also for the experienced data analyst. Reproducibility is supported via state saving and automated report generation.

### Citation info

If you use *WEIN* for your analysis, please cite the original *ideal* package as well as acknowledge the development of *WEIN*:

```r
citation("ideal") # Original package by Federico Marini
```

## Using the application

There are different ways to use `WEIN` for interactive differential expression analysis.

### Launching `WEIN` locally

First load the library

```r
library("WEIN")
```

and then launch the app with the `WEIN` function. This takes the following essential parameters as input:

- `dds_obj` - a `DESeqDataSet` object. If not provided, then a `countmatrix` and a `expdesign` need to be provided. If none of the above is provided, it is possible to upload the data during the execution of the Shiny App
- `res_obj` -  a `DESeqResults` object. If not provided, it can be computed during the execution of the application
- `annotation_obj` - a `data.frame` object, with row.names as gene identifiers (e.g. ENSEMBL ids) and a column, `gene_name`, containing e.g. HGNC-based gene
symbols. If not provided, it can be constructed during the execution via the `org.eg.XX.db` packages
- `countmatrix` - a count matrix, with genes as rows and samples as columns. If not provided, it is possible to upload the data during the execution of the Shiny App
- `expdesign` -a `data.frame` containing the info on the experimental covariates of each sample. If not provided, it is possible to upload the data during the execution of the Shiny App

Different modalities are supported to launch the application:

- `WEIN(dds_obj = dds, res_obj = res, annotation_obj = anno)`, where the objects are precomputed in the current session and provided as parameters
- `WEIN(dds_obj = dds)`, as in the command above, but where the result object is assembled at runtime 
- `WEIN(countmatrix = countmatrix, expdesign = expdesign)`, where instead of passing the defined `DESeqDataSet` object, its components are given, namely the count matrix (e.g. generated after a run of featureCounts or HTSeq-count) and a data frame with the experimental covariates. The design formula can be constructed interactively at runtime
- `WEIN()`, where the count matrix and experimental design can simply be uploaded at runtime, where all the derived objects can be extracted and computed live. These files have to be formatted as tabular text files, and a function in the package tries to guess the separator, based on heuristics of occurrencies per line of commonly used characters

## Getting to know the user interace and the functionality

The user interface is dynamically displayed according to the provided and computed objects, with tabs that are actively usable only once the required input is effectively available.

Moreover, for some relevant UI widgets, the user can receive additional information by hovering over with the mouse, with the functionality powered by the *[shinyBS](http://cran.fhcrc.org/web/packages/shinyBS/index.html)* package.

For the user which is either new with the app UI/functionality, or not extensively familiar with the topic of differential expression, it is possible to obtain a small *guided tour* of the App by clicking on the respective help buttons, marked in the app like this.

<button id="btn" type="button" class="btn btn-default action-button" style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4">
  <i class="fa fa-info"></i>
  Click me for a quick tour
</button>

These trigger the start of a step-by-step guide and feature introduction, powered by the *[rintrojs](http://cran.fhcrc.org/web/packages/rintrojs/index.html)* package.

### The controls sidebar

Some of the input controls which affect different tabs are located in the sidebar, while others are as well in the individual tabs of the app. By changing one or more of the input parameters, the user can get a fine control on what is computed and displayed.

#### App settings

- **Group/color by** - Select the group of samples to stratify the analysis for plotting. Can also assume multiple values.
- **Select the gene(s) of interest - ids** - Select a subset of genes for deeper analysis. If an annotation object is provided, the user can handily select the genes e.g. based on their HGNC symbol
- **False Discovery Rate** - Set as default to 0.05, it is the FDR value for the Benjamini-Hochberg procedure for adjusting p-values in the multiple testing comparison scenario

#### Plot export settings   

**Width** and **Height** for the figures to export are input here in cm.

#### Quick viewer

This displays a list of the underlying objects with which basically all of the analysis can be performed. A green tick icon appears close to each when the respective component is either provided or calculated. For obtaining the best analysis experience in `WEIN`, it is recommended to provide all of them.

#### First steps help

Clicking on this button activated the `intro.js` based tour for getting to know the components and the structure of the app. Dedicated step-by-step procedures are also available in each individual tab.

### The task menu

The task menu, accessible by clicking on the cog icon in the upper right part of the application, provides two functionalities:

- `Exit WEIN & save` will close the application and store the content of the `input` and `values` reactive objects in a list of two elements in the `ideal_env` environment, respectively called `ideal_inputs_YYYYMMDD_HHMMSS` and `ideal_values_YYYYMMDD_HHMMSS`
- `Save State as .RData` will similarly store `LiveInputs` and `r_data` in a binary file named `idealState_YYYYMMDD_HHMMSS.Rdata`, without closing the application 

## The main app panels

The *WEIN* app is a one-paged dashboard, structured in different panels, where each of them is focused on a different aspect of the data exploration. 

On top of the panels, three `valueBox` objects serve as guiding elements for having an overview of the data at hand: how many genes and samples are in the data, how many entries are in the annotation object, and how many genes were found to be differentially expressed in the results. Whenever each of the underlying objects is available, the background color turns from red to green.

For the main analysis, the available panels are described in the following subsections.

### Welcome!

The landing page for the app is also where you might likely be reading this text (otherwise in the package vignette).

### Data Setup

The Data Setup panel is where you can upload or inspect the required inputs for running the app. This builds on the primary idea used by *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)* and extends it with the following aspects:

- the panel structure appears dynamically in three consecutive mandatory steps, marked with color from red to yellow to green, with optional steps in light blue.
- the optional step of retrieving the annotation on the fly relieves the user from the task of composing the `data.frame` in advance, and is based on the widely adopted `org.XX.eg.db` Bioconductor packages.
- when the objects are already passed as parameters, or computed, a brief overview/summary for them is displayed
- to tighten the concert operations between similar tools with different scope (as *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)* and the original *[ideal](http://bioconductor.org/packages/ideal)* are), the information flow can move from the data exploration to decisions taken at the moment of testing

A diagnostic mean-dispersion plot is also provided in a collapsible element at the bottom of the panel, shown when the `DESeqDataSet` is generated and the `DESeq` command from the `DESeq2` package has been applied.

### Counts Overview

As in pcaExplorer, interactive tables for the raw, normalized or (r)log-transformed counts are shown in this tab. The user can also generate a sample-to-sample correlation scatter plot with the selected data.

Additionally, `WEIN` has an option to include a filter step at the gene level by removing genes with low absolute or averages low values. After this, it might be possible to have to re-run the analysis in step 3 from the Data Setup panel. 

### Extract Results

This tab is an interface for generating the summary tables after testing for DE. It is usually based on the Wald test, as implemented in DESeq2, but when the factor of interest is assuming more than two levels, the user can also perform an ANOVA-like test across the groups with the likelihood ratio test. Options for enabling/disabling automated independent filtering, adding the additional column of unshrunken log2 fold change values (instead of the moderated estimates used by default), as well as using the Independent Hypothesis Weighting (IHW) framework, are provided.

The False Discovery Rate (FDR) can be set from the sidebar panel, and a couple of diagnostic plots, such as the histogram of raw p-values and the distribution of log2fc, are shown below the interactive enhanced version of the table - with clickable elements to link to ENSEMBL database and NCBI website.

### Summary Plots

In this tab an interactive MA plot for the contrast selected in the Extract Results tab is displayed. Clicking on a single gene in the zoomed plot (enabled by brushing in the main plot), it is possible to obtain a boxplot for its expression values, flanked by an overview of information accessed live from the Entrez database. Alternatively, a volcano plot of -log10(p-value) versus log fold change can provide a slightly different perspective. The subset of selected genes are also here presented in static and interactive heatmaps, with the underlying data accessible from the collapsible box element.

### Gene Finder

The functionality in the Gene Finder builds upon the one provided by `pcaExplorer`, and allows to query up to four genes in the same view, which can here be selected from a dropdown input list which supports autocompletion. 

A combined summary table (with both normalized counts and results statistics) is located below an MA plot where the selected genes are marked and annotated on the plot. To avoid repeating this manually, the user can also quickly upload a list of genes as text file (one gene identifier per line), such as members of gene families (e.g. all cytokines, all immunoglobines, ...) or defined by common function (e.g. all housekeeping genes, or others based on any annotation).

### Functional Analysis

The Functional Analysis tab takes the user from the simple lists of DE genes to insight on the affected biological pathways, with three approaches based on the Gene Ontology (GO) databases. This panel of WEIN has a slim interface to 

- `limma::goana` for the quick yet standard implementation
- `topGO`, particularly valuable for pruning terms which are topologically less meaningful than their specific nodes
- `goseq`, which accounts for the specific length bias intrinsic in RNA-Seq assays (longer genes have higher chances of being called DE).

*WEIN* allows the user to work simultaneously with more gene lists, two of which can be uploaded in a custom way (e.g. list of gene families, or extracted from other existing publications). 

The interaction among these lists can be visually represented in Venn diagrams, as well as with the appealing alternative from the UpSetR package, where all combination of sets are explicitly shown. 

Each of the methods for GO enrichment delivers its own interactive `DT`-based table, which can then be explored interactively with the display of a heatmap for all the (DE) genes annotated to a particular term, picking the normalized transformed values for comparing robustly the expression values. This is simply triggered by clicking any of the rows for the results tables. Another useful feature is provided by the clickable link to the AmiGO database on each of the GO term identifiers.

### Report Editor

The Report Editor tab works in the same way of `pcaExplorer`, with the scope of providing an interface to full computational reproducibility of the analyses.

General `Markdown options` and `Editor options` are available, and the text editor, based on the `shinyAce` package, contains a comprehensive template report, that can be edited to the best convenience of the user.

The code contained in the template report fetches the latest state of the reactive values in the ongoing session, and its output is a comprehensive HTML file that can be expanded, edited, previewed in the tab itself, downloaded, and shared with a few mouse clicks.

### About

The About tab contains the output of `sessionInfo`, plus general information on *WEIN*, including the link to the Github development version. If requested, the modular structure of the app can be easily expanded, and many new operations on the same set of input data and derived results can be embedded in the same framework. 

## Functions exported by the package for standalone usage

The functions exported by the *WEIN* package can be also used in a standalone scenario, provided the required objects are in the working environment. They are listed here for an overview, but please refer to the documentation for additional details. Where possible, for each function a code snippet will be provided for its typical usage.

### `deseqresult2DEgenes` and `deseqresult2tbl`

`deseqresult2DEgenes` and `deseqresult2tbl` generate a tidy table with the results of DESeq2, sorted by the values in the `padj` column.

### `ggplotCounts`

`ggplotCounts` extends the functionality of the `plotCounts` function of *[DESeq2](http://bioconductor.org/packages/DESeq2)*, and plots the normalized counts of a single gene as a boxplot, with jittered points superimposed.

If an `annotation_obj` is provided, their gene name can also be included in the title.

When used in the context of the app, it is possible to seamless search for genes also by their gene name, making exploration even more immediate.

### `goseqTable`

`goseqTable` is a wrapper to extract the functional GO terms enriched in  in a list of (DE) genes, based on the algorithm and the implementation in the *[goseq](http://bioconductor.org/packages/goseq)* package.

Its counterpart, based on the *[topGO](http://bioconductor.org/packages/topGO)* package, can be found in the *[pcaExplorer](http://bioconductor.org/packages/pcaExplorer)* package.

### `plot_ma`

The MA plot provided by *WEIN* displays the gene-wise log2-fold changes (logFCs) versus the average expression value. As a main input parameter, a `DESeqResults` object is required. Control on the appearance of the plot can be applied by selecting the False Discovery Rate (`FDR`), the point transparency (`point_alpha`), adding horizontal lines at particular logFC values (`hlines`). The user can also decide to add rug plots in the margins (setting `add_rug` to `TRUE`).

To facilitate the inspection of a particular gene or gene set, `intgenes` can assume the value of a vector of genes (either the IDs or the gene symbols if `symbol` column is provided in `res_obj`. Labels can be added via `labels_intgenes`, while classical labels/title can be also edited as preferred (see `plot_ma` for all details). 

### `plot_volcano`

The volcano plot gives a different flavor for the gene overview, displaying log2-fold changes and log p-values

In a way similar to `plot_ma`, genes can be annotated with `intgenes`, and vertical lines can be added via `vlines`. `ylim_up` controls the y axis upper limit to visualize better the bulk of genes - keep in mind that very small p-values due to robust differences/large effect sizes can be "cut out".

### `sepguesser`   

`sepguesser` makes an educated guess on the separator character for the input text file (`file`). The separator list can be provided as a vector in `sep_list` (defaults to comma, tab, semicolon, and whitespace - which ideally could cover most of the cases). The heuristics is based on the number of occurrencies of each separator in each line.

## Creating and sharing output objects

While running the app, the user can

- generate and save graphics
- create and export tables
- generate, preview, download/export an HTML report
- save the values of the `reactiveValues` in an environment, or in binary format

This functionality to retrieve and share the output is provided by action buttons that are placed close to each element of interest. 

## Further development

Additional functionality for *WEIN* will be added in the future, as it continues to evolve from the original *ideal* package with enhancements for immunology-focused transcriptomics analysis.

Improvements, suggestions, bugs, issues and feedback of any type can be sent to bernd.jagla@pasteur.fr.

## Session Info {.unnumbered}

```r
sessionInfo()
# devtools::session_info()
