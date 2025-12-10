# Project Documentation Rules (Non-Obvious Only)

- Package is based on the ideal project by Federico Marini, extended by Bernd Jagla at Institut Pasteur
- Main function is idealImmunoTP() which launches the Shiny app
- App has modular structure with separate UI and server components
- UI modules are in R/*_ui.R files and server modules in R/*_server.R files
- Core data structures are DESeqDataSet, DESeqResults, and annotation data frames
- Gene signatures can be provided as lists of vectors (e.g., from read_gmt function)
- App supports bookmarking to save and restore state
- State can be saved to binary .RData files or to global R environment
- App uses reactiveValues for state management across modules
- Modules communicate through the shared 'values' reactiveValues object
- UI components use namespacing via NS() function for ID isolation
- App supports multiple species through annoSpecies_df lookup table
- Functional analysis uses GOseq, topGO, and other Bioconductor packages
- Report generation uses R Markdown templates stored in inst/extdata/irt.Rmd