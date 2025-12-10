# AGENTS.md

This file provides guidance to agents when working with code in this repository.

## Stack
- R package using Bioconductor, Shiny, and Shiny Dashboard
- Dependencies include DESeq2, ggplot2, DT, plotly, and many Bioconductor packages
- Tests use testthat framework

## Build/Test Commands
- Run tests: `Rscript -e "devtools::test()"`
- Check package: `Rscript -e "devtools::check()"`
- Install package: `R CMD INSTALL .` or `devtools::install()`
- Run app: `idealImmunoTP()` after loading library
- Run tests with data: Source runTest.R after ensuring testData.RData exists

## Project Structure
- Modular Shiny app with separate UI and server components
- UI modules in R/*_ui.R files
- Server modules in R/*_server.R files
- Core app functions in R/idealImmunoTP.R
- Helper functions in R/helpers.R
- Tests in tests/testthat/

## Critical Patterns
- Uses Shiny modules with namespace pattern (e.g., ui_setup_server/ui_setup_ui)
- Heavy use of reactiveValues for state management
- Bookmarking support via onBookmark/onRestore callbacks
- Extensive use of shinydashboard for UI layout
- All server modules follow pattern: function(input, output, session, values, ...)
- UI components use NS() for namespacing
- Global environment ideal_env used for state persistence
- Extensive use of Bioconductor packages for genomic analysis

## Code Style
- Roxygen2 documentation for all exported functions
- Function names use snake_case
- Module functions follow naming pattern *_server/*_ui
- Heavy use of require() inside functions rather than top-level library()
- Extensive import declarations in NAMESPACE file
- Comments use # for single line, longer comments for sections

## Testing Specifics
- Tests use testthat framework
- Test files named test_*.R in tests/testthat/
- Tests typically create small example datasets with DESeq2::makeExampleDESeqDataSet()
- Integration tests in test_integration.R
- Component tests in test_components.R
- Utility function tests in test_utils.R
- Shiny-specific tests in test_shiny.R

## Gotchas
- Many Bioconductor dependencies must be installed separately
- testData.RData required for some test scripts but not committed to repo
- Global variable ideal_env used for state persistence across sessions
- Extensive use of require() inside functions means dependencies checked at runtime
- Package assumes many org.*.eg.db packages may be installed