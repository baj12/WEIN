# Project Debug Rules (Non-Obvious Only)

- testData.RData required for some test scripts but not committed to repo
- Global variable ideal_env used for state persistence across sessions
- Bookmarking state is stored in values$restoreBookmark flag
- Shiny app state can be saved to global environment with "task_exit_and_save" button
- Reactive values are flushed to ideal_env when exiting the app
- Session information can be viewed with output$sessioninfo
- State saving creates timestamped objects in ideal_env (ideal_inputs_TIMESTAMP, ideal_values_TIMESTAMP)
- Bookmarking saves all reactive values in state$values$ValList
- Export plots are stored in exportPlots reactiveValues container
- Debugging requires understanding of Shiny's reactive programming model
- Many Bioconductor dependencies must be installed separately for full functionality