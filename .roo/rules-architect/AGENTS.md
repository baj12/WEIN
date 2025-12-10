# Project Architecture Rules (Non-Obvious Only)

- App follows modular Shiny architecture with clear separation of UI and server components
- All modules use namespacing via NS() function for ID isolation to prevent conflicts
- Reactive state is centralized in a single 'values' reactiveValues object passed to all modules
- Modules follow strict naming convention: *_ui for UI components and *_server for server logic
- Server modules must accept input, output, session, and values parameters in that order
- UI components are pure functions that return UI definitions without side effects
- Server logic encapsulates all reactive behavior and data manipulation
- Bookmarking is implemented at the main server level with onBookmark/onRestore callbacks
- State persistence uses global environment 'ideal_env' for cross-session storage
- Export functionality stores plots in exportPlots reactiveValues container
- Species annotation mapping is centralized in annoSpecies_df data frame
- Functional analysis modules integrate with multiple Bioconductor packages (GOseq, topGO, etc.)
- Report generation uses R Markdown templates with dynamic content injection
- Data flow follows Shiny reactive principles with explicit dependency chains
- Memory management uses temporary directories (tempfile()) for user-specific files