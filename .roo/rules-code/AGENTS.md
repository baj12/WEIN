# Project Coding Rules (Non-Obvious Only)

- All Shiny modules follow strict naming convention: *_ui.R for UI components and *_server.R for server logic
- UI functions must be exported and follow pattern function_name_ui()
- Server functions must be exported and follow pattern function_name_server()
- All modules use namespacing via NS() function for ID isolation
- Reactive values are passed as 'values' parameter to all server modules
- Heavy use of require() inside functions rather than top-level library() calls
- Global environment 'ideal_env' is used for state persistence across sessions
- Bookmarking is implemented via onBookmark/onRestore callbacks in main server function
- Exported functions must have Roxygen documentation with @export tag
- Helper functions should be placed in R/helpers.R file
- UI components are organized in separate files following *_ui.R pattern
- Server logic is organized in separate files following *_server.R pattern
- All server modules must accept input, output, session, and values parameters
- Use reactiveValues() for state management instead of regular variables