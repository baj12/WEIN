
tourModuleServer <- function(id, input, session) {
  moduleServer(id, function(input, output, session) {
    # Tour configurations: input_id -> filename
    tour_configs <- list(
      btn = "intro_firsttour.txt",
      tour_datasetup = "intro_datasetup.txt",
      tour_countsoverview = "intro_countsoverview.txt",
      tour_results = "intro_results.txt",
      tour_plots = "intro_plots.txt",
      tour_genefinder = "intro_genefinder.txt",
      tour_funcanalysis = "intro_funcanalysis.txt",
      tour_signatureexplorer = "intro_signatureexplorer.txt",
      tour_report = "intro_report.txt"
    )
    
    # Create observers for file-based tours
    lapply(names(tour_configs), function(tour_id) {
      observeEvent(input[[tour_id]], {
        tour_data <- read.delim(
          system.file("extdata", tour_configs[[tour_id]], package = "idealImmunoTP"),
          sep = ";", stringsAsFactors = FALSE
        )
        introjs(session, options = list(steps = tour_data))
      }, ignoreNULL = TRUE)
    })
    
    # Special observer for introexample
    observeEvent(input$introexample, {
      intro_example <- data.frame(
        element = c("#introexample", "#introexample"),
        intro = c(
          "Tour elements can be anchored to elements of the UI that are intended to be highlighted...",
          "Well done. This is how a tour can look like. Click outside of this window..."
        )
      )
      introjs(session, options = list(steps = intro_example))
    }, ignoreNULL = TRUE)
    
    output$ui_instructions <- renderUI({
      box(width = 12, 
          title = "Instructions", status = "info", solidHeader = TRUE, 
          collapsible = TRUE, collapsed = TRUE,
          includeMarkdown(system.file("extdata", "instructions.md",package = "idealImmunoTP"))
      )
    })
    
  })
}
