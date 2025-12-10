about_ui <- function(){
  tabPanel(
    "About", icon = icon("institution", verify_fa = FALSE),
    
    # headerPanel("Information on ideal/session"),
    
    fluidRow(
      column(
        width = 8,
        includeMarkdown(system.file("extdata", "about.md",package = "WEIN")),
        
        verbatimTextOutput("sessioninfo")
      )
    )
  )
}