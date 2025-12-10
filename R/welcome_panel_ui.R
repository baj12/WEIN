welcome_panel_ui = function(){
  nsTour <- NS("tour_manager")
  tabPanel(
    title = "Welcome!",  icon = icon("house"), value="tab-welcome",
    
    fluidRow(
      column(
        width = 8,
        includeMarkdown(system.file("extdata", "welcome.md",package = "idealImmunoTP")),
        br(),br(),
        p("If you see a grey box like this one open below..."),
        
        shinyBS::bsCollapse(
          id = "help_welcome",open = "Help", 
          shinyBS::bsCollapsePanel(
            "Help", 
            includeMarkdown(system.file("extdata", "help_welcome.md",package = "idealImmunoTP"))
          )
        ),
        
        actionButton(nsTour("introexample"), "If you see a button like this...", icon("info"),
                     style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
        p("... you can click on that to start a tour based on introJS"),
        br(),br(),
        
        uiOutput(nsTour("ui_instructions"))
      )
    )
  )
}