#' UI Components for WEIN
#'
#' Contains all UI-related functions for the WEIN Shiny application.
#'
#' @name WEIN-ui
#' @docType package
#' @keywords internal
"_PACKAGE"

#' Create the UI for WEIN
#'
#' Defines the complete UI structure for the WEIN Shiny application.
#'
#' @param req The Shiny request object
#' @param values Reactive values containing application state
#' @param annoSpecies_df Data frame with species annotation information
#' 
#' @return A Shiny UI definition
#' @export
WEIN_ui <- function() {
  require(shiny)
  require(shinydashboard)
  nsUI <- NS("ui_setup")
  shinydashboard::dashboardPage(
    title = "WEIN - Web-based Engine for Interactive Next-generation sequencing analysis",
    dash_Header_ui(),
    
    sideBar_ui(), 
    
    dashboardBody(
      introjsUI(),
      
      ## Define output size and style of error messages, and also the style of the icons e.g. check
      ## plus, define the myscrollbox div to prevent y overflow when page fills up
      tags_head_ui(),
      
      # value boxes to always have an overview on the available data
      fluidRow(
        valueBoxOutput(nsUI("box_ddsobj")),
        valueBoxOutput(nsUI("box_annobj")),
        valueBoxOutput(nsUI("box_resobj"))
      ),
      
      ## main structure of the body for the dashboard
      div(
        id = "myScrollBox", # trick to have the y direction scrollable
        tabBox(
          width=12,
          welcome_panel_ui(), 
          data_setup_ui(),
          panel_count_overview_ui(), 
          extract_results_ui(),
          summaryPlots_ui(), 
          geneFinder_ui(), 
          panalAnalysis_ui(), 
          signature_ui(),
          report_editor_ui(),
          
          about_ui() 
          # end of box
        )
      ) # end of myScrollBox
      ,footer()
    ), # end of dashboardBody
    skin="black"
  ) # end of dashboardPage
  
}
