dash_Header_ui <- function(){
  shinydashboard::dashboardHeader(
    title = tags$span(
      img(src = "WEIN_logo_v2.png", height = "50px"),
      paste0("WEIN - Web-based Engine for Interactive Next-generation sequencing analysis ",
             packageVersion("WEIN"))),
    titleWidth = 600,
    
    # Logo implementation reference:
    # http://stackoverflow.com/questions/31440564/adding-a-company-logo-to-shinydashboard-header
    # Current implementation uses text title with embedded image
    
    # task menu for saving state to environment or binary data
    shinydashboard::dropdownMenu(
      type = "tasks",icon = icon("cog", verify_fa = FALSE),
      badgeStatus = NULL,
      headerText = "WEIN Tasks menu",
      notificationItem(
        text = actionButton("task_exit_and_save","Exit WEIN & save",
                            class = "btn_no_border",
                            onclick = "setTimeout(function(){window.close();}, 100); "),
        icon = icon("right-from-bracket"),status = "primary"),
      menuItem(
        text = downloadButton("task_state_save","Save State as .RData"))
    )
  )
}