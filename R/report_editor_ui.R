report_editor_ui <- function() {NULL}
# # ui panel report editor -----------------------------------------------------------
# tabPanel(
#   "Report Editor",
#   icon = icon("pencil"),
#   headerPanel("Create, view and export a report of your analysis"),
#   fluidRow(
#     column(
#       width = 8,
#       shinyBS::bsCollapse(
#         id = "help_reporteditor",open = NULL, 
#         shinyBS::bsCollapsePanel(
#           "Help",
#           includeMarkdown(system.file("extdata", "help_report.md",package = "idealImmunoTP")))
#       )
#     )
#   ),
#   
#   actionButton("tour_report", "Click me for a quick tour of the section", icon("info"),
#                style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),
#   
#   fluidRow(
#     column(
#       width = 6,
#       box(
#         title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
#         id = "md_opts",
#         radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = TRUE),
#         textInput("report_title", "Title: "),
#         textInput("report_author", "Author: "),
#         radioButtons("report_toc", "Table of Contents", choices = list("Yes" = "true", "No" = "false")),
#         radioButtons("report_ns", "Number sections", choices = list("Yes" = "true", "No" = "false")),
#         selectInput("report_theme", "Theme", choices = list("Default" = "default", "Cerulean" = "cerulean",
#                                                             "Journal" = "journal", "Flatly" = "flatly",
#                                                             "Readable" = "readable", "Spacelab" = "spacelab",
#                                                             "United" = "united", "Cosmo" = "cosmo")),
#         radioButtons("report_echo", "Echo the commands in the output", choices = list("Yes" = "TRUE", "No" = "FALSE")))),
#     column(
#       width = 6,
#       box(
#         title = "editor options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
#         id = "editor_opts",
#         checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
#         conditionalPanel(
#           "input.enableAutocomplete",
#           wellPanel(
#             checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
#             checkboxInput("enableRCompletion", "R code completion", TRUE)
#           )
#         ),
#         
#         selectInput("mode", "Mode: ", choices=shinyAce::getAceModes(), selected="markdown"),
#         selectInput("theme", "Theme: ", choices=shinyAce::getAceThemes(), selected="solarized_light"))
#     )
#   ),
#   fluidRow(
#     column(3,
#            actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
#     ),
#     column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success")),
#     column(3, uiOutput("ui_iSEEexport"))
#   ),
#   
#   tabBox(
#     width = NULL,
#     id="report_tabbox",
#     tabPanel("Report preview",
#              icon = icon("file-lines"),
#              htmlOutput("knitDoc")
#     ),
#     
#     tabPanel("Edit report",
#              icon = icon("pen-to-square"),
#              aceEditor("acereport_rmd", mode="markdown",theme = "solarized_light",autoComplete = "live",
#                        value="_Initialization of the_ `idealImmunoTP` _report generation..._",
#                        placeholder = "You can enter some code and text in R Markdown format",
#                        height="800px"))
#   )
# ), # end of Report Editor panel