sideBar_ui <- function(){
  dashboardSidebar(
    width = 340,
    {
      nsSetup <- NS("ui_setup")
    sidebarMenu(style = "position: fixed; overflow: visible;",
                id = "sidebarMenu",
                menuItem("App settings",
                         icon = icon("cogs", verify_fa = FALSE),
                         startExpanded = TRUE,
                         uiOutput("color_by"),
                         shinyBS::bsTooltip(
                           "color_by", 
                           paste0("Select the group(s) of samples to stratify the analysis, and ideally match the contrast of interest. Can also assume multiple values, in this case the interaction of the factors is used."),
                           "right", options = list(container = "body")),
                         uiOutput("available_genes"),
                         shinyBS::bsTooltip(
                           "available_genes", 
                           paste0("Select one or more features (genes) from the list to inspect. Autocompletion is provided, so you can easily find your genes of interest by started typing their names. Defaults to the row names if no annotation object is provided."),
                           "right", options = list(container = "body")),
                         numericInput("FDR","False Discovery Rate",value = 0.05, min = 0, max = 1, step = 0.01),
                         shinyBS::bsTooltip(
                           "FDR", 
                           paste0("Select the alpha level at which you would like to control the FDR (False Discovery Rate) for the set of multiple tests in your dataset. The sensible choice of 0.05 is provided as default, 0.1 is more liberal, while 0.01 is more stringent - keep in mind this does not tell anything on the effect size for the expression change."),
                           "right", options = list(container = "body"))
                         
                ),
                menuItem("Plot export settings", 
                         icon = icon("paint-brush", verify_fa = FALSE),
                         startExpanded = TRUE,
                         numericInput("export_width",label = "Width of exported figures (cm)",value = 16,min = 2),
                         shinyBS::bsTooltip(
                           "export_width", paste0("Width of the figures to export, expressed in cm"),
                           "right", options = list(container = "body")),
                         numericInput("export_height",label = "Height of exported figures (cm)",value = 10,min = 2),
                         shinyBS::bsTooltip(
                           "export_height", paste0("Height of the figures to export, expressed in cm"),
                           "right", options = list(container = "body"))
                ),
                menuItem("Quick viewer", 
                         icon = icon("bolt", verify_fa = FALSE), 
                         startExpanded = TRUE,
                         id = "qvmenu",
                         fluidRow(column(12,
                                         fluidRow(column(6,p("Count matrix")), column(6,  uiOutput(nsSetup("ok_cm")))),
                                         fluidRow(column(6,p("Experimental design")), column(6,uiOutput(nsSetup("ok_ed")))),
                                         fluidRow(column(6,p("Annotation")), column(6,uiOutput(nsSetup("ok_anno")))),
                                         # Maybe we need this when having preoaded count matrix?
                                         # fluidRow(column(6,p("DESeqDataset")), column(6,uiOutput("ok_dds"))),
                                         fluidRow(column(6,p("DESeqRun")), column(6,uiOutput(nsSetup("ok_ddsRun")))),
                                         fluidRow(column(6,p("Results")), column(6,uiOutput(nsSetup("ok_resu"))))
                         ))),
                menuItem("First steps help", 
                         icon = icon("question-circle", verify_fa = FALSE),
                         startExpanded = TRUE,
                         actionButton("btn", "Click me for a quick tour", icon("info"),
                                      style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4")
                ),
                bookmarkButton()
    )
    }
  )
}