geneFinder_ui <- function(){
  nsGF <- NS("gene_finder")
  tabPanel(
    "Gene Finder", icon = icon("crosshairs"),
    conditionalPanel(
      condition="!output.checkdds",
      headerPanel("Find your gene(s) of interest"),
      fluidRow(
        column(
          width = 8,
          shinyBS::bsCollapse(
            id = "help_genefinder",open = NULL,
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_genefinder.md",package = "WEIN")))
          )
        )
      ),
      actionButton("tour_genefinder", "Click me for a quick tour of the section", icon("info"),
                   style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
      br(),
      fluidRow(
        column(6,checkboxInput(nsGF("ylimZero_genefinder"),"Set y axis limit to 0",value=TRUE))),
      fluidRow(
        column(6,
               shinyjqui::jqui_resizable(plotOutput(nsGF("bp1"))),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsGF("download_plotbp1"), "Download Plot"),
                   textInput(nsGF("filename_plotbp1"),label = "Save as...",value = "plotbp1.pdf"))
        ),
        column(6,
               shinyjqui::jqui_resizable(plotOutput(nsGF("bp2"))),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsGF("download_plotbp2"), "Download Plot"),
                   textInput(nsGF("filename_plotbp2"),label = "Save as...",value = "plotbp2.pdf")))
      ),
      fluidRow(
        column(6,
               shinyjqui::jqui_resizable(plotOutput(nsGF("bp3"))),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsGF("download_plotbp3"), "Download Plot"),
                   textInput(nsGF("filename_plotbp3"),label = "Save as...",value = "plotbp3.pdf"))
        ),
        column(6,
               shinyjqui::jqui_resizable(plotOutput(nsGF("bp4"))),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsGF("download_plotbp4"), "Download Plot"),
                   textInput(nsGF("filename_plotbp4"),label = "Save as...",value = "plotbp4.pdf")))
      ),
      
      fluidRow(
        column(6,
               shinyjqui::jqui_resizable(plotOutput(nsGF("plotCoefficients")))
        )
      ),
      fluidRow(
        column(width = 10,offset = 1,
               plotOutput(nsGF("ma_highlight")),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsGF("download_plot_mahighlight"), "Download Plot"),
                   textInput(nsGF("filename_plot_mahighlight"),label = "Save as...",value = "plot_mahighlight.pdf")),
               DT::dataTableOutput(nsGF("table_combi")),
               downloadButton(nsGF("downloadTblCombi"),"Download", class = "btn btn-success"),
               
               fileInput(inputId = nsGF("gl_ma"),
                         label = "Upload a gene list file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE),
               plotOutput(nsGF("ma_hl_list")),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsGF("download_plot_mahllist"), "Download Plot"),
                   textInput(nsGF("filename_plot_mahllist"),label = "Save as...",value = "plot_mahllist.pdf")),
               DT::dataTableOutput(nsGF("table_combi_list")),
               downloadButton(nsGF("downloadTblCombiList"),"Download", class = "btn btn-success")
        )
      )
    ),
    conditionalPanel(
      condition="output.checkdds",
      h2("You did not create the dds object yet. Please go the main tab and generate it")
    )
  )
}