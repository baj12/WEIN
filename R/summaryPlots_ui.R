summaryPlots_ui <-  function(){
  nsSP <- NS("summary_plots")
  nsTR <- NS("tour_manager")
  tabPanel(
    "Summary Plots", icon = icon("camera", verify_fa = FALSE),
    conditionalPanel(
      condition="!output.checkresu",
      headerPanel("Interactive graphical exploration of the results"),
      fluidRow(
        column(
          width = 8,
          shinyBS::bsCollapse(
            id = nsTR("help_summaryplots"),open = NULL, 
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_plots.md",package = "WEIN")))
          )
        )
      ),
      
      actionButton(nsSP("tour_plots"), "Click me for a quick tour of the section", icon("info"),
                   style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
      
      br(),
      fluidRow(
        column(12,
               sliderInput(nsSP("max_points"), 
                           "Maximum points to display (for performance):",
                           min = 1000, max = 20000, value = 5000, step = 1000,
                           width = "100%")
        )
      ),
      br(),
      fluidRow(
        column(6,
               h4("MA plot - Interactive!"),
               plotOutput(nsSP('plotma'), brush = nsSP('ma_brush')),
               # plotOutput(nsSP('plotma')),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsSP("download_plot_ma"), "Download Plot"),
                   textInput(nsSP("filename_plot_ma"),label = "Save as...",value = "plot_ma.pdf"))),
        column(6,
               h4("Zoomed section"),
               plotOutput(nsSP("mazoom"),click= nsSP('mazoom_click')),
               numericInput(nsSP('size_genelabels'), label = 'Labels size: ', value = 4,min = 1,max = 8),
               textAreaInput(nsSP('gene_list_input'), label = 'Add genes (space-separated):', value = '', placeholder = 'Paste gene names here...', rows = 3),
               actionButton(nsSP('add_genes_button'), 'Add Genes to Selection'),
               actionButton(nsSP('clear_genes_button'), 'Clear Gene Selection'),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsSP("download_plot_mazoom"), "Download Plot"),
                   textInput(nsSP("filename_plot_mazoom"),label = "Save as...",value = "plot_mazoom.pdf")))
      ),
      fluidRow(
        column(6,
               h4("Selected gene"),
               checkboxInput(nsSP("ylimZero_genes"),"Set y axis limit to 0",value=TRUE),
               plotOutput(nsSP("genefinder_plot")),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsSP("download_plot_genefinder"), "Download Plot"),
                   textInput(nsSP("filename_plot_genefinder"),label = "Save as...",value = "plot_genefinder.pdf"))
        ),
        column(6,
               h4("Gene infobox"),
               htmlOutput(nsSP("rentrez_infobox")))
      ),
      
      fluidRow(
        column(6,
               h4("volcano plot"),
               plotOutput(nsSP("volcanoplot")),
               checkboxInput(nsSP("show_gene_names"), "Show gene names", value = TRUE),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsSP("download_plot_volcanoplot"), "Download Plot"),
                   textInput(nsSP("filename_plot_volcanoplot"),label = "Save as...",value = "plot_volcanoplot.pdf"))
        )),
      
      fluidRow(radioButtons(nsSP("heatmap_colv"),"Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
      fluidRow(
        column(4,
               checkboxInput(nsSP("rowscale"),label = "Scale by rows",value = TRUE)),
        column(4,
               checkboxInput(nsSP("pseudocounts"),"use log2(1+counts)",value = TRUE))
      ),
      fluidRow(
         column(12,
               shinyjqui::jqui_resizable(
                 plotlyOutput(nsSP("hpi_brush")))
        )
      ),
      
      box(
        title = "Brushed table", status = "primary", solidHeader = TRUE,
        id = "box_brushedtbl",
        collapsible = TRUE, collapsed = TRUE, width = 12,
        fluidRow(DT::dataTableOutput(nsSP("ma_brush_out")),
                 downloadButton(nsSP("downloadTblMabrush"),"Download", class = "btn btn-success")),
        fluidRow(
          verbatimTextOutput(nsSP("selectedGenesMAplot"))
        ))
      
    ),
    conditionalPanel(
      condition="output.checkresu",
      h2("You did not create the result object yet. Please go the dedicated tab and generate it")
    )
  )
}