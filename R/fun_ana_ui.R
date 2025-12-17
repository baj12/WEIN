
panalAnalysis_ui <- function(){
  nsFA <- NS("fun_ana")
  nsTR <- NS("tour_manager")
  
  tabPanel(
    "Functional Analysis", icon = icon("rectangle-list"),
    conditionalPanel(
      condition="!output.checkresu",
      headerPanel("Find functions enriched in gene sets"),
      fluidRow(
        column(
          width = 8,
          shinyBS::bsCollapse(
            id = nsTR("help_functionalanalysis"), open = NULL,
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_funcanalysis.md", package = "WEIN")))
          )
        )
      ),
      actionButton(nsTR("tour_funcanalysis"), "Click me for a quick tour of the section", icon("info"),
                   style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),
      
      selectInput(nsFA("go_cats"), label = "Select the GO category(ies) of interest",
                  choices = list("GO Biological Process" = "BP", "GO Molecular Function" = "MF", "GO Cellular Component" = "CC"),
                  selected = "BP", multiple = TRUE
      ),
      
      div(
        id = "myAnchorBox",
        tabBox(
          width = NULL,
          id = "gse_tabbox",
          tabPanel("UPregu", icon = icon("arrow-circle-up", verify_fa = FALSE),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrUP"), "Perform gene set enrichment analysis on the upregulated genes", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrUP_goseq"), "Perform gene set enrichment analysis on the upregulated genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrUP_topgo"), "Perform gene set enrichment analysis on the upregulated genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_up")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_up_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_up_topgo")),
                   downloadButton(nsFA("downloadGOTbl_up"), "Download", class = "btn btn-success"),
                   fluidRow(h3("Heatmap of selected TopGo entry")),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_up_topgo")))))
          ),
          tabPanel("DOWNregu", icon = icon("arrow-circle-down", verify_fa = FALSE),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrDOWN"), "Perform gene set enrichment analysis on the downregulated genes", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrDOWN_goseq"), "Perform gene set enrichment analysis on the downregulated genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrDOWN_topgo"), "Perform gene set enrichment analysis on the downregulated genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_down")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_down_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_down_topgo")),
                   downloadButton(nsFA("downloadGOTbl_down"), "Download", class = "btn btn-success"),
                   fluidRow(h3("Heatmap of selected TopGo entry")),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_down_topgo")))))
          ),
          tabPanel("UPDOWN", icon = icon("up-down"),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrUPDOWN"), "Perform gene set enrichment analysis on the up- and downregulated genes", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrUPDOWN_goseq"), "Perform gene set enrichment analysis on the up- and downregulated genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrUPDOWN_topgo"), "Perform gene set enrichment analysis on the up- and downregulated genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_updown")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_updown_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_updown_topgo")),
                   downloadButton(nsFA("downloadGOTbl_updown"), "Download", class = "btn btn-success"),
                   fluidRow(h3("Heatmap of selected TopGo entry")),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_updown_topgo")))))
          ),
          tabPanel("List1", icon = icon("list"),
                   fileInput(inputId = nsFA("gl1"),
                             label = "Upload a gene list file",
                             accept = c("text/csv", "text/comma-separated-values",
                                        "text/tab-separated-values", "text/plain",
                                        ".csv", ".tsv"), multiple = FALSE),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST1"), "Perform gene set enrichment analysis on the genes in list1", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST1_goseq"), "Perform gene set enrichment analysis on the list1 genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST1_topgo"), "Perform gene set enrichment analysis on the list1 genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list1")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list1_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_list1_topgo")),
                   downloadButton(nsFA("downloadGOTbl_l1"), "Download", class = "btn btn-success"),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_list1_topgo")))))
          ),
          tabPanel("List2", icon = icon("rectangle-list"),
                   fileInput(inputId = nsFA("gl2"),
                             label = "Upload a gene list file",
                             accept = c("text/csv", "text/comma-separated-values",
                                        "text/tab-separated-values", "text/plain",
                                        ".csv", ".tsv"), multiple = FALSE),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST2"), "Perform gene set enrichment analysis on the genes in list2", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST2_goseq"), "Perform gene set enrichment analysis on the list2 genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST2_topgo"), "Perform gene set enrichment analysis on the list2 genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list2")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list2_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_list2_topgo")),
                   downloadButton(nsFA("downloadGOTbl_l2"), "Download", class = "btn btn-success"),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_list2_topgo")))))
          ),
          tabPanel("List3", icon = icon("rectangle-list"),
                   fileInput(inputId = nsFA("gl3"),
                             label = "Upload a gene list file",
                             accept = c("text/csv", "text/comma-separated-values",
                                        "text/tab-separated-values", "text/plain",
                                        ".csv", ".tsv"), multiple = FALSE),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST3"), "Perform gene set enrichment analysis on the genes in list3", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST3_goseq"), "Perform gene set enrichment analysis on the list3 genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST3_topgo"), "Perform gene set enrichment analysis on the list3 genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list3")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list3_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_list3_topgo")),
                   downloadButton(nsFA("downloadGOTbl_l3"), "Download", class = "btn btn-success"),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_list3_topgo")))))
          ),
          tabPanel("List4", icon = icon("rectangle-list"),
                   fileInput(inputId = nsFA("gl4"),
                             label = "Upload a gene list file",
                             accept = c("text/csv", "text/comma-separated-values",
                                        "text/tab-separated-values", "text/plain",
                                        ".csv", ".tsv"), multiple = FALSE),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST4"), "Perform gene set enrichment analysis on the genes in list4", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST4_goseq"), "Perform gene set enrichment analysis on the list4 genes - goseq", class = "btn btn-primary"))),
                   fluidRow(column(width = 6, actionButton(nsFA("button_enrLIST4_topgo"), "Perform gene set enrichment analysis on the list4 genes - topGO", class = "btn btn-primary"))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list4")))),
                   fluidRow(column(width = 12, uiOutput(nsFA("ui_DT_gse_list4_goseq")))),
                   uiOutput(nsFA("ui_DT_gse_list4_topgo")),
                   downloadButton(nsFA("downloadGOTbl_l4"), "Download", class = "btn btn-success"),
                   fluidRow(column(width = 12, (uiOutput(nsFA("goterm_heatmap_list4_topgo")))))
          )
        )
      ),
      
      verbatimTextOutput(nsFA("debuggls")),
      
      h2("Intersection of gene sets"),
      
      fluidRow(
        column(width = 4,
               checkboxInput(nsFA("toggle_updown"), "Use up and down regulated genes", TRUE),
               checkboxInput(nsFA("toggle_up"), "Use up regulated genes", FALSE),
               checkboxInput(nsFA("toggle_down"), "Use down regulated genes", FALSE)
        ),
        column(width = 4,
               checkboxInput(nsFA("toggle_list1"), "Use list1 genes", TRUE),
               checkboxInput(nsFA("toggle_list2"), "Use list2 genes", FALSE),
               checkboxInput(nsFA("toggle_list3"), "Use list3 genes", FALSE),
               checkboxInput(nsFA("toggle_list4"), "Use list4 genes", FALSE)
        )
      ),
      
      fluidRow(
        column(width = 6, plotOutput(nsFA("vennlists")),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsFA("download_plot_vennlists"), "Download Plot"),
                   textInput(nsFA("filename_plot_vennlists"), label = "Save as...", value = "plot_vennlists.pdf")),
               offset = 3)),
      fluidRow(
        column(width = 6, plotOutput(nsFA("upsetLists")),
               div(align = "right", style = "margin-right:15px; margin-bottom:10px",
                   downloadButton(nsFA("download_plot_upsetlists"), "Download Plot"),
                   textInput(nsFA("filename_plot_upsetlists"), label = "Save as...", value = "plot_upsetlists.pdf")),
               offset = 3)),
      box(
        title = "Gene lists", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12, collapsed = TRUE,
        id = "gll_lists",
        fluidRow(
          column(width = 12,
                 verbatimTextOutput(nsFA("debuglists"))
          )
        ),
        fluidRow(
          column(width = 12,
                 DT::dataTableOutput(nsFA("debugTable"))
          )
        ),
        fluidRow(
          column(width = 12,
                 verbatimTextOutput(nsFA("debugTableSelected")))
        )
      )
    ),
    conditionalPanel(
      condition = "output.checkresu",
      h2("You did not create the result object yet. Please go the dedicated tab and generate it")
    )
  )
}
