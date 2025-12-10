signature_ui <- function(){
  nsSE <- NS("sig_expl")
  tabPanel(
    "Signatures Explorer",
    icon = icon("map"),
    conditionalPanel(
      condition="!output.checkdds",
      fluidRow(
        column(
          width = 8,
          shinyBS::bsCollapse(
            id = "help_signatureexplorer",open = NULL,
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_signatureexplorer.md",package = "WEIN")))
          )
        )
      ),
      actionButton("tour_signatureexplorer", "Click me for a quick tour of the section", icon("info"),
                   style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),
      
      fluidRow(
        column(
          width = 6,
          h4("Setup options"),
          wellPanel(
            uiOutput(nsSE("sig_ui_gmtin")),
            uiOutput(nsSE("sig_ui_nrsigs")),
            actionButton(nsSE("sig_button_computevst"),
                         label = "Compute the variance stabilized transformed data", 
                         icon = icon("spinner"), class = "btn btn-success")
          )
        ),
        column(
          width = 6,
          h4("Conversion options"),
          wellPanel(
            uiOutput(nsSE("sig_ui_id_data")),
            uiOutput(nsSE("sig_ui_id_sigs")),
            uiOutput(nsSE("sig_ui_orgdbpkg")),
            actionButton(nsSE("sig_convert_setup"),
                         label = "Apply id conversion between data and signatures")
          ),
          verbatimTextOutput(nsSE("sig_convcheck"))
          
        )
      ),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            uiOutput(nsSE("sig_ui_selectsig")),
            uiOutput("sig_ui_annocoldata"),
            checkboxInput(nsSE("sig_useDEonly"),
                          label = "Use only DE genes in the signature",value = FALSE)
          )
         ),
        column(
          width = 6,
          wellPanel(
            checkboxInput(nsSE("sig_clusterrows"),label = "Cluster rows", value = TRUE),
            checkboxInput(nsSE("sig_clustercols"), label = "Cluster columns"),
            checkboxInput(nsSE("sig_centermean"), label = "Center mean",value = TRUE),
            checkboxInput(nsSE("sig_scalerow"), label = "Standardize by row")
          )
        )
      ),
      fluidRow(
        column(
          width = 8, offset = 2,
          shinyjqui::jqui_resizable(plotlyOutput(nsSE("sig_heat")))
        )
      ),
      fluidRow(
        column(
          width = 8, offset = 2,
          verbatimTextOutput(nsSE("sig_heat_genes"))
        )
      )
    ),
    conditionalPanel(
      condition="output.checkdds",
      h2("You did not create the dds object yet. Please go the main tab and generate it")
    )
  )
}