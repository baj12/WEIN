tags_head_ui <- function(){
  tags$head(
    tags$style(HTML("
                        .shiny-output-error-validation {
                        font-size: 15px;
                        color: forestgreen;
                        text-align: center;
                        }
                        .icon-done {
                        color: green;
                        }
                        #myScrollBox{
                        overflow-y: scroll;

                        .dataTables_wrapper{
                        overflow-x: scroll;
                        }
                        }
                        #myAnchorBox{}
                        "))
  )
}