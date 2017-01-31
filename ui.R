
library(shiny)

shinyUI(
    tabsetPanel(id="Unsup",
                tabPanel(title="PCA animation",
                         fluidPage(
                           fluidRow(
                             column(12, align="left",
                                    sliderInput("animation", "Probe set variance cutoff:", 0.3, 5, 0.3,step = 0.1, animate=animationOptions(interval=500, loop=FALSE)))
                           ),
                           fluidRow(
                             column(6,plotOutput("vardensity",height = 350, width = 420)),
                             column(6,plotOutput("percentVar",height = 350, width = 420))
                           ),
                           fluidRow(
                             column(6,plotOutput("pc1pc2",height = 350, width = 420)),
                             column(6,plotOutput("pc1pc3",height = 350, width = 420))
                           )
                         )),
                tabPanel(title="Probe set table",
                         fluidPage(
                           fluidRow(
                             column(12,DT::dataTableOutput("geneTable"))
                           )
                         )
                ),
                tabPanel(title="Heat map",
                         fluidPage(
                           fluidRow(
                             column(4, align="left",uiOutput("slideHeatmap")),
                             column(2, align="center",verbatimTextOutput("numGene")),
                             column(6,align="center",imageOutput("legend"))
                           ),
                           fluidRow(
                             column(12,align="center",imageOutput("heatmap"))
                           ),
                           tags$style(type='text/css', "#numGene { width:100%; margin-top: 25px;}")
                         )
                )
    )
)
    