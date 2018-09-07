########################
# Load the data
########################

data(named_tf)
data(named_metadata_tf)
data(encode690)

########################
# Define UI for the application
########################

shinyUI(fluidPage(
  
  titlePanel("FIMObyGene"),
  sidebarLayout(position = "left",
                sidebarPanel(width=3,
                             selectInput("transcriptionFactor", "Select Transcription Factor", named_tf),
                             selectInput("encodeTF", "Select Encode Transcription Factor",encode690$target ),
                             textInput("geneName", "Enter gene of interest", value="ORMDL3"),
                             textInput("downloadName", "Download Name"),
                             downloadButton("downloadData", "Download Data")
                ),
                
                
                mainPanel(
                  tabsetPanel(id="inTabset",
                              tabPanel("Scored Motifs in Transcribed Region", value="panel1", DT::dataTableOutput("mytable1")),
                              tabPanel("TF Model", value="panel2", plotOutput("tfplot")),
                              #tabPanel("Gene Model", value="panel3", plotlyOutput("ggPlot")),
                              tabPanel("GeneTF Model", value="panel4", plotOutput("tfgenePlot")),
                              tabPanel("Metadata",value="panel3", DT::dataTableOutput("mytable2"))
                  )
                )
  )
)
)
