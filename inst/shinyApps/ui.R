
########################
# Load the data
########################

data(chr)
data(tf)
data(named_tf)
data(named_metadata_tf)

########################
# Define UI for the application
########################

shinyUI(fluidPage(
  
  titlePanel("FIMObyGene"),
  sidebarLayout(position = "left",
                sidebarPanel(width=4,
                             selectInput("transcriptionFactor", "Select Transcription Factor", named_tf),
                             textInput("geneName", "Enter gene of interest"),
                             actionButton("geneModel", "Gene Model"),
                             #actionButton("subsetByoverlaps", "Subset By Overlaps"),
                             actionButton("tfmodel", "TF Model"),
                             textInput("downloadName", "Download Name"),
                             downloadButton("downloadData", "Download Data")
                ),
                
                
                mainPanel(
                  tabsetPanel(id="inTabset",
                              tabPanel("Subset By Overlaps", value="panel1", DT::dataTableOutput("mytable1")),
                              tabPanel("Gene Model", value="panel2", plotOutput("geneplot")),
                              tabPanel("TF Model", value="panel3", plotOutput("tfplot")),
                              tabPanel("Metadata",value="panel4", DT::dataTableOutput("mytable2"))
                  )
                )
  )
)
)