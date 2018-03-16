######## Shiny APP ########
### FIMObyGene ###

#'shiny interface to Kimbie's CISBP FIMO scan TF bed files
#'@import shiny
#'@param host character string with inet ip address
#'@export
runFIMObyGene<-function(host){
  myapp = system.file("shinyApps",package="FIMObyGene")
  runApp(myapp, host=host)
}
