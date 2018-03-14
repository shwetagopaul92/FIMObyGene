######## Shiny APP ########
### FIMObyGene ###

#'shiny interface to Kimbie's CISBP FIMO scan TF bed files
#'@import shiny
#'@param ip character string with inet ip address
#'@export
runFIMObyGene<-function(ip){
  runApp('/udd/reshg/github/FIMObyGene/inst/shinyApps/', host=ip)
}
