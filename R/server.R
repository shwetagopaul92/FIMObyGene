
########################
# Functions
########################

plotTF<-function(mytf, mysymbol){
  require(Gviz)
  require(TFutils)
  require(GenomicRanges)
  mygene = genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = seqnames(mygene)[1]
  tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",mytf,".02_sort.bed.gz"))
  mygeneRange = GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mysub = importFIMO(tf, mygeneRange)
  genome(mysub)="hg19"
  atrack1 <- AnnotationTrack(mysub, name = mytf)
  gtrack <- GenomeAxisTrack()
  dtrack=DataTrack(data=mysub$score, start=start(mysub), end=end(mysub),
                   chromosome = mychr, genome = "hg19", name = mytf)
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  tfp = plotTracks(list(dtrack, gtrack, biomtrack), sizes=c(2,2,5))
}


########################
# Define the server logic
########################

shinyServer(function(input, output, session) {
  
  output$mytable1 = renderDataTable({
    require(GenomicRanges)
    require(TFutils)
    require(Rsamtools)
    genemodel.df = genemodelDF(input$geneName, EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
    grGene = makeGRangesFromDataFrame(genemodel.df)
    seqlevels(grGene) <- sub("","chr",seqlevels(grGene))
    chromosome = grGene@seqnames@values
    tf = Rsamtools::TabixFile(paste0("/udd/reshg/tbifiles/tabix_all_tf_new/",input$transcriptionFactor,".02_sort.bed.gz"))
    outputOverlaps = importFIMO(tf, range(grGene))
    overlaps.df = as.data.frame(outputOverlaps)
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$downloadName, Sys.Date(), ".csv", sep="")},
      content = function(file) {
        write.csv(overlaps.df, file, row.names=F)
      })
    overlaps.df
  })
  
  observeEvent(input$tfmodel, {
    updateTabsetPanel(session, "inTabset", selected = "panel2")
    output$tfplot = renderPlot(plotTF(input$transcriptionFactor, input$geneName))
  })
  
  output$mytable2 = renderDataTable(as.data.frame(metadata_tf), options=list(scrollX=TRUE,pageLength=25))
  
})
