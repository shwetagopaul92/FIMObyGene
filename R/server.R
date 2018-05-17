########################
# Load the data
########################

data(named_tf)
data(named_metadata_tf)

plotTF<-function(mytf, mysymbol){
  require(Gviz)
  require(TFutils)
  mygene = genemodForGviz(mysymbol, resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  #mygene = txmodels(sym=mysymbol, genome="hg19")
  mychr = seqnames(mygene)[1]
  myfile = paste0("/udd/reshg/chrfiles_tf/",mytf,"/",mychr,".bed")
  mydat = read.csv(myfile, sep="\t",header=FALSE)
  mydat = mydat[,c(1,2,3,5,6)]
  names(mydat) = c("chrom","start","end","score","strand")
  mydatrange = GRanges(mychr, IRanges(start=mydat$start,end=mydat$end), mcols=mydat[,4:5])
  mygeneRange = GRanges(mychr, IRanges(start=min(start(mygene)),end=max(end(mygene))))
  mysub = subsetByOverlaps(mydatrange, mygeneRange)
  rm(mydat)
  rm(mydatrange)
  genome(mysub)="hg19"
  atrack1 <- AnnotationTrack(mysub, name = mytf)
  gtrack <- GenomeAxisTrack()
  #itrack <- IdeogramTrack(genome = "hg19", chromosome = mychr)
  dtrack=DataTrack(data=mysub$mcols.score, start=start(mysub), end=end(mysub),
                   chromosome = mychr, genome = "hg19", name = mytf)
  biomtrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = mychr, symbol=mysymbol, name=mysymbol)
  #tfp = plotTracks(list(gtrack, biomtrack, dtrack))
  tfp = plotTracks(list(dtrack, gtrack, biomtrack), sizes=c(2,2,5))
}


########################
# Define the server logic
########################

shinyServer(function(input, output, session) {
  
  output$mytable1 = renderDataTable({
    #grGene = geneModel(input$geneName)
    grGene=genemodel(input$geneName)
    chromosome = grGene@seqnames@values
    importFIMO_local_split(input$transcriptionFactor, chromosome)
    chr = plyr::rename(chrbed, c("V1"= "chr","V2"="start","V3"="end","V4"="interval","V5"="score","V6"="strand","V7"="pvalue"))
    tfGRanges = GRanges(chr$chr, IRanges(chr$start, chr$end))
    strand(tfGRanges) = c(chr$strand)
    mcols(tfGRanges) = DataFrame(score=chr$score, pvalue=chr$pvalue)
    #grTf = tfDrill(input$transcriptionFactor, chromosome)
    outputOverlaps = subsetByoverlaps(tfGRanges, grGene)
    overlaps.df = as.data.frame(outputOverlaps)
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$downloadName, Sys.Date(), ".csv", sep="")},
      content = function(file) {
        write.csv(overlaps.df, file, row.names=F)
      })
    overlaps.df
  })
  
  observeEvent(input$tfmodel, {
    updateTabsetPanel(session, "inTabset", selected = "panel3")
    output$tfplot = renderPlot(plotTF(input$transcriptionFactor, input$geneName))
  })
  
  output$mytable2 = renderDataTable(as.data.frame(metadata_tf), options=list(scrollX=TRUE,pageLength=25))
  
})