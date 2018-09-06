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

ggvisForSymbol = function (sym, resource = EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79, 
                           columnsKept = c("gene_id", "tx_id"), yval = 1, arrmm=1.5, viewtype="transcripts", ...) 
{
  exs = GenomicFeatures::exons(resource, filter = SymbolFilter(sym), columns = columnsKept, 
                               ...)
  if (viewtype == "exons") exs = unique(exs)
  rd = reduce(exs)
  fo = findOverlaps(rd, exs)
  gr = split(subjectHits(fo), queryHits(fo))
  pp = function(n) (seq_len(n)-1)/n
  st = start(exs)
  en = end(exs)
  if (viewtype == "exons") {
    ys = lapply(gr, function(x) pp(length(x)))
    yvs = unlist(ys) #1+(0:(nel-1))/nel
  }
  else if (viewtype == "transcripts") {
    tnms = exs$tx_id
    ft = factor(tnms)
    yvs = (as.numeric(ft)-1)/length(levels(ft))
  }
  else stop("viewtype not %in% c('exons', 'transcripts')")
  newdf = data.frame(st, en, yv = yvs, sym = sym)
  rng = range(exs)
  df = data.frame(range = c(start(rng), end(rng)), yval = rep(yval,2)) 
  strn = as.character(strand(exs)[1])
  ardir = ifelse(strn=="+", "last", "first")
  pl = ggplot(df, aes(x = range, y = yval)) + 
    geom_segment(aes(x = st, y = yv, xend = en, yend = yv, colour = sym),       data = newdf, arrow=arrow(ends=ardir, length=unit(arrmm, "mm")))
  ggplotly(pl)
  #pl + xlab(as.character(seqnames(exs)[1]))
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
  
  #observeEvent(input$tfmodel, {
   # updateTabsetPanel(session, "inTabset", selected = "panel2")
  
  output$tfplot = renderPlot(plotTF(input$transcriptionFactor, input$geneName))
  require(plotly)
  output$ggPlot = renderPlotly(ggvisForSymbol(input$geneName))
  
  output$mytable2 = renderDataTable(as.data.frame(metadata_tf), options=list(scrollX=TRUE,pageLength=25))
  
})
