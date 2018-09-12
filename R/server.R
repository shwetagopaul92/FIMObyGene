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
    geom_segment(aes(x = st, y = yv, xend = en, yend = yv, colour = sym), data = newdf, arrow=arrow(ends=ardir, length=unit(arrmm, "mm")))
  #ggplotly(pl)
  pl + xlab(as.character(seqnames(exs)[1]))
}

enc690ByFactor = function (factor=myTF,sym=mysymbol,filtrange=NULL) 
{
  data(encode690)
  encode690 = as.data.frame(encode690)
  stopifnot(factor %in% encode690$target)
  tmp = dplyr::filter(encode690, target == factor)[, c("AHID", "cell")]
  ids = tmp[["AHID"]]
  cells = tmp[["cell"]]
  message(paste("retrieving", length(ids), "AnnotationHub elements."))
  suppressWarnings({
    suppressMessages({
      l1 = lapply(seq_len(length(ids)), function(x) {
        cat(".")
        tmp = AnnotationHub::AnnotationHub()[[ids[x]]]
        if (!is.null(filtrange)) {
          seqlevelsStyle(filtrange) = seqlevelsStyle(tmp)
          tmp = subsetByOverlaps(tmp, filtrange)
        }
        metadata(tmp) = c(metadata(tmp), cell = cells[x])
        tmp
      })
    })
  })
  cat("\n")
  names(l1) = ids
  sapply(l1,length)
  #gm = GRanges("chr17", IRanges(38077296, 38083884))
  #l2 = lapply(l1, function(x) subsetByOverlaps(x, gm+10000))
  #l3 = enc690ByFactor(filtrange=gm+10000)
  mygene = genemodelDF(sym, EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = paste0("chr",mygene$seqnames[[2]])
  br = GRanges(mychr, IRanges(start=min(mygene$start),end=max(mygene$end)))
  #br = GRanges("chr13", IRanges(32889617, 32973809)) # BRCA2
  dd = lapply(l1, function(x) subsetByOverlaps(x, br+10000))
  lens = sapply(dd,length)
  cls = sapply(dd, function(x) metadata(x)$cell)
  cls = rep(cls, lens)
  ee = do.call(rbind, lapply(dd, as.data.frame))
  ee$cell = factor(as.character(cls))
  ee$yval = 1+(as.numeric(factor(as.character(cls)))-1)/length(unique(cls))
  edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  ggvisForSymbol(sym, resource=edb) +
    geom_segment(aes(x=start, xend=end, y=yval, yend=yval,
                     group=cell, colour=cell), data=ee, size=2.5) +
    theme(axis.text.y = element_blank(), axis.title.y=element_blank()) + 
    ylim(-.5,2) + ggtitle(paste0(factor," binding near ", sym))
  
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
  
  output$tfgenePlot = renderPlot(enc690ByFactor(input$encodeTF, input$geneName))
  
  require(plotly)
  output$ggPlot = renderPlotly(ggvisForSymbol(input$geneName))
  
  output$mytable2 = renderDataTable(as.data.frame(metadata_tf), options=list(scrollX=TRUE,pageLength=25))
  
})
