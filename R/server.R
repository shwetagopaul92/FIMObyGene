#not uptodate
########################
# Functions
########################

#txmodels = function(sym, gr,
#                    edb="EnsDb.Hsapiens.v75",
#                    plot.it=FALSE, radius=0, ...) {
#  if (sum(c(missing(sym), missing(gr)))!=1) stop("exactly one of sym, gr must be nonmissing")
#  require(edb, character.only=TRUE)
#  if (!missing(sym)) {
#    data(genesymbol, package="biovizBase")
#     gr = genesymbol[sym]
#   }
#   seqlevelsStyle(gr) = "Ensembl"
#   seq = as.character(seqnames(gr))[1]
#   na = SeqNameFilter(c(seq))
#   gf = GeneIdFilter("ENS","startsWith")
#   filter = AnnotationFilter::AnnotationFilterList(na,gf)
#   eg <- ensembldb::exonsBy(get(edb), by = "tx",filter = filter,columns = c("gene_biotype", "gene_name", "tx_id"))
#   # names(eg) = make.names(names(eg), unique=TRUE)
#   # stopifnot(sum(duplicated(names(eg)))==0)
#   eg2 = as.data.frame(subsetByOverlaps(eg, gr+radius))
#   eg2$chromosome = eg2$seqnames
#   eg2$symbol = eg2$gene_name
#   eg2$exon = eg2$exon_id
#   eg2$transcript = eg2$transcript_id
#   eg2$feature=eg2$gene_biotype
#   eg2$gene = eg2$gene_name
#   eg2$transcript=eg2$tx_id
#   grt = GeneRegionTrack(eg2, ...)
#   if (plot.it) plotTracks(list(GenomeAxisTrack(),grt), showId=TRUE)
#   invisible(grt)
# }

plotTF<-function(mytf, mysymbol){
  require(Gviz)
  require(TFutils)
  require(GenomicRanges)
  #mygene = txmodels(sym=mysymbol, genome="hg19")
  mygene = genemodForGviz(sym=mysymbol,resource= EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  mychr = seqnames(mygene)[1]
  myfile = paste0("/udd/reshg/chrfiles_tf/",mytf,"/",mychr,".bed")
  #myfile = paste0("~/Desktop/",mytf,"/",mychr,".bed")
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
    require(GenomicRanges)
    require(TFutils)
    #grGene = genemodel(input$geneName)
    genemodel.df = genemodelDF(input$geneName, EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
    grGene = makeGRangesFromDataFrame(genemodel.df)
    seqlevels(grGene) <- sub("","chr",seqlevels(grGene))
    chromosome = grGene@seqnames@values
    #grTf = tfDrill(input$transcriptionFactor, chromosome)
    #chrbed = importFIMO_local_split(input$transcriptionFactor,chromosome)
    myfile = paste0("/udd/reshg/chrfiles_tf/",input$transcriptionFactor,"/",chromosome,".bed")
    chrbed = fread(myfile)
    chr = plyr::rename(chrbed, c("V1"= "chr","V2"="start","V3"="end","V4"="interval","V5"="score","V6"="strand","V7"="pvalue"))
    tfGRanges = GRanges(chr$chr, IRanges(chr$start, chr$end))
    strand(tfGRanges) = c(chr$strand)
    mcols(tfGRanges) = DataFrame(score=chr$score, pvalue=chr$pvalue)
    outputOverlaps = subsetByOverlaps(tfGRanges, grGene)
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
