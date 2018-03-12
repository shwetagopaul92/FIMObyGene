
#.libPaths("/udd/reshg/R/x86_64-pc-linux-gnu-library/3.4")
########################
# Load the libraries
########################

library(shiny)
library(GenomicRanges)
library(erma)
library(data.table)
library(DT)
library(Sushi)
library(Gviz)
library(biovizBase)
library(ensembldb)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(AnnotationFilter)


########################
# Load the data
########################

data(chr)
data(tf)
data(named_tf)
data(named_metadata_tf)

########################
# Functions
########################

tfDrill <- function(tfname, chromosome){
  setwd(paste0("~/Desktop/",tfname,"/"))
  chrbed = fread(paste0(chromosome,".bed"))
  chr = rename(chrbed, c("V1"= "chr","V2"="start","V3"="end","V4"="interval","V5"="score","V6"="strand","V7"="pvalue"))
  tfGRanges = GRanges(chr$chr, IRanges(chr$start, chr$end))
  strand(tfGRanges) = c(chr$strand)
  mcols(tfGRanges) = DataFrame(score=chr$score, pvalue=chr$pvalue)
  return(tfGRanges)
}

geneModel <- function(geneName){
  geneM = genemodel(geneName)
  return(geneM)
}

overlaps <- function(tfGRanges, geneM){
  output = subsetByOverlaps(tfGRanges, geneM)
  return(output)
}

plotGeneModel<-function(symbol){
  
  require(Sushi)
  require(erma)
  mygene = genemodel(symbol)
  
  chrom = as.character(seqnames(mygene)[1])
  chromstart = min(start(mygene))-250
  chromend = max(end(mygene))+250
  mybed = data.frame(chrom=chrom, start=start(mygene),  stop=end(mygene), gene=as.character(mygene$SYMBOL), score=".",  strand=as.numeric(paste0(as.character(strand(mygene)),1)),stringsAsFactors=F)
  
  pg = plotGenes(mybed, chrom, chromstart, chromend ,maxrows=1, bheight=0.2,
                 plotgenetype="arrow", bentline=FALSE,
                 labeloffset=.4, fontsize=1.2, arrowlength = 0.025,
                 labeltext=TRUE)
  labelgenome(chrom, chromstart, chromend, n=8, scale="Mb")
}

txmodels = function(sym, gr,
                    edb="EnsDb.Hsapiens.v75",
                    plot.it=FALSE, radius=0, ...) {
  if (sum(c(missing(sym), missing(gr)))!=1) stop("exactly one of sym, gr must be nonmissing")
  require(edb, character.only=TRUE) # FIXME
  if (!missing(sym)) {
    data(genesymbol, package="biovizBase")
    gr = genesymbol[sym]
  }
  seqlevelsStyle(gr) = "Ensembl"
  seq = as.character(seqnames(gr))[1]
  na = SeqNameFilter(c(seq))
  gf = GeneIdFilter("ENS","startsWith")
  filter = AnnotationFilter::AnnotationFilterList(na,gf)
  eg <- ensembldb::exonsBy(get(edb), by = "tx",filter = filter,columns = c("gene_biotype", "gene_name", "tx_id"))
  # names(eg) = make.names(names(eg), unique=TRUE)
  # stopifnot(sum(duplicated(names(eg)))==0)
  eg2 = as.data.frame(subsetByOverlaps(eg, gr+radius))
  eg2$chromosome = eg2$seqnames
  eg2$symbol = eg2$gene_name
  eg2$exon = eg2$exon_id
  eg2$transcript = eg2$transcript_id
  eg2$feature=eg2$gene_biotype
  eg2$gene = eg2$gene_name
  eg2$transcript=eg2$tx_id
  grt = GeneRegionTrack(eg2, ...)
  if (plot.it) plotTracks(list(GenomeAxisTrack(),grt), showId=TRUE)
  invisible(grt)
}

plotTF<-function(mytf, mysymbol){
  require(Gviz)
  mygene = txmodels(sym=mysymbol, genome="hg19")
  mychr = seqnames(mygene)[1]
  myfile = paste0("~/Desktop/",mytf,"/",mychr,".bed")
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
  tfp = plotTracks(list(gtrack, biomtrack, dtrack))
}


########################
# Define the server logic
########################

shinyServer(function(input, output, session) {
  
  output$mytable1 = renderDataTable({
    grGene = geneModel(input$geneName)
    chromosome = grGene@seqnames@values
    grTf = tfDrill(input$transcriptionFactor, chromosome)
    outputOverlaps = overlaps(grTf, grGene)
    overlaps.df = as.data.frame(outputOverlaps)
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$downloadName, Sys.Date(), ".csv", sep="")},
      content = function(file) {
        write.csv(overlaps.df, file, row.names=F)
      })
    overlaps.df
  })
 
  observeEvent(input$geneModel, {
    updateTabsetPanel(session, "inTabset", selected ="panel2")
    output$geneplot = renderPlot(plotGeneModel(input$geneName))
  })
  
  observeEvent(input$tfmodel, {
    updateTabsetPanel(session, "inTabset", selected = "panel3")
    output$tfplot = renderPlot(plotTF(input$transcriptionFactor, input$geneName))
  })

  output$mytable2 = renderDataTable(as.data.frame(metadata_tf), options=list(scrollX=TRUE,pageLength=25))
  
})