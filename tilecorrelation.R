library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(IRanges)

#make GRanges object from BSgenome (accessed with "Celegans") (i remove the mito chr)
genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6],ranges=IRanges(start=1, end=seqlengths(Celegans)[1:6]), strand="*")

# get a list of all the files with the enrichment data
Datafiles<-list.files("/Users/imac/Desktop/Bolaji/PhD/analysis/Bioinformatics_scripts/Bigwiggle/","bw")

#create GRangesList object with different sized tiles along genome
tile1Mb<-unlist(tile(x=genomeGR,width=1000000))
tile100kb<-unlist(tile(x=genomeGR,width=100000))
tile10kb<-unlist(tile(x=genomeGR,width=10000))
tile1kb<-unlist(tile(x=genomeGR,width=1000))
tile100bp<-unlist(tile(x=genomeGR,width=100))
tileList<-list("tile1Mb"=tile1Mb,"tile100kb"=tile100kb,"tile10kb"=tile10kb,"tile1kb"=tile1kb,"tile100bp"=tile100bp)

# read in the bedgraph files in pairs to calculate enrichment counts at different scales
for (f in 1:length(Datafiles)) {
  #get normalised counts
  Data<-import(paste0("/Users/imac/Desktop/Bolaji/PhD/analysis/Bioinformatics_scripts/Bigwiggle/",Datafiles[f]))

  #make a GRangesList from GRanges
  for (i in 1:length(tileList)) {
    tiles<-split(tileList[[i]],seqnames(tileList[[i]]))
    
    # convert GRanges back to coverage
    rle<-coverage(Data,weight="score")
    if (length(rle)>6){
      j<-grep("M",names(rle))
      rle<-coverage(Data,weight="score")[-j]
    }
    
    names(rle)<-gsub("chr","", names(rle))
    names(rle)<-paste0("chr",names(rle))
    #seqlevels(rle)<-paste0("chr", seqlevels(rle))
    
    #create views corresponding to the tiles on the coverage rle and average counts
    summedCov<-viewMeans(RleViewsList(rangesList=tiles,rleList = rle))
    
    #now save the values into tiles
    mcols(tileList[[i]])[,Datafiles[f]]<-as.data.frame(unlist(summedCov))
  }
}

BiocManager::install("GGally")
ggpairs(as.matrix(mcols(tileList[[i]])))

