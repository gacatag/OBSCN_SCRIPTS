# It is assumed that the working directory is where the files of the scripts are
workDir= "/PATH/TO/WORKFOLDER/AND/SCRIPTS/FILES"
# bamFiles.txt should incldue the paths of the input bam files
bamFilePaths="./bamFiles.txt"
# Number of computing cores to use for analysis
ncore=20

setwd(workDir)

load(file="./refUncol.rda")
load(file="./limitRange.rda")


library(BiocParallel)
library(IntEREst)

#Prepare bamfiles
bamF1<-scan(bamFilePaths, what="character", sep='\n')
inBam<- unique(bamF1)

names(inBam)<- gsub("(.*/|_.*)","",inBam)


dir.create("./IntEREst_res/")
setwd("./IntEREst_res/")

# The methods "IntSpan" are used for interest()
for(i in 1:length(inBam) ){
  
  print(paste(i, length(inBam), sep="/"))
  bf<- as.character(inBam[i])
  print(paste(i, length(inBam), sep="/"))
  
  interest(
    bamFileYieldSize=100000,
    bpparam=SnowParam(workers = ncore, type = "SOCK"),
    isPaired=TRUE,
    isPairedDuplicate=NA,
    isSingleReadDuplicate=NA,
    referenceGeneNames=refUncol$transcript_id,
    referenceIntronExon=refUncol$int_ex,
    outFile=paste(names(inBam)[i], "IS.tsv", sep="_"),
    bamFile=bf,
    reference= refUncol,
    junctionReadsOnly=TRUE,
    method=c("IntSpan"),
    scaleLength= c(TRUE, TRUE),
    scaleFragment=c(TRUE, TRUE),
    limitRanges=limitRange,
    loadLimitRangesReads=FALSE,
    excludeFusionReads=TRUE,
    logFile=paste("tmpFile_reverse_", names(inBam)[i],".txt", sep=""),
    strandSpecific="reverse"
  )
  
}



