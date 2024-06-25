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

# Read bam file paths
bamF1<- scan(bamFilePaths, what="character", sep='\n')

inBam<- unique(c(bamF1))

# This assumes that the sample name is the bam file path without the 
# folder and subfolder names and the .bam in the end of the file name. 
names(inBam)<- gsub("(.*/|_.*)","",inBam)


# Create folder to store the result files
# fro each bam file a result file is produced.

 dir.create("./IntEREst_res/")
 setwd("./IntEREst_res/")
# The methods "ExEx" and "ExSkip" are used for interest()
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
      outFile=paste(names(inBam)[i], "EEES.tsv", sep="_"),
      bamFile=bf,
      reference= refUncol,
      junctionReadsOnly=TRUE,
      method=c("ExEx","ExSkip"),
      scaleLength= c(TRUE,TRUE),
      scaleFragment=c(TRUE,TRUE),
      limitRanges=limitRange,
      loadLimitRangesReads=FALSE,
      excludeFusionReads=TRUE,
      logFile=paste("tmpFile_reverse_", names(inBam)[i],".txt", sep=""),
      strandSpecific="reverse"
    )
    
}
