# It is assumed that the working directory is where the files of the scripts are
workDir= "/PATH/TO/WORKFOLDER/AND/SCRIPTS/FILES"
refAnnoGtf<- 'PATH/TO/GENE/ANNOTATION/GTF/FILE'


# Running on SERVER hackman-2-20
setwd(workDir)


library(BiocParallel)
library(IntEREst)
library(utils)
library(R.utils)

refUncol<- referencePrepare (sourceBuild="file", filePath= refAnnoGtf, fileFormat="gtf", collapseExons=FALSE, ignore.strand=FALSE)
unique(refUncol$chr)
unique(refUncol$strand)
# Make reference only from genes on Chr1 to Chr22 and chrX, chrY and chrM.
refUncol<- refUncol[refUncol$chr%in%paste("chr",c(1:22,"X","Y","M"), sep=""),]
save(refUncol, file="refUncol.rda")


selectTr<- unique(refUncol$transcript_id)
minSt<- tapply(refUncol$begin, refUncol$transcript_id, min)
maxEn<- tapply(refUncol$end, refUncol$transcript_id, max)
uniCh<- tapply(refUncol$chr, refUncol$transcript_id, head, n=1)
limitRange<- GRanges(
  seqnames= as.character(unlist(uniCh)),
  IRanges::IRanges(
    start=as.numeric(unlist(minSt))-5 ,
    end=as.numeric(unlist(maxEn))+5
  ))
save(limitRange, file="limitRange.rda")


refCol<- referencePrepare (sourceBuild="file", filePath= refAnnoGtf, fileFormat="gtf", collapseExons=TRUE, ignore.strand=FALSE)
refCol<- refCol[refCol$chr%in%paste("chr",c(1:22,"X","Y","M"), sep=""),]
save(refCol, file="refCol.rda")






