measurePSIUnadj<- function(dataPath,
                outPath,
                uniqExGr){
# Running on SERVER hackman-2-20
setwd(outPath)

  for(f in dir("./","\\.rda$", full.names = TRUE))
    load(f)
  for(f in dir("./shinyData","\\.rda$", full.names = TRUE))
    load(f)

resFilesPath<- paste(dataPath,"/res/",sep="/")



esFiles<- dir(resFilesPath, "_EEES.tsv")
esNames<- gsub("_EEES.tsv","", esFiles)
names(esFiles)<- esNames

eeFiles<- esFiles
eeNames<- esNames

table(eeNames==esNames)
# TRUE 
# 45 

library(openxlsx)
anno<- read.xlsx("./annoFile.xlsx", 1)


table(anno$sample_ID%in%eeNames)
# FALSE  TRUE 
# 36     7 

if(!all(anno$sample_ID%in%eeNames))
  stop(paste(
    "The sample names in sample_ID column don't match the _EEES.tsv files in", 
             paste(dataPath,"/res/",sep="/")))

indMatch<-lapply(anno$sample_ID, function(x) {return(which( (eeNames==x) | (eeNames==paste("FM",x, sep=""))))})

group= anno$group
names(eeFiles)<- eeNames
names(esFiles)<- esNames

# Reorder files to match the excel sample order
eeFiles<- eeFiles[match(anno$sample_ID, eeNames)]
esFiles<- esFiles[match(anno$sample_ID, esNames)]

# Read batch information
batchVec<- scan("./batchVec.txt", what = "character", sep = "\n")

# Building SummarizedExperimnt objects from summarization results
library(IntEREst)
resExEx<- readInterestResults(
  resultFiles=paste(resFilesPath, as.character(eeFiles[match(as.character(mapNamesfil), names(eeFiles))]), sep=""),
  sampleNames=as.character(mapNamesfil), 
  sampleAnnotation=data.frame(sampleNames=as.character(mapNamesfil), 
                              group=group,
                              batch=batchVec), 
  commonColumns=1:8, freqCol=9, scaledRetentionCol=10, 
  scaleLength=FALSE, scaleFragment=FALSE, reScale=FALSE
)
resExSk<- readInterestResults(
  resultFiles=paste(resFilesPath,  as.character(esFiles[match(as.character(mapNamesfil), names(esFiles))]), sep=""),
  sampleNames=as.character(mapNamesfil),
  sampleAnnotation=data.frame(sampleNames=as.character(mapNamesfil), 
                              group=group,
                              batch=batchVec), 
  commonColumns=1:8, freqCol=11, scaledRetentionCol=12, 
  scaleLength=FALSE, scaleFragment=TRUE, reScale=FALSE
)


exObj<- cbind(resExEx[rowData(resExEx)$int_ex=="exon",],resExSk[rowData(resExSk)$int_ex=="exon",])

exObj<- addAnnotation(x=exObj,
                         sampleAnnotationType="EEES",
                         sampleAnnotation=factor(c(rep("ee",ncol(resExEx)), rep("es",ncol(resExSk))), levels=c("es","ee"))
)

save(resExSk, file="resExSk.rda")
save(resExEx, file="resExEx.rda")
save(exObj, file="exObj.rda")

######### Make the data for Shiny FOR UNADJUSTED PSI levels





allExObj<- exObj

exGr<- GRanges(
  seqnames= as.character(rowData(allExObj)[,"chr"]),
  IRanges::IRanges(
    start=as.numeric(rowData(allExObj)[,"begin"]) ,
    end=as.numeric(rowData(allExObj)[,"end"])
  ))

# Get result obejct with unique exons
uniqAllExGr<- unique(exGr)

uniqInd<- findOverlaps(uniqAllExGr, exGr, type="equal", select= "first")

uniqExObj<- allExObj[uniqInd,]

table((start(uniqAllExGr))==(rowData(uniqExObj)[,"begin"]))


uniqEeExObj<- uniqExObj[,which(colData(uniqExObj)$EEES=="ee")]
uniqEsExObj<- uniqExObj[,which(colData(uniqExObj)$EEES=="es")]

library(sva)


# Building annotation vector for groups of samples

groupSam<- anno$group

groupSam<- factor(groupSam, levels=c("FETAL MUSCLE","ADULT MUSCLE", "FETAL HEART", "ADULT HEART"))

save(groupSam, file="shinyData/groupSam.rda")

# Creating matrix that includes average PSIs within groups
psimean<- t(aggregate(psidat, by=list(groupSam), mean))
psimeanCol<- psimean[1,]
psimean<- apply(psimean[-1,],2,as.numeric)
colnames(psimean)<- psimeanCol

uniqExObjEnsembl<- uniqExObj
allExObjEnsembl<- allExObj

#### Select Reliable REFSEQ annotated transcript and the isoforms of OBSCN that we want to study

selectTr<- scan(trFile, what="character", sep = "\n")
selectGen<- rep(geneName, length(selectTr))
selectTr<- gsub("\\..*","", selectTr)

require('biomaRt')

mart <- useMart('ENSEMBL_MART_ENSEMBL', host="http://jul2023.archive.ensembl.org")
mart <- useDataset('hsapiens_gene_ensembl', mart)
annoTmp <- getBM(
  mart = mart,
  attributes = c(
    'ensembl_transcript_id', "ensembl_gene_id",
    "refseq_mrna", "refseq_ncrna"),
  filters = "ensembl_transcript_id",
  values = unique(gsub("\\..*","",rowData(allExObj)[,"transcript_id"])),
  uniqueRows = TRUE)

# Choosing transcripts with RefSeq ID
inTr<- unique(c(annoTmp[which(annoTmp$refseq_mrna!="" | annoTmp$refseq_ncrna!=""),"ensembl_transcript_id"], selectTr))

allExObj<- allExObj[which(gsub("\\..*","",rowData(allExObj)[,"transcript_id"])%in%inTr),]

exGr<- GRanges(
  seqnames= as.character(rowData(allExObj)[,"chr"]),
  IRanges::IRanges(
    start=as.numeric(rowData(allExObj)[,"begin"]) ,
    end=as.numeric(rowData(allExObj)[,"end"])
  ), strand=rowData(allExObj)[,"strand"])

### Uncomment next if Adjusted data
uniqAllExGr<- unique(exGr)


uniqInd<- findOverlaps(uniqAllExGr, exGr, type="equal", select= "first")

uniqExObj<- allExObj[uniqInd,]


save(psidat, file="shinyData/psidat.rda")
save(allExObjEnsembl, file="allExObjEnsembl.rda")
save(allExObj, file="allExObj.rda")
save(psimean, file="shinyData/psimean.rda")
save(uniqExObjEnsembl, file="uniqExObjEnsembl.rda")
save(uniqExObj,file="uniqExObj.rda")

}