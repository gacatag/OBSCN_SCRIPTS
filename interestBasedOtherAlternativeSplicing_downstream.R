library(IntEREst)
workDir= "/PATH/TO/WORKFOLDER/"
setwd(workDir)

load(file="./shinyData/uniqExGr.rda")
load(file="./shinyData/dupInd.rda")
load(file="./shinyData/firstExInd.rda")
load(file="./shinyData/groupSam.rda")
load(file="./shinyData/lastExInd.rda")
load(file="./allExObj.rda")
load(file="./shinyData/psidat.rda")
load(file="./shinyData/difMuscleSel.rda")
load(file="./shinyData/difAdultSel.rda")
load(file="./shinyData/difListSel.rda")
load(file="./shinyData/psimean.rda")
load("./exObj.rda")


exObjGr<- GRanges(  seqnames= as.character(rowData(exObj)$chr),
                         IRanges::IRanges(
                           start=as.numeric(rowData(exObj)$begin),
                           end=as.numeric(rowData(exObj)$end)
                         ), strand=as.character(rowData(exObj)$strand))
uniqExOrig<- unique(exObjGr)
uniqExOrig$strand="+"



exGr<- GRanges(
  seqnames= as.character(rowData(exObj)[,"chr"]),
  IRanges::IRanges(
    start=as.numeric(rowData(exObj)[,"begin"]) ,
    end=as.numeric(rowData(exObj)[,"end"])
  ), strand= as.character(rowData(exObj)[,"strand"]))



uniExIndEqu<- findOverlaps(unique(exGr[rowData(exObj)[,"int_ex_num"]==1,]), exGr, type="equal")
uniExIndAny<- findOverlaps(unique(exGr[rowData(exObj)[,"int_ex_num"]==1,]), exGr, type="any")

altFirst<- unique(exGr[subjectHits(uniExIndAny)[!(subjectHits(uniExIndAny)%in%subjectHits(uniExIndEqu))]])

uniExIndChk<- uniExIndAny[
  !(apply(as.data.frame(uniExIndAny), 1, paste, collapse="/") %in%
      apply(as.data.frame(uniExIndEqu), 1, paste, collapse="/")),
  ]



allFirst<- unique(exGr[rowData(exObj)[,"int_ex_num"]==1,])

allFirst<- allFirst[findOverlaps(uniqExGr[firstExInd,],allFirst, type="equal", select="first"),]

(allFirst<- allFirst[order(start(allFirst), decreasing = F),])

altFirst<- altFirst[order(start(altFirst),decreasing=F),]
altFirst



uniExIndEqu<- findOverlaps(unique(exGr[c(which(rowData(exObj)[,"int_ex_num"]==1)[-1]-1,length(exGr)),]), exGr, type="equal")
uniExIndAny<- findOverlaps(unique(exGr[c(which(rowData(exObj)[,"int_ex_num"]==1)[-1]-1,length(exGr)),]), exGr, type="any")

altLast<- unique(exGr[subjectHits(uniExIndAny)[!(subjectHits(uniExIndAny)%in%subjectHits(uniExIndEqu))]])


uniExIndChk<- uniExIndAny[
  !(apply(as.data.frame(uniExIndAny), 1, paste, collapse="/") %in%
      apply(as.data.frame(uniExIndEqu), 1, paste, collapse="/")), 
]

rowData(exObj)[subjectHits(uniExIndChk),]

uniqExGr[lastExInd,]


allLast<- unique(exGr[as.numeric(tapply(1:length(exGr), rowData(exObj)[,"transcript_id"], function(x) {return(x[which(rowData(exObj)[x,"int_ex_num"]==max(rowData(exObj)[x,"int_ex_num"]))])})),])

# The 126 th exon was ommited from the exon inclusion plot because its start position is similar to the 125th exon (alo a last exon)
load("./obscnUniExGr.rda")
load("./omittedExon.rda")

allLast<- allLast[findOverlaps(obscnUniExGr[c(lastExInd,omittedExon),],allLast, type="equal", select="first"),]


(allLast<- allLast[order(start(allLast), decreasing = F),])


  

save(allLast, file="allLast.rda")
save(allFirst, file="allFirst.rda")
save(altLast, file="altLast.rda")
save(altFirst, file="altFirst.rda")
### Now Introin Spanning results

resFilesPath="./IntEREst_res/"
isFilesPath="./IntEREst_res/"

isFiles<- dir(isFilesPath, "_IS.tsv")
isNames<- gsub("_IS.tsv","", isFiles)
names(isFiles)<- isNames


library(openxlsx)
anno<- read.xlsx("./annoFile.xlsx", 1)


if(!all(anno$sample_ID%in%isNames))
  stop(paste(
    "The sample names in sample_ID column don't match the _IS.tsv files in", 
    paste(dataPath,"/res/",sep="/")))


indMatch<-match(anno$sample_ID, isNames)
isFiles<- isFiles[indMatch]
isNames<- isNames[indMatch]
group<- anno$group


names(isFiles)<- isNames

library(IntEREst)
resIntSp<- readInterestResults(
  resultFiles=paste(isFilesPath, as.character(isFiles[match(as.character(mapNamesfil), names(isFiles))]), sep=""),
  sampleNames=as.character(mapNamesfil), 
  sampleAnnotation=data.frame(sampleNames=as.character(mapNamesfil), group=group), 
  commonColumns=1:8, freqCol=9, scaledRetentionCol=10, 
  scaleLength=FALSE, scaleFragment=FALSE, reScale=FALSE
)

isObj<-resIntSp

isObj<- isObj[rowData(isObj)$int_ex=="intron",]
save(isObj, file="isObj.rda")


allIsObj<- isObj
save(allIsObj, file="allIsObj.rda")

unique(rowData(allIsObj)$transcript_id[which(gsub("\\..*","",rowData(allIsObj)$gene_id)=="ENSG00000154358")])

allIsObjBkup<- allIsObj

trFile="./OBSCN_ENSID.txt"
selectTr<- scan(trFile, what="character", sep = "\n")
selTr<- gsub("\\..*","",selectTr)

allIsObj<- allIsObj[which((gsub("\\..*","",rowData(allIsObj)$gene_id)!="ENSG00000154358")|
                            ((gsub("\\..*","",rowData(allIsObj)$gene_id)=="ENSG00000154358") & (gsub("\\..*","",rowData(allIsObj)$transcript_id)%in%selTr))),]

unique(rowData(allIsObj)$transcript_id[which(gsub("\\..*","",rowData(allIsObj)$gene_id)=="ENSG00000154358")])

save(allIsObj, file="allIsObj.rda")

allIsObjEnsembl<- allIsObj
save(allIsObjEnsembl, file="allIsObjEnsembl.rda")

#### Select Refseq and the chosen transcripts for OBSCN 
selectTr<- scan(trFile, what="character", sep = "\n")
#selectGen<- rep(geneName, length(selectTr))
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
  values = unique(gsub("\\..*","",rowData(allIsObj)[,"transcript_id"])),
  uniqueRows = TRUE)

inTr<- unique(c(annoTmp[which(annoTmp$refseq_mrna!="" | annoTmp$refseq_ncrna!=""),"ensembl_transcript_id"], selectTr))

allIsObj<- allIsObjEnsembl[which(gsub("\\..*","",rowData(allIsObjEnsembl)[,"transcript_id"])%in%inTr),]
#allIsObj<- allIsObjBkup[which(gsub("\\..*","",rowData(allIsObjBkup)[,"transcript_id"])%in%inTr),]
save(allIsObj, file="allIsObj.rda")
dim(allIsObj)
# [1] 463501     75
dim(allIsObjEnsembl)
# [1] 1307462      75

  
batchVec<- scan("./batchVec.txt", what="character", sep="\n")


## Get info for Unique coordinates. 

allIsGr<- GRanges(as.character(rowData(allIsObj)$chr), 
                  IRanges(start=as.numeric(rowData(allIsObj)$begin),
                          end=as.numeric(rowData(allIsObj)$end)),
                  strand=as.character(rowData(allIsObj)$strand))

allIsUniqGr<- unique(allIsGr)

selRows<- findOverlaps(allIsUniqGr, allIsGr, type="equal", select="first")

uniqIsObj<- allIsObj[selRows,]
uniqIsGr<- allIsUniqGr

## Check PSI of 1st and last exons.rda
allLast
allFirst

# #### estimate Pvalues





#load(file="allIsObj.rda")
load(file="./allIsObj.rda")

require('biomaRt')
print(Sys.setenv(HTTP_PROXY= "http://www-cache.cs.helsinki.fi:3128",
                 HTTPS_PROXY= "http://www-cache.cs.helsinki.fi:3128",
                 https_proxy= "http://www-cache.cs.helsinki.fi:3128",
                 http_proxy= "http://www-cache.cs.helsinki.fi:3128"))
mart <- useMart('ENSEMBL_MART_ENSEMBL', host="http://jul2023.archive.ensembl.org")
mart <- useDataset('hsapiens_gene_ensembl', mart)
annoTr <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id',
                 'ensembl_transcript_id', "transcript_tsl", "transcript_gencode_basic" ,"transcript_mane_select", "refseq_mrna"),
  filter='ensembl_transcript_id',
  values = gsub("\\..*", "", unique(rowData(allIsObj)$transcript_id)),
  uniqueRows = TRUE)

save(annoTr, file="annoTr.rda")


firstInd<- c( which(rowData(allIsObj)$strand=="+" & rowData(allIsObj)$int_ex_num==2),
               as.numeric(tapply(which(rowData(allIsObj)$strand=="-"),
                      rowData(allIsObj)$transcript_id[rowData(allIsObj)$strand=="-"],
                      max))) 
# Incorrect 1st and last exons of OBSCN were detected in DifExInc_obscn.R. It is stored in rmZeroGr
load(file="./rmZeroGr.rda")
allIsEndGr<- GRanges(seqnames(allIsGr),
                     IRanges(start=end(allIsGr)+1,
                             end=end(allIsGr)+1),
                     strand=strand(allIsGr))
allIsStaGr<- GRanges(seqnames(allIsGr),
                     IRanges(start=start(allIsGr)-1,
                             end=start(allIsGr)-1),
                     strand=strand(allIsGr))

firstInd<-firstInd[-c(queryHits(findOverlaps(allIsStaGr[firstInd,],rmZeroGr, type="end")))]
firstObj<- allIsObj[firstInd,]

lastInd<- c( which(rowData(allIsObj)$strand=="-" & rowData(allIsObj)$int_ex_num==2),
              as.numeric(tapply(which(rowData(allIsObj)$strand=="+"),
                                rowData(allIsObj)$transcript_id[rowData(allIsObj)$strand=="+"],
                                max))) 
lastInd<- lastInd[-c(queryHits(findOverlaps(allIsEndGr[lastInd,],rmZeroGr, type="start")))]
lastObj<- allIsObj[lastInd,]

firstGr<- GRanges(
  seqnames = rowData(firstObj)[,"gene_id"],
  ranges = IRanges(start=rowData(firstObj)[,"begin"],
                   end=rowData(firstObj)[,"end"]),
  strand = rowData(firstObj)[,"strand"])
lastGr<- GRanges(
  seqnames = rowData(lastObj)[,"gene_id"],
  ranges = IRanges(start=rowData(lastObj)[,"begin"],
                   end=rowData(lastObj)[,"end"]),
  strand = rowData(lastObj)[,"strand"])

firstUniGr<- unique(firstGr)
lastUniGr<- unique(lastGr)

firstIn<- names(table(seqnames(firstUniGr))[which(table(seqnames(firstUniGr))>1)])
lastIn<- names(table(seqnames(lastUniGr))[which(table(seqnames(lastUniGr))>1)])

firstUniGr<- firstUniGr[which(seqnames(firstUniGr) %in% firstIn),]
lastUniGr<- lastUniGr[which(seqnames(lastUniGr) %in% lastIn),]

firstUniMap<- findOverlaps(firstUniGr, firstGr, type="equal", select="first")
lastUniMap<- findOverlaps(lastUniGr, lastGr, type="equal", select="first")

firstUniObj<- firstObj[firstUniMap,]
lastUniObj<- lastObj[lastUniMap,]



firstUniPlus<- which(as.character(strand(firstUniGr))=="+")
firstUniMinus<- which(as.character(strand(firstUniGr))=="-")
lastUniPlus<- which(as.character(strand(lastUniGr))=="+")
lastUniMinus<- which(as.character(strand(lastUniGr))=="-")

firstUniPlusClass<- findOverlaps(firstUniGr[firstUniPlus],firstUniGr[firstUniPlus], type="start")
firstUniMinusClass<- findOverlaps(firstUniGr[firstUniMinus],firstUniGr[firstUniMinus], type="end")

lastUniPlusClass<- findOverlaps(lastUniGr[lastUniPlus],lastUniGr[lastUniPlus], type="end")
lastUniMinusClass<- findOverlaps(lastUniGr[lastUniMinus],lastUniGr[lastUniMinus], type="start")

firstUniClass<- cbind(c(firstUniPlus[queryHits(firstUniPlusClass)],firstUniMinus[queryHits(firstUniMinusClass)]), 
                      c(firstUniPlus[subjectHits(firstUniPlusClass)],firstUniMinus[subjectHits(firstUniMinusClass)]))
lastUniClass<- cbind(c(lastUniPlus[queryHits(lastUniPlusClass)],lastUniMinus[queryHits(lastUniMinusClass)]), 
                      c(lastUniPlus[subjectHits(lastUniPlusClass)],lastUniMinus[subjectHits(lastUniMinusClass)]))

########!!!!!!!!!!!!!!!
#### FIRST EXONS FOR ALL TRANSCRIPTS

classesFirst<- tapply(firstUniClass[,2],firstUniClass[,1],c)
classeslast<- tapply(lastUniClass[,2], lastUniClass[,1],c)


firstClassCnt<- aggregate(counts(firstUniObj)[firstUniClass[,2],,drop=F],list(firstUniClass[,1]),sum)
lastClassCnt<- aggregate(counts(lastUniObj)[lastUniClass[,2],,drop=F],list(lastUniClass[,1]),sum)



firstUniGrStr<- firstUniGr[firstUniClass[,2],]
firstUniGrStr<- paste(seqnames(firstUniGrStr),":",start(firstUniGrStr),"-", end(firstUniGrStr), ":", strand(firstUniGrStr), sep="")
firstClassGr<- tapply(firstUniGrStr,firstUniClass[,1], paste, collapse="/")
firstClassSeq<- gsub("\\:.*","", firstClassGr)

as.character(firstClassGr[grep("ENSG00000154358\\.",firstClassSeq)])

firstClassInd<- as.numeric(tapply(1:length(firstClassGr), firstClassGr, head, 1))
firstClassCnt<- firstClassCnt[firstClassInd,]
firstClassSeq<- firstClassSeq[firstClassInd]
firstClassGr<- firstClassGr[firstClassInd]

firstClassSum<- aggregate(firstClassCnt[,-1],list(firstClassSeq),sum)

firstClassCntNam<- firstClassCnt[,1]
firstClassCnt<- firstClassCnt[,-1]
rownames(firstClassCnt)<- firstClassGr

firstClassSumNam<- firstClassSum[,1]
firstClassSum<- firstClassSum[,-1]
rownames(firstClassSum)<- firstClassSumNam

firstClassSum<- firstClassSum[gsub("\\:.*","", rownames(firstClassCnt)),]

firstClassCnt[grep("ENSG00000154358\\..*",firstClassSeq),]

firstClassPsi<- 100*firstClassCnt/(firstClassSum+1)

table(rownames(psidat)==colnames(firstClassCnt))
# TRUE 
# 75

tempPlot<- aggregate(t(firstClassPsi[grep("ENSG00000154358\\.",rownames(firstClassPsi)),]), list(groupSam), c)

listPlot<- unlist(c(tempPlot[,-1]), recursive=FALSE)
samTypes<- tempPlot[,1]

listPlotLab<- unlist(lapply(names(listPlot), function(x){
  stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
  return(paste("chr1:", as.numeric(unlist(strsplit(x,split=":|-"))[2])-1, " ",stype, sep=""))
}))

listPlotLabType<- unlist(lapply(names(listPlot), function(x){
  stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
  return(stype)
}))

listPlotLabCoord<- unlist(lapply(names(listPlot), function(x){
  return(paste("chr1:",as.numeric(unlist(strsplit(x,split=":|-"))[2])-1, sep=""))
}))


mapFirstGr<-findOverlaps(GRanges(seqnames=sapply(strsplit(unique(listPlotLabCoord), split=":"), head, 1),
                     IRanges(start=as.numeric(sapply(strsplit(unique(listPlotLabCoord), split=":"), tail, 1)),
                             end=as.numeric(sapply(strsplit(unique(listPlotLabCoord), split=":"), tail, 1))),
             strand="+"),
             obscnUniExGr, type='end')

firstExInd<- subjectHits(mapFirstGr)


labFirstCoords<-tapply(subjectHits(mapFirstGr),as.character(unique(listPlotLabCoord)[queryHits(mapFirstGr)]),function(x){return(paste(as.character(obscnUniExGr)[x], collapse="\n"))})

firstClassNorm<- firstClassSum-firstClassCnt
library(DESeq2)
annoDf<- data.frame(group=c(groupSam,groupSam),
                    LA=c(rep("l",length(groupSam)), rep("a",length(groupSam))),
                    BATCH=c(batchVec,batchVec))
difCnt<- cbind(firstClassCnt,firstClassNorm)

dds<- DESeq2::DESeqDataSetFromMatrix(countData = difCnt, 
                                     colData = annoDf, design = ~BATCH+group+group:LA)

dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=20))

DESeq2::resultsNames(dds)

ddsDiff_first_AM_vs_AH<- DESeq2::results(dds, 
                                         contrast=list("groupADULT.MUSCLE.LAl","groupADULT.HEART.LAl"), 
                                         pAdjustMethod="BH", parallel = TRUE, 
                                         BPPARAM = SnowParam(workers=20))


ddsDiff_first_EM_vs_EH<- DESeq2::results(dds, 
                                         contrast=list("groupFETAL.MUSCLE.LAl", "groupFETAL.HEART.LAl"), 
                                         pAdjustMethod="BH", parallel = TRUE, 
                                         BPPARAM = SnowParam(workers=20))


save(ddsDiff_first_AM_vs_AH, file="ddsDiff_first_AM_vs_AH.rda")
save(ddsDiff_first_EM_vs_EH, file="ddsDiff_first_EM_vs_EH.rda")

(ddsDiff_first_AM_VS_AH_Sel<- ddsDiff_first_AM_VS_AH[grep("ENSG00000154358\\.",rownames(ddsDiff_first_AM_VS_AH)),])
(ddsDiff_first_EM_VS_EH_Sel<- ddsDiff_first_EM_VS_EH[grep("ENSG00000154358\\.",rownames(ddsDiff_first_EM_VS_EH)),])
ddsDiff_first_EM_VS_EH_Sel_End<-sapply(strsplit(rownames(ddsDiff_first_EM_VS_EH_Sel), split="\\:|\\-"), function(x){return(as.numeric(x[2])-1)})
ddsDiff_first_EM_VS_EH_Sel_Gr<- GRanges(seqnames="chr1",
                                        IRanges(start=ddsDiff_first_EM_VS_EH_Sel_End, 
                                                end=ddsDiff_first_EM_VS_EH_Sel_End),
                                        strand="+")


listPlotLabCoordEnd<- sapply(strsplit(listPlotLabCoord[seq(from=1, by=4, length.out=length(listPlotLabCoord)/4)], split="\\:|\\-"), function(x){return(as.numeric(x[2]))})
listPlotLabCoordSeqn<- sapply(strsplit(listPlotLabCoord[seq(from=1, by=4, length.out=length(listPlotLabCoord)/4)], split="\\:|\\-"), function(x){return(as.character(x[1]))})
listPlotLabCoordGr<-GRanges(seqnames=listPlotLabCoordSeqn,
                            IRanges(start=listPlotLabCoordEnd, 
                                    end=listPlotLabCoordEnd),
                            strand="+")



ddsDiff_first_EM_VS_EH_Sel_Ind<- subjectHits(findOverlaps(listPlotLabCoordGr, 
                                                          ddsDiff_first_EM_VS_EH_Sel_Gr, 
                                                          type="equal"))


#### LAST EXONS FOR ALL TRANSCRIPTS
lastUniGrStr<- lastUniGr[lastUniClass[,2],]
lastUniGrStr<- paste(seqnames(lastUniGrStr),":",start(lastUniGrStr),"-", end(lastUniGrStr), ":", strand(lastUniGrStr), sep="")
lastClassGr<- tapply(lastUniGrStr,lastUniClass[,1], paste, collapse="/")
lastClassSeq<- gsub("\\:.*","", lastClassGr)




lastClassCnt[grep("ENSG00000154358\\.",lastClassSeq),]
lastClassGr[grep("ENSG00000154358\\.",lastClassSeq)]

lastClassInd<- as.numeric(tapply(1:length(lastClassGr), lastClassGr, head, 1))
lastClassCnt<- lastClassCnt[lastClassInd,]
lastClassSeq<- lastClassSeq[lastClassInd]
lastClassGr<- lastClassGr[lastClassInd]

lastClassSum<- aggregate(lastClassCnt[,-1],list(lastClassSeq),sum)

lastClassCntNam<- lastClassCnt[,1]
lastClassCnt<- lastClassCnt[,-1]
rownames(lastClassCnt)<- lastClassGr

lastClassSumNam<- lastClassSum[,1]
lastClassSum<- lastClassSum[,-1]
rownames(lastClassSum)<- lastClassSumNam

lastClassSum<- lastClassSum[gsub("\\:.*","", rownames(lastClassCnt)),]

lastClassCnt[grep("ENSG00000154358\\.",lastClassSeq),]

lastClassPsi<- 100*lastClassCnt/(lastClassSum+1)

table(rownames(psidat)==colnames(lastClassCnt))
# TRUE 
# 75

tempPlot<- aggregate(t(lastClassPsi[grep("ENSG00000154358\\.",rownames(lastClassPsi)),]), list(groupSam), c)

listPlot<- unlist(c(tempPlot[,-1]), recursive=FALSE)
samTypes<- tempPlot[,1]

listPlotLab<- unlist(lapply(names(listPlot), function(x){
  stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
  return(paste("chr1:", as.numeric(unlist(strsplit(x,split=":|-"))[3])+1, " ",stype, sep=""))
}))

lastClassNorm<- lastClassSum-lastClassCnt
library(DESeq2)
annoDf<- data.frame(group=c(groupSam,groupSam),
                    LA=c(rep("l",length(groupSam)), rep("a",length(groupSam))),
                    BATCH=c(batchVec,batchVec))
difCnt<- cbind(lastClassCnt,lastClassNorm)
indChooseCol= which(annoDf$group%in%c("ADULT MUSCLE", "ADULT HEART"))
dds<- DESeq2::DESeqDataSetFromMatrix(countData = difCnt[,indChooseCol], 
                                     colData = annoDf[indChooseCol,], design = ~BATCH+group+group:LA)

dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=20))

DESeq2::resultsNames(dds)


ddsDiff_last_AM_vs_AH<- DESeq2::results(dds, 
                                        contrast=list("groupADULT.MUSCLE.LAl", "groupADULT.HEART.LAl"), 
                                        pAdjustMethod="BH", parallel = TRUE, 
                                        BPPARAM = SnowParam(workers=15))


indChooseCol= which(annoDf$group%in%c("FETAL MUSCLE", "FETAL HEART"))
dds<- DESeq2::DESeqDataSetFromMatrix(countData = difCnt[,indChooseCol], 
                                     colData = annoDf[indChooseCol,], design = ~BATCH+group+group:LA)

dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=20))

DESeq2::resultsNames(dds)

ddsDiff_last_EM_vs_EH<- DESeq2::results(dds, 
                                        contrast=list("groupFETAL.MUSCLE.LAl", "groupFETAL.HEART.LAl"), 
                                        pAdjustMethod="BH", parallel = TRUE, 
                                        BPPARAM = SnowParam(workers=15))

save(ddsDiff_last_AM_vs_AH,file="ddsDiff_last_AM_vs_AH.rda")
save(ddsDiff_last_EM_vs_EH,file="ddsDiff_last_EM_vs_EH.rda")


indChooseCol= which(annoDf$group%in%c("ADULT MUSCLE", "FETAL MUSCLE"))
dds<- DESeq2::DESeqDataSetFromMatrix(countData = difCnt[,indChooseCol], 
                                     colData = annoDf[indChooseCol,], design = ~BATCH+group+group:LA)

dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=20))

ddsDiff_last_AM_vs_EM<- DESeq2::results(dds, 
                                        contrast=list("groupADULT.MUSCLE.LAl", "groupFETAL.MUSCLE.LAl"), 
                                        pAdjustMethod="BH", parallel = TRUE, 
                                        BPPARAM = SnowParam(workers=15))

indChooseCol= which(annoDf$group%in%c("ADULT HEART", "FETAL HEART"))
dds<- DESeq2::DESeqDataSetFromMatrix(countData = difCnt[,indChooseCol], 
                                     colData = annoDf[indChooseCol,], design = ~group+group:LA)

dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=20))


ddsDiff_last_AH_vs_EH<- DESeq2::results(dds, 
                                        contrast=list("groupADULT.HEART.LAl", "groupFETAL.HEART.LAl"), 
                                        pAdjustMethod="BH", parallel = TRUE, 
                                        BPPARAM = SnowParam(workers=15))


save(ddsDiff_last_AM_vs_EM,file="ddsDiff_last_AM_vs_EM.rda")
save(ddsDiff_last_AH_vs_EH,file="ddsDiff_last_AH_vs_EH.rda")




(ddsDiff_last_AM_vs_AH_Sel<- ddsDiff_last_AM_vs_AH[grep("ENSG00000154358\\.",rownames(ddsDiff_last_AM_vs_AH)),])

(ddsDiff_last_EM_vs_EH_Sel<- ddsDiff_last_EM_vs_EH[grep("ENSG00000154358\\.",rownames(ddsDiff_last_EM_vs_EH)),])

(ddsDiff_last_AM_vs_EM_Sel<- ddsDiff_last_AM_vs_EM[grep("ENSG00000154358\\.",rownames(ddsDiff_last_AM_vs_EM)),])

(ddsDiff_last_AH_vs_EH_Sel<- ddsDiff_last_AH_vs_EH[grep("ENSG00000154358\\.",rownames(ddsDiff_last_AH_vs_EH)),])

ddsDiff_last_EM_vs_EH_Sel_Beg<-sapply(strsplit(rownames(ddsDiff_last_EM_vs_EH_Sel), split="\\:|\\-"), function(x){return(as.numeric(x[3])+1)})

ddsDiff_last_EM_vs_EH_Sel_Gr<- GRanges(seqnames="chr1",
                                       IRanges(start=ddsDiff_last_EM_vs_EH_Sel_Beg, 
                                               end=ddsDiff_last_EM_vs_EH_Sel_Beg),
                                       strand="+")

ddsDiff_last_EM_vs_EH_Beg<-sapply(strsplit(rownames(ddsDiff_last_EM_vs_EH), split="\\:|\\-"), function(x){return(as.numeric(x[3])+1)})
ddsDiff_last_EM_vs_EH_Gr<- GRanges(seqnames="chr1",
                                   IRanges(start=ddsDiff_last_EM_vs_EH_Beg, 
                                           end=ddsDiff_last_EM_vs_EH_Beg),
                                   strand="+")

listPlotLab<- unlist(lapply(names(listPlot), function(x){
  stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
  return(paste("chr1:", as.numeric(unlist(strsplit(x,split=":|-"))[3])+1, " ",stype, sep=""))
}))

listPlotLabType<- unlist(lapply(names(listPlot), function(x){
  stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
  return(stype)
}))

listPlotLabCoord<- unlist(lapply(names(listPlot), function(x){
  return(paste("chr1:",as.numeric(unlist(strsplit(x,split=":|-"))[3])+1, sep=""))
}))

mapLastGr<-findOverlaps(GRanges(seqnames=sapply(strsplit(unique(listPlotLabCoord), split=":"), head, 1),
                                IRanges(start=as.numeric(sapply(strsplit(unique(listPlotLabCoord), split=":"), tail, 1)),
                                        end=as.numeric(sapply(strsplit(unique(listPlotLabCoord), split=":"), tail, 1))),
                                strand="+"),
                        obscnUniExGr, type='start', select="all")

labLastCoords<-tapply(subjectHits(mapLastGr),as.character(unique(listPlotLabCoord)[queryHits(mapLastGr)]),function(x){return(paste(as.character(obscnUniExGr)[x], collapse="\n"))})

names(labLastCoords)<- unique(listPlotLabCoord)
listPlotLabCoordEnd<- sapply(strsplit(listPlotLabCoord[seq(from=1, by=4, length.out=length(listPlotLabCoord)/4)], split="\\:|\\-"), function(x){return(as.numeric(x[2]))})
listPlotLabCoordSeqn<- sapply(strsplit(listPlotLabCoord[seq(from=1, by=4, length.out=length(listPlotLabCoord)/4)], split="\\:|\\-"), function(x){return(as.character(x[1]))})
listPlotLabCoordGr<-GRanges(seqnames=listPlotLabCoordSeqn,
                            IRanges(start=listPlotLabCoordEnd, 
                                    end=listPlotLabCoordEnd),
                            strand="+")


ddsDiff_last_EM_vs_EH_Sel_Ind<- subjectHits(findOverlaps(listPlotLabCoordGr, 
                                                         ddsDiff_last_EM_vs_EH_Sel_Gr, 
                                                         type="equal"))

drawBracket<- function(x=c(), y=c(), l=2, col=1){
  lines(x=c(x[1],x[2]),y=c(y,y),col="black")
  lines(x=c(x[2],x[2]),y=c(y,y+l),col="black")
  lines(x=c(x[1],x[1]),y=c(y,y+l),col="black")
}

tmpOrd<-4:1
listPlot<- listPlot[c(tmpOrd, (4+tmpOrd), (8+tmpOrd))]
listPlotLabType<-listPlotLabType[c(tmpOrd, (4+tmpOrd), (8+tmpOrd))]






pdf("fig5.pdf", width=20, height=14, pointsize=14)
par(mar=c(25.1, 4.1, 4.1, 2.1))
ordIdx<-c(3,1,4,2)
tmpOrdIdx<-c(2,4,1,3)
cols<- rainbow(length(samTypes))
cols[4]<- "#FF80FF"
listPlotTmp<- listPlot


resPlot<-c()

resPlot$stats<- lapply(listPlot, function(x){return(as.matrix(summary(x)))})
resPlot$stats<- do.call(cbind.data.frame, resPlot$stats)[-1,]

lowSamInd<- as.numeric(which(unlist(lapply(listPlotTmp, length)<5)))
listPlotTmp[lowSamInd]<- NA

resPlot2<- boxplot(listPlotTmp, main="", 
                  names=listPlotLabType, border=rep(cols[ordIdx],3), las=2, ylim=c(0,100),
                  ylab=bquote(Psi ~ .("(%)")))

for(k in lowSamInd)
  points(x=rep(k,length(listPlot[[k]])), y=listPlot[[k]], col=rep(cols[ordIdx],3)[k], pch=1)


for(k in 1:length(listPlot))
  points(x=k, y=mean(listPlot[[k]]), col= rep(cols[ordIdx],3)[k], pch=18, cex=1.5)


lapply(seq(from=4.5, by=4, length.out=length(listPlot)/4-1), function(x){abline(v=x, col="grey")})


mains<- labLastCoords[listPlotLabCoord][seq(from=1, by=4, length.out=length(listPlotLabCoord)/4)]

# Add exon numbers from mapLastGr I generated earlier
namMains<- names(mains)
mains[1:2]<- paste("Exon ", subjectHits(mapLastGr)[1:2], ": ", mains[1:2], sep="")
mains[3]<- paste("Exon ", subjectHits(mapLastGr)[3], ": ", 
                 gsub("\n",paste("\nExon ", subjectHits(mapLastGr)[4], ": ", sep=""), mains[3]), sep="")
names(mains)<- namMains

for(i in 1:length(mains))
  mtext(text=mains[i], side = 3, line = 1, outer = FALSE, at = (i-1)*4+2.5, font=2)



FDRAMAH<- ddsDiff_last_AM_vs_AH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"padj"]
PAMAH<- ddsDiff_last_AM_vs_AH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"pvalue"]

FDREMEH<- ddsDiff_last_EM_vs_EH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"padj"]
PEMEH<- ddsDiff_last_EM_vs_EH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"pvalue"]

FDRAMEM<- ddsDiff_last_AM_vs_EM_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"padj"]
PAMEM<- ddsDiff_last_AM_vs_EM_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"pvalue"]

FDRAHEH<- ddsDiff_last_AH_vs_EH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"padj"]
PAHEH<- ddsDiff_last_AH_vs_EH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"pvalue"]

delPsi<- as.numeric(resPlot$stats[3,seq(from=1,by=2,length.out=length(listPlot)/2)]-resPlot$stats[3,seq(from=2,by=2,length.out=length(listPlot)/2)])
labs<- list()
for(i in 1:(length(delPsi))){
  if(i %% 2 == 1){
    labs<- c(labs, list(FDRAHEH[trunc(i/2)+1],PAHEH[trunc(i/2)+1],delPsi[i]) )
  } else {
    labs<- c(labs, list(FDRAMEM[trunc(i/2)],PAMEM[trunc(i/2)],delPsi[i]) )
  }
}
labs<- unlist(labs)

tmpNumb<- format(round(labs, 3), nsmall=3, digits=3)
tmpNumb[which(round(labs, 4)==0)]<- format(labs[which(round(labs, 4)==0)], nsmall=4, digits=4, scientific=T)
labelsPlot<- 
  paste(rep(c("FDR(AH/FH)","P(AH/FH)", "DeltaPsi(AH/FH)", 
              "FDR(AM/FM)","P(AM/FM)", "DeltaPsi(AM/FM)"),(length(labs)/6)),
        gsub(" ","",tmpNumb), sep="=")
labPlot=c()
for(i in 1:(length(labelsPlot)/3)){
  labPlot<- c(labPlot,paste(labelsPlot[(i-1)*3+(1:3)], collapse="\n"))
}


sig<- F
for(i in 1:length(labelsPlot)){
  if((i%%3) == 1)
    sig<- (as.numeric(unlist(strsplit(labelsPlot[i], split="="))[2]) < 0.05 &
             abs(as.numeric(unlist(strsplit(labelsPlot[i+2], split="="))[2])) > 10)
  if((i%%3) != 0){
    if(!sig){
      mtext(text=labelsPlot[i], side = 1, line = 8+((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=1.5, by=2, length.out=(2*length(mains)))[trunc((i-1)/3)+1], col="black")
    } else {
      mtext(text=paste("*",labelsPlot[i],sep=""), side = 1, line = 8+((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=1.5, by=2, length.out=(2*length(mains)))[trunc((i-1)/3)+1], col="red")
    }
  } else {
    if(!sig){
      mtext(text=bquote(Delta*Psi ~ .(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 8+ ((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=1.5, by=2, length.out=(2*length(mains)))[trunc((i-1)/3)+1], col="black")
    } else{
      mtext(text=bquote(.("*")*Delta*Psi*.(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 8+ ((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=1.5, by=2, length.out=(2*length(mains)))[trunc((i-1)/3)+1], col="red")
    }
  }
}

par(xpd=TRUE)
for(i in seq(from=1, by=2, length.out=6))
  drawBracket(x=c(i,i+1), y=-34)





# Plot second set of x axis info P and FDR and etc
drawBracket(x=c(1,3), y=-55)
drawBracket(x=c(5,7), y=-55)
drawBracket(x=c(9,11), y=-55)

delPsi<- resPlot$stats[3,seq(from=1,by=4,length.out=length(listPlot)/4)] - resPlot$stats[3,seq(from=3,by=4,length.out=length(listPlot)/4)] 
labs<- list()
for(i in 1:(length(delPsi))){
  labs<- c(labs, list(FDRAMAH[i],PAMAH[i],delPsi[i]) )
}
labs<- unlist(labs)
tmpNumb<- format(round(labs, 3), nsmall=3, digits=3)
tmpNumb[which(round(labs, 4)==0)]<- format(labs[which(round(labs, 4)==0)], nsmall=4, digits=4, scientific=T)

labelsPlot<- 
  paste(rep(c("FDR(AM/AH)","P(AM/AH)", "DeltaPsi(AM/AH)"),(length(labs)/3)),
        gsub(" ","",tmpNumb), sep="=")

labPlot=c()
for(i in 1:(length(labelsPlot)/3)){
  labPlot<- c(labPlot,paste(labelsPlot[(i-1)*3+(1:3)], collapse="\n"))
}


sig<- F
for(i in 1:length(labelsPlot)){
  if((i%%3) == 1)
    sig<- (as.numeric(unlist(strsplit(labelsPlot[i], split="="))[2]) < 0.05 &
             abs(as.numeric(unlist(strsplit(labelsPlot[i+2], split="="))[2])) > 10)
  if((i%%3) != 0){
    if(!sig){
      mtext(text=labelsPlot[i], side = 1, line = 14+((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=2, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="black")
    } else {
      mtext(text=paste("*",labelsPlot[i],sep=""), side = 1, line = 14+((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=2, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="red")
    }
  } else {
    if(!sig){
      mtext(text=bquote(Delta*Psi ~ .(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 14+ ((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=2, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="black")
    } else{
      mtext(text=bquote(.("*")*Delta*Psi*.(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 14+ ((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=2, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="red")
    }
  }
}


# Plot third set of x axis info P and FDR and etc
drawBracket(x=c(2,4), y=-76)
drawBracket(x=c(6,8), y=-76)
drawBracket(x=c(10,12), y=-76)

FDRAHEH<- ddsDiff_last_AH_vs_EH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"padj"]
PAHEH<- ddsDiff_last_AH_vs_EH_Sel[ddsDiff_last_EM_vs_EH_Sel_Ind,"pvalue"]

delPsi<- resPlot$stats[3,seq(from=2,by=4,length.out=length(listPlot)/4)] - resPlot$stats[3,seq(from=4,by=4,length.out=length(listPlot)/4)]
labs<- list()
for(i in 1:(length(delPsi))){
  labs<- c(labs, list(FDREMEH[i],PEMEH[i],delPsi[i]) )
}
labs<- unlist(labs)
tmpNumb<- format(round(labs, 3), nsmall=3, digits=3)
tmpNumb[which(round(labs, 4)==0)]<- format(labs[which(round(labs, 4)==0)], nsmall=4, digits=4, scientific=T)

labelsPlot<- 
  paste(rep(c("FDR(FM/FH)","P(FM/FH)", "DeltaPsi(FM/FH)"),(length(labs)/3)),
        gsub(" ","",tmpNumb), sep="=")

labPlot=c()
for(i in 1:(length(labelsPlot)/3)){
  labPlot<- c(labPlot,paste(labelsPlot[(i-1)*3+(1:3)], collapse="\n"))
}


sig<- F
for(i in 1:length(labelsPlot)){
  if((i%%3) == 1)
    sig<- (as.numeric(unlist(strsplit(labelsPlot[i], split="="))[2]) < 0.05 &
             abs(as.numeric(unlist(strsplit(labelsPlot[i+2], split="="))[2])) > 10)
  if((i%%3) != 0){
    if(!sig){
      mtext(text=labelsPlot[i], side = 1, line = 20+((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=3, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="black")
    } else {
      mtext(text=paste("*",labelsPlot[i],sep=""), side = 1, line = 20+((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=3, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="red")
    }
  } else {
    if(!sig){
      mtext(text=bquote(Delta*Psi ~ .(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 20+ ((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=3, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="black")
    } else{
      mtext(text=bquote(.("*")*Delta*Psi*.(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 20+ ((i-1)%%3)+1, outer = FALSE, 
            at = seq(from=3, by=4, length.out=(1*length(mains)))[trunc((i-1)/3)+1], col="red")
    }
  }
}

dev.off()



##### Big Boxplot 

pdf("Incusion_levels_of_all_exons.pdf", width=18, height=138, pointsize=14)
ordIdx<-c(3,1,4,2)
colInd<- c(2,4,1,3)
cols<- rainbow(length(samTypes))
cols[4]<- "#FF80FF"
  par(mar=c(31.1, 5.1, 5.1, 2.1))
  par(cex.lab=1.7)
  par(cex.axis=1.7)
  par(cex.main=1.7)
  xLim<- c(0,Inf)
  yLim<- c(1,100)
  if(yLim[1]==1){
    limPlot<- 125
    if(length(xLim[1])>0){
      strTmp<- "+"
      ordInd<- order(start(obscnUniExGr), decreasing=FALSE)
      selInd<- which(start(obscnUniExGr)>=xLim[1] & (start(obscnUniExGr)<=xLim[2]))
      exonNo<- selInd
      if(length(selInd)>0){
        if( strTmp== "-")
          exonNo<- length(start(obscnUniExGr))-selInd+1
        if(length(selInd)>limPlot)
          selInd<- selInd[1:limPlot]
        psidat[,selInd]
        exonNo<- exonNo[1:limPlot]
        
        par(mfrow=c(21,6))
        
        for (i in 1:length(selInd)){
          plotList<- tapply(100*psidat[,selInd[i]], groupSam, c)
          grp<- rep("A",length(groupSam))
          grp[which(groupSam!=levels(groupSam)[1] & groupSam!=levels(groupSam)[2] )]<- "O"
          dattmp<- data.frame(exinc=psidat[,selInd[i]], group=grp)
          
          grp2<- rep("A",length(groupSam))
          grp2[grep("HEART",groupSam)]<- "O"
          dattmp2<- data.frame(exinc=psidat[,selInd[i]], group=grp2)
          
 
          pvalTxt1<- paste("FDR(M) = ", format(difMuscleSel[selInd[i],"padj"], digits =3), sep="")
          pvalTxt2<- paste("P(M) = ", format(difMuscleSel[selInd[i],"pvalue"], digits =3), sep="")
          pvalTxt3<- paste("(M) = ", format(100*(mean(psimean[selInd[i], c("ADULT MUSCLE", "FETAL MUSCLE")])-mean(psimean[selInd[i], c("ADULT HEART", "FETAL HEART")])), digits =3), sep="")
          
          pvalTxt4<- paste("FDR(A) = ", format(difAdultSel[selInd[i],"padj"], digits =3), sep="")
          pvalTxt5<- paste("P(A) = ", format(difAdultSel[selInd[i],"pvalue"], digits =3), sep="")
          pvalTxt6<- paste("(A) = ", format(100*(mean(psimean[selInd[i], c("ADULT MUSCLE", "ADULT HEART")])-mean(psimean[selInd[i], c("FETAL MUSCLE", "FETAL HEART")])), digits =3), sep="")
          
          pvalTxt7<- paste("FDR(",names(difListSel)[1],") = ", format(difListSel[[1]][selInd[i],"padj"], digits =3), sep="")
          pvalTxt8<- paste("P(",names(difListSel)[1],") = ", format(difListSel[[1]][selInd[i],"pvalue"], digits =3), sep="")
          pvalTxt9<- paste("(",names(difListSel)[1],") = ", format(100*(psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("FETAL MUSCLE")]), digits =3), sep="")
          
          pvalTxt10<- paste("FDR(",names(difListSel)[2],") = ", format(difListSel[[2]][selInd[i],"padj"], digits =3), sep="")
          pvalTxt11<- paste("P(",names(difListSel)[2],") = ", format(difListSel[[2]][selInd[i],"pvalue"], digits =3), sep="")
          pvalTxt12<- paste("(",names(difListSel)[2],") = ", format(100*(psimean[selInd[i], c("ADULT HEART")]-psimean[selInd[i], c("FETAL HEART")]), digits =3), sep="")
          
          pvalTxt13<- paste("FDR(",names(difListSel)[3],") = ", format(difListSel[[3]][selInd[i],"padj"], digits =3), sep="")
          pvalTxt14<- paste("P(",names(difListSel)[3],") = ", format(difListSel[[3]][selInd[i],"pvalue"], digits =3), sep="")
          pvalTxt15<- paste("(",names(difListSel)[3],") = ", format(100*(psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("ADULT HEART")]), digits =3), sep="")
          
          pvalTxt16<- paste("FDR(",names(difListSel)[4],") = ", format(difListSel[[4]][selInd[i],"padj"], digits =3), sep="")
          pvalTxt17<- paste("P(",names(difListSel)[4],") = ", format(difListSel[[4]][selInd[i],"pvalue"], digits =3), sep="")
          pvalTxt18<- paste("(",names(difListSel)[4],") = ", format(100*(psimean[selInd[i], c("FETAL MUSCLE")]-psimean[selInd[i], c("FETAL HEART")]), digits =3), sep="")
          
          
          pvalTxt1[difMuscleSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt1[difMuscleSel[selInd[i],"padj"]<0.05], sep="")
          pvalTxt2[difMuscleSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt2[difMuscleSel[selInd[i],"padj"]<0.05], sep="")
          #pvalTxt3[difMuscleSel[selInd[i],"padj"]<0.05]<- bquote(.("*") ~ pvalTxt3[difMuscleSel[selInd[i],"padj"]<0.05])
          
          pvalTxt4[difAdultSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt4[difAdultSel[selInd[i],"padj"]<0.05], sep="")
          pvalTxt5[difAdultSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt5[difAdultSel[selInd[i],"padj"]<0.05], sep="")
          #pvalTxt6[difAdultSel[selInd[i],"padj"]<0.05]<- bquote(.("*") ~  pvalTxt6[difAdultSel[selInd[i],"padj"]<0.05])
          
          pvalTxt7[difListSel[[1]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt7[difListSel[[1]][selInd[i],"padj"]<0.05], sep="")
          pvalTxt8[difListSel[[1]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt8[difListSel[[1]][selInd[i],"padj"]<0.05], sep="")
          #pvalTxt9[difListSel[[1]][selInd[i],"padj"]<0.05]<- bquote(.("*") ~  pvalTxt9[difListSel[[1]][selInd[i],"padj"]<0.05])
          
          pvalTxt10[difListSel[[2]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt10[difListSel[[2]][selInd[i],"padj"]<0.05], sep="")
          pvalTxt11[difListSel[[2]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt11[difListSel[[2]][selInd[i],"padj"]<0.05], sep="")
          #pvalTxt12[difListSel[[2]][selInd[i],"padj"]<0.05]<- bquote(.("*") ~  pvalTxt12[difListSel[[2]][selInd[i],"padj"]<0.05])
          
          pvalTxt13[difListSel[[3]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt13[difListSel[[3]][selInd[i],"padj"]<0.05], sep="")
          pvalTxt14[difListSel[[3]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt14[difListSel[[3]][selInd[i],"padj"]<0.05], sep="")
          #pvalTxt15[difListSel[[3]][selInd[i],"padj"]<0.05]<- bquote(.("*") ~  pvalTxt15[difListSel[[3]][selInd[i],"padj"]<0.05])
          
          pvalTxt16[difListSel[[4]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt16[difListSel[[4]][selInd[i],"padj"]<0.05], sep="")
          pvalTxt17[difListSel[[4]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt17[difListSel[[4]][selInd[i],"padj"]<0.05], sep="")
          #pvalTxt18[difListSel[[4]][selInd[i],"padj"]<0.05]<- bquote(.("*") ~  pvalTxt18[difListSel[[4]][selInd[i],"padj"]<0.05])
          
          pvals<- c(difMuscleSel[selInd[i],"padj"],
                    difAdultSel[selInd[i],"padj"],
                    difListSel[[1]][selInd[i],"padj"],
                    difListSel[[2]][selInd[i],"padj"],
                    difListSel[[3]][selInd[i],"padj"],
                    difListSel[[4]][selInd[i],"padj"])
          
          delpsis<- 100*c(mean(psimean[selInd[i], c("ADULT MUSCLE", "FETAL MUSCLE")])-mean(psimean[selInd[i], c("ADULT HEART", "FETAL HEART")]),
                   mean(psimean[selInd[i], c("ADULT MUSCLE", "ADULT HEART")])-mean(psimean[selInd[i], c("FETAL MUSCLE", "FETAL HEART")]),
                   psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("FETAL MUSCLE")],
                   psimean[selInd[i], c("ADULT HEART")]-psimean[selInd[i], c("FETAL HEART")],
                   psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("ADULT HEART")],
                   psimean[selInd[i], c("FETAL MUSCLE")]-psimean[selInd[i], c("FETAL HEART")])
          

          
          titleList<- list(pvalTxt1,
                           pvalTxt2,
                           pvalTxt3,
                           pvalTxt4,
                           pvalTxt5,
                           pvalTxt6,
                           pvalTxt7,
                           pvalTxt8,
                           pvalTxt9,
                           pvalTxt10,
                           pvalTxt11,
                           pvalTxt12,
                           pvalTxt13,
                           pvalTxt14,
                           pvalTxt15,
                           pvalTxt16,
                           pvalTxt17,
                           pvalTxt18)
          
          colTitles<- rep("black", length(titleList))
          colTitles[c(((which(pvals<0.05 & abs(delpsis)>10))*3),((which(pvals<0.05 & abs(delpsis)>10))*3-1),((which(pvals<0.05 & abs(delpsis)>10))*3-2))]<- "red"
          
 
          plotList2<- plotList
          lowSamInd<- which(unlist(lapply(plotList,length))<5)
          plotList2[lowSamInd]<- NA
          
          boxplot(plotList2[c("ADULT MUSCLE","ADULT HEART","FETAL MUSCLE","FETAL HEART")[ordIdx]], main=paste(paste("Exon", exonNo[selInd[i]]), 
                                                                                                             paste(paste(seqnames(obscnUniExGr)[selInd[i]], start(obscnUniExGr)[selInd[i]], sep=":"),
                                                                                                                   end(obscnUniExGr)[selInd[i]], sep="-\n"),sep="\n"), 
                  ylim=c(0,100),
                  ylab=bquote(Psi*.(" (%)")),
                  xlab="",border=cols[colInd], 
                  names=c("AM","AH","FM","FH")[ordIdx], las=2)
          
          for(k in lowSamInd)
            points(x=rep(k,length(plotList[[k]])), y=plotList[[k]], col=cols[colInd][k], pch=1)
          
          
          for(k in 1:length(plotList))
            points(x=k, y=mean(plotList[[k]]), col= cols[colInd][k], pch=18, cex=1.5)
          
          #title(xlab =pvalTxt, line = 20)
          for(l in 1:length(titleList)){
            if(l %% 3 == 0){
              if(length(grep("^\\*",titleList[[l-1]]))==0){
                title(xlab=bquote(Delta*Psi*.(titleList[[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
              } else {
                title(xlab=bquote(.("*")*Delta*Psi*.(titleList[[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
              }
            } else {
              title(xlab=titleList[[l]], line = 2+(l*1.5), col.lab = colTitles[l])
            }
          }
          
          
        }
      }
    }
  }

dev.off()

save.image(file="FigS1Fig3.rda")

#####
pdf("fig4.pdf", width=20, height=30, pointsize=16)

ordIdx<-c(3,1,4,2)
colInd<- c(2,4,1,3)

cols<- rainbow(length(samTypes))
cols[4]<- "#FF80FF"

par(mar=c(31.1, 5.1, 6.1, 2.1))
par(cex.lab=1.5)
par(cex.axis=1.5)
par(cex.main=1.5)

xLim<- c(0,Inf)
yLim<- c(1,100)
if(yLim[1]==1){
  limPlot<- 125
  if(length(xLim[1])>0){
    strTmp<- "+"
    ordInd<- order(start(obscnUniExGr), decreasing=FALSE)
    selInd<- which(start(obscnUniExGr)>=xLim[1] & (start(obscnUniExGr)<=xLim[2]))
    exonNo<- selInd
    if(length(selInd)>0){
      if( strTmp== "-")
        exonNo<- length(start(obscnUniExGr))-selInd+1
      if(length(selInd)>limPlot)
        selInd<- selInd[1:limPlot]
      psidat[,selInd]
      exonNo<- exonNo[1:limPlot]
      
      par(mfrow=c(3,5))
      
      
      pvalsList<- list(difMuscleSel[selInd,"padj"],
                       difAdultSel[selInd,"padj"],
                       difListSel[[1]][selInd,"padj"],
                       difListSel[[2]][selInd,"padj"],
                       difListSel[[3]][selInd,"padj"],
                       difListSel[[4]][selInd,"padj"])
      
      delpsisList<- list(100*(rowMeans(psimean[selInd, c("ADULT MUSCLE", "FETAL MUSCLE")])-rowMeans(psimean[selInd, c("ADULT HEART", "FETAL HEART")])),
                         100*(rowMeans(psimean[selInd, c("ADULT MUSCLE", "ADULT HEART")])-rowMeans(psimean[selInd, c("FETAL MUSCLE", "FETAL HEART")])),
                         100*(psimean[selInd, c("ADULT MUSCLE")]-psimean[selInd, c("FETAL MUSCLE")]),
                         100*(psimean[selInd, c("ADULT HEART")]-psimean[selInd, c("FETAL HEART")]),
                         100*(psimean[selInd, c("ADULT MUSCLE")]-psimean[selInd, c("ADULT HEART")]),
                         100*(psimean[selInd, c("FETAL MUSCLE")]-psimean[selInd, c("FETAL HEART")]))
      
      selInd<- selInd[which( ((pvalsList[[1]]<0.05 & abs(delpsisList[[1]])>10) +
                                (pvalsList[[2]]<0.05 & abs(delpsisList[[2]])>10)+
                                (pvalsList[[3]]<0.05 & abs(delpsisList[[3]])>10)+
                                (pvalsList[[4]]<0.05 & abs(delpsisList[[4]])>10)+
                                (pvalsList[[5]]<0.05 & abs(delpsisList[[5]])>10)+
                                (pvalsList[[6]]<0.05 & abs(delpsisList[[6]])>10)) >= 2
      )]
      subFig<- LETTERS[1:length(selInd)]      
      for (i in 1:length(selInd)){
        plotList<- tapply(100*psidat[,selInd[i]], groupSam, c)
        grp<- rep("A",length(groupSam))
        grp[which(groupSam!=levels(groupSam)[1] & groupSam!=levels(groupSam)[2] )]<- "O"
        dattmp<- data.frame(exinc=psidat[,selInd[i]], group=grp)
        
        grp2<- rep("A",length(groupSam))
        grp2[grep("HEART",groupSam)]<- "O"
        dattmp2<- data.frame(exinc=psidat[,selInd[i]], group=grp2)
        
        pvalTxt1<- paste("FDR(M) = ", format(difMuscleSel[selInd[i],"padj"], digits =3), sep="")
        pvalTxt2<- paste("P(M) = ", format(difMuscleSel[selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt3<- paste("(M) = ", format(100*(mean(psimean[selInd[i], c("ADULT MUSCLE", "FETAL MUSCLE")])-mean(psimean[selInd[i], c("ADULT HEART", "FETAL HEART")])), digits =3), sep="")
        
        pvalTxt4<- paste("FDR(A) = ", format(difAdultSel[selInd[i],"padj"], digits =3), sep="")
        pvalTxt5<- paste("P(A) = ", format(difAdultSel[selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt6<- paste("(A) = ", format(100*(mean(psimean[selInd[i], c("ADULT MUSCLE", "ADULT HEART")])-mean(psimean[selInd[i], c("FETAL MUSCLE", "FETAL HEART")])), digits =3), sep="")
        
        pvalTxt7<- paste("FDR(",names(difListSel)[1],") = ", format(difListSel[[1]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt8<- paste("P(",names(difListSel)[1],") = ", format(difListSel[[1]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt9<- paste("(",names(difListSel)[1],") = ", format(100*(psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("FETAL MUSCLE")]), digits =3), sep="")
        
        pvalTxt10<- paste("FDR(",names(difListSel)[2],") = ", format(difListSel[[2]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt11<- paste("P(",names(difListSel)[2],") = ", format(difListSel[[2]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt12<- paste("(",names(difListSel)[2],") = ", format(100*(psimean[selInd[i], c("ADULT HEART")]-psimean[selInd[i], c("FETAL HEART")]), digits =3), sep="")
        
        pvalTxt13<- paste("FDR(",names(difListSel)[3],") = ", format(difListSel[[3]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt14<- paste("P(",names(difListSel)[3],") = ", format(difListSel[[3]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt15<- paste("(",names(difListSel)[3],") = ", format(100*(psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("ADULT HEART")]), digits =3), sep="")
        
        pvalTxt16<- paste("FDR(",names(difListSel)[4],") = ", format(difListSel[[4]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt17<- paste("P(",names(difListSel)[4],") = ", format(difListSel[[4]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt18<- paste("(",names(difListSel)[4],") = ", format(100*(psimean[selInd[i], c("FETAL MUSCLE")]-psimean[selInd[i], c("FETAL HEART")]), digits =3), sep="")
        
        pvalTxt1[difMuscleSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt1[difMuscleSel[selInd[i],"padj"]<0.05], sep="")
        pvalTxt2[difMuscleSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt2[difMuscleSel[selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt4[difAdultSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt4[difAdultSel[selInd[i],"padj"]<0.05], sep="")
        pvalTxt5[difAdultSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt5[difAdultSel[selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt7[difListSel[[1]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt7[difListSel[[1]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt8[difListSel[[1]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt8[difListSel[[1]][selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt10[difListSel[[2]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt10[difListSel[[2]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt11[difListSel[[2]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt11[difListSel[[2]][selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt13[difListSel[[3]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt13[difListSel[[3]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt14[difListSel[[3]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt14[difListSel[[3]][selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt16[difListSel[[4]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt16[difListSel[[4]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt17[difListSel[[4]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt17[difListSel[[4]][selInd[i],"padj"]<0.05], sep="")
        
        
        
        pvals<- c(difMuscleSel[selInd[i],"padj"],
                  difAdultSel[selInd[i],"padj"],
                  difListSel[[1]][selInd[i],"padj"],
                  difListSel[[2]][selInd[i],"padj"],
                  difListSel[[3]][selInd[i],"padj"],
                  difListSel[[4]][selInd[i],"padj"])
        
        delpsis<- 100*c(mean(psimean[selInd[i], c("ADULT MUSCLE", "FETAL MUSCLE")])-mean(psimean[selInd[i], c("ADULT HEART", "FETAL HEART")]),
                        mean(psimean[selInd[i], c("ADULT MUSCLE", "ADULT HEART")])-mean(psimean[selInd[i], c("FETAL MUSCLE", "FETAL HEART")]),
                        psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("FETAL MUSCLE")],
                        psimean[selInd[i], c("ADULT HEART")]-psimean[selInd[i], c("FETAL HEART")],
                        psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("ADULT HEART")],
                        psimean[selInd[i], c("FETAL MUSCLE")]-psimean[selInd[i], c("FETAL HEART")])
        
        titleList<- list(pvalTxt1,
                         pvalTxt2,
                         pvalTxt3,
                         pvalTxt4,
                         pvalTxt5,
                         pvalTxt6,
                         pvalTxt7,
                         pvalTxt8,
                         pvalTxt9,
                         pvalTxt10,
                         pvalTxt11,
                         pvalTxt12,
                         pvalTxt13,
                         pvalTxt14,
                         pvalTxt15,
                         pvalTxt16,
                         pvalTxt17,
                         pvalTxt18)
        
        colTitles<- rep("black", length(titleList))
        colTitles[c(((which(pvals<0.05 & abs(delpsis)>10))*3),((which(pvals<0.05 & abs(delpsis)>10))*3-1),((which(pvals<0.05 & abs(delpsis)>10))*3-2))]<- "red"
        
        plotList2<- plotList
        lowSamInd<- which(unlist(lapply(plotList,length))<5)
        plotList2[lowSamInd]<- NA

        boxplot(plotList2[c("ADULT MUSCLE","ADULT HEART","FETAL MUSCLE","FETAL HEART")[ordIdx]], main=paste(paste("Exon", exonNo[selInd[i]]), 
          paste(paste(seqnames(obscnUniExGr)[selInd[i]], start(obscnUniExGr)[selInd[i]], sep=":"),
          end(obscnUniExGr)[selInd[i]], sep="-\n"),sep="\n"), 
                ylim=c(0,100),
                ylab=bquote(Psi*.(" (%)")),
                xlab="",border=cols[colInd], 
                names=c("AM","AH","FM","FH")[ordIdx], las=2)
        
        
        for(k in lowSamInd)
          points(x=rep(k,length(plotList[[k]])), y=plotList[[k]], col=cols[colInd][k], pch=1)
        
        
        for(k in 1:length(plotList))
          points(x=k, y=mean(plotList[[k]]), col= cols[colInd][k], pch=18, cex=1.5)
        
        
        mtext(subFig[i],  2, adj=5, padj = -15.5, las=1, line=0, font=2)
        
        for(l in 1:length(titleList)){
          if(l %% 3 == 0){
            if(length(grep("^\\*",titleList[[l-1]]))==0){
              title(xlab=bquote(Delta*Psi*.(titleList[[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
            } else {
              title(xlab=bquote(.("*")*Delta*Psi*.(titleList[[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
            }
          } else {
            title(xlab=titleList[[l]], line = 2+(l*1.5), col.lab = colTitles[l])
          }
        }
        
        
      }
    }
  }
}

dev.off()


####

pdf("figS1.pdf", width=21, height=35, pointsize=16)

ordIdx<-c(3,1,4,2)
colInd<- c(2,4,1,3)

cols<- rainbow(length(samTypes))
cols[4]<- "#FF80FF"

par(mar=c(31.1, 5.1, 6.1, 2.1))
par(cex.lab=1.7)
par(cex.axis=1.7)
par(cex.main=1.7)

xLim<- c(0,Inf)
yLim<- c(1,100)
if(yLim[1]==1){
  limPlot<- 125
  if(length(xLim[1])>0){
    strTmp<- "+"
    ordInd<- order(start(obscnUniExGr), decreasing=FALSE)
    selInd<- which(start(obscnUniExGr)>=xLim[1] & (start(obscnUniExGr)<=xLim[2]))
    exonNo<- selInd
    if(length(selInd)>0){
      if( strTmp== "-")
        exonNo<- length(start(obscnUniExGr))-selInd+1
      if(length(selInd)>limPlot)
        selInd<- selInd[1:limPlot]
      psidat[,selInd]
      exonNo<- exonNo[1:limPlot]
      
      par(mfrow=c(4,6))
      
      
      pvalsList<- list(difMuscleSel[selInd,"padj"],
                       difAdultSel[selInd,"padj"],
                       difListSel[[1]][selInd,"padj"],
                       difListSel[[2]][selInd,"padj"],
                       difListSel[[3]][selInd,"padj"],
                       difListSel[[4]][selInd,"padj"])
      
      delpsisList<- list(100*(rowMeans(psimean[selInd, c("ADULT MUSCLE", "FETAL MUSCLE")])-rowMeans(psimean[selInd, c("ADULT HEART", "FETAL HEART")])),
                         100*(rowMeans(psimean[selInd, c("ADULT MUSCLE", "ADULT HEART")])-rowMeans(psimean[selInd, c("FETAL MUSCLE", "FETAL HEART")])),
                         100*(psimean[selInd, c("ADULT MUSCLE")]-psimean[selInd, c("FETAL MUSCLE")]),
                         100*(psimean[selInd, c("ADULT HEART")]-psimean[selInd, c("FETAL HEART")]),
                         100*(psimean[selInd, c("ADULT MUSCLE")]-psimean[selInd, c("ADULT HEART")]),
                         100*(psimean[selInd, c("FETAL MUSCLE")]-psimean[selInd, c("FETAL HEART")]))
      
      selInd<- selInd[which( ((pvalsList[[1]]<0.05 & abs(delpsisList[[1]])>10) +
                                (pvalsList[[2]]<0.05 & abs(delpsisList[[2]])>10)+
                                (pvalsList[[3]]<0.05 & abs(delpsisList[[3]])>10)+
                                (pvalsList[[4]]<0.05 & abs(delpsisList[[4]])>10)+
                                (pvalsList[[5]]<0.05 & abs(delpsisList[[5]])>10)+
                                (pvalsList[[6]]<0.05 & abs(delpsisList[[6]])>10)) >= 1
      )]
      subFig<- LETTERS[1:length(selInd)]      
      for (i in 1:length(selInd)){
        plotList<- tapply(100*psidat[,selInd[i]], groupSam, c)
        grp<- rep("A",length(groupSam))
        grp[which(groupSam!=levels(groupSam)[1] & groupSam!=levels(groupSam)[2] )]<- "O"
        dattmp<- data.frame(exinc=psidat[,selInd[i]], group=grp)
        
        grp2<- rep("A",length(groupSam))
        grp2[grep("HEART",groupSam)]<- "O"
        dattmp2<- data.frame(exinc=psidat[,selInd[i]], group=grp2)
        

        pvalTxt1<- paste("FDR(M) = ", format(difMuscleSel[selInd[i],"padj"], digits =3), sep="")
        pvalTxt2<- paste("P(M) = ", format(difMuscleSel[selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt3<- paste("(M) = ", format(100*(mean(psimean[selInd[i], c("ADULT MUSCLE", "FETAL MUSCLE")])-mean(psimean[selInd[i], c("ADULT HEART", "FETAL HEART")])), digits =3), sep="")
        
        pvalTxt4<- paste("FDR(A) = ", format(difAdultSel[selInd[i],"padj"], digits =3), sep="")
        pvalTxt5<- paste("P(A) = ", format(difAdultSel[selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt6<- paste("(A) = ", format(100*(mean(psimean[selInd[i], c("ADULT MUSCLE", "ADULT HEART")])-mean(psimean[selInd[i], c("FETAL MUSCLE", "FETAL HEART")])), digits =3), sep="")
        
        pvalTxt7<- paste("FDR(",names(difListSel)[1],") = ", format(difListSel[[1]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt8<- paste("P(",names(difListSel)[1],") = ", format(difListSel[[1]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt9<- paste("(",names(difListSel)[1],") = ", format(100*(psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("FETAL MUSCLE")]), digits =3), sep="")
        
        pvalTxt10<- paste("FDR(",names(difListSel)[2],") = ", format(difListSel[[2]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt11<- paste("P(",names(difListSel)[2],") = ", format(difListSel[[2]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt12<- paste("(",names(difListSel)[2],") = ", format(100*(psimean[selInd[i], c("ADULT HEART")]-psimean[selInd[i], c("FETAL HEART")]), digits =3), sep="")
        
        pvalTxt13<- paste("FDR(",names(difListSel)[3],") = ", format(difListSel[[3]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt14<- paste("P(",names(difListSel)[3],") = ", format(difListSel[[3]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt15<- paste("(",names(difListSel)[3],") = ", format(100*(psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("ADULT HEART")]), digits =3), sep="")
        
        pvalTxt16<- paste("FDR(",names(difListSel)[4],") = ", format(difListSel[[4]][selInd[i],"padj"], digits =3), sep="")
        pvalTxt17<- paste("P(",names(difListSel)[4],") = ", format(difListSel[[4]][selInd[i],"pvalue"], digits =3), sep="")
        pvalTxt18<- paste("(",names(difListSel)[4],") = ", format(100*(psimean[selInd[i], c("FETAL MUSCLE")]-psimean[selInd[i], c("FETAL HEART")]), digits =3), sep="")
        
        pvalTxt1[difMuscleSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt1[difMuscleSel[selInd[i],"padj"]<0.05], sep="")
        pvalTxt2[difMuscleSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt2[difMuscleSel[selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt4[difAdultSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt4[difAdultSel[selInd[i],"padj"]<0.05], sep="")
        pvalTxt5[difAdultSel[selInd[i],"padj"]<0.05]<- paste("*", pvalTxt5[difAdultSel[selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt7[difListSel[[1]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt7[difListSel[[1]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt8[difListSel[[1]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt8[difListSel[[1]][selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt10[difListSel[[2]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt10[difListSel[[2]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt11[difListSel[[2]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt11[difListSel[[2]][selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt13[difListSel[[3]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt13[difListSel[[3]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt14[difListSel[[3]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt14[difListSel[[3]][selInd[i],"padj"]<0.05], sep="")
        
        
        pvalTxt16[difListSel[[4]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt16[difListSel[[4]][selInd[i],"padj"]<0.05], sep="")
        pvalTxt17[difListSel[[4]][selInd[i],"padj"]<0.05]<- paste("*", pvalTxt17[difListSel[[4]][selInd[i],"padj"]<0.05], sep="")
        
        
        
        pvals<- c(difMuscleSel[selInd[i],"padj"],
                  difAdultSel[selInd[i],"padj"],
                  difListSel[[1]][selInd[i],"padj"],
                  difListSel[[2]][selInd[i],"padj"],
                  difListSel[[3]][selInd[i],"padj"],
                  difListSel[[4]][selInd[i],"padj"])
        
        delpsis<- 100*c(mean(psimean[selInd[i], c("ADULT MUSCLE", "FETAL MUSCLE")])-mean(psimean[selInd[i], c("ADULT HEART", "FETAL HEART")]),
                        mean(psimean[selInd[i], c("ADULT MUSCLE", "ADULT HEART")])-mean(psimean[selInd[i], c("FETAL MUSCLE", "FETAL HEART")]),
                        psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("FETAL MUSCLE")],
                        psimean[selInd[i], c("ADULT HEART")]-psimean[selInd[i], c("FETAL HEART")],
                        psimean[selInd[i], c("ADULT MUSCLE")]-psimean[selInd[i], c("ADULT HEART")],
                        psimean[selInd[i], c("FETAL MUSCLE")]-psimean[selInd[i], c("FETAL HEART")])
        
        titleList<- list(pvalTxt1,
                         pvalTxt2,
                         pvalTxt3,
                         pvalTxt4,
                         pvalTxt5,
                         pvalTxt6,
                         pvalTxt7,
                         pvalTxt8,
                         pvalTxt9,
                         pvalTxt10,
                         pvalTxt11,
                         pvalTxt12,
                         pvalTxt13,
                         pvalTxt14,
                         pvalTxt15,
                         pvalTxt16,
                         pvalTxt17,
                         pvalTxt18)
        
        colTitles<- rep("black", length(titleList))
        colTitles[c(((which(pvals<0.05 & abs(delpsis)>10))*3),((which(pvals<0.05 & abs(delpsis)>10))*3-1),((which(pvals<0.05 & abs(delpsis)>10))*3-2))]<- "red"
        

        
        plotList2<- plotList
        lowSamInd<- which(unlist(lapply(plotList,length))<5)
        plotList2[lowSamInd]<- NA
        
        boxplot(plotList2[c("ADULT MUSCLE","ADULT HEART","FETAL MUSCLE","FETAL HEART")[ordIdx]], main=paste(paste("Exon", exonNo[selInd[i]]), 
                                                                                                           paste(paste(seqnames(obscnUniExGr)[selInd[i]], start(obscnUniExGr)[selInd[i]], sep=":"),
                                                                                                                 end(obscnUniExGr)[selInd[i]], sep="-\n"),sep="\n"), 
                ylim=c(0,100),
                ylab=bquote(Psi*.(" (%)")),
                xlab="",border=cols[colInd], 
                names=c("AM","AH","FM","FH")[ordIdx], las=2)
        
        for(k in lowSamInd)
          points(x=rep(k,length(plotList[[k]])), y=plotList[[k]], col=cols[colInd][k], pch=1)
        
        
        for(k in 1:length(plotList))
          points(x=k, y=mean(plotList[[k]]), col= cols[colInd][k], pch=18, cex=1.5)
        
        
        
        mtext(subFig[i],  2, adj=5, padj = -12, las=1, line=0, font=2)
        
        for(l in 1:length(titleList)){
          if(l %% 3 == 0){
            if(length(grep("^\\*",titleList[[l-1]]))==0){
              title(xlab=bquote(Delta*Psi*.(titleList[[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
            } else {
              title(xlab=bquote(.("*")*Delta*Psi*.(titleList[[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
            }
          } else {
            title(xlab=titleList[[l]], line = 2+(l*1.5), col.lab = colTitles[l])
          }
        }
        
        
      }
    }
  }
}

dev.off()

######
#Line plot



for(f in dir("./shinyData/", full.names = T))
  load(f)

pdf("fig3.pdf", width=15, height=10, pointsize=14)

ordIdx<-c(3,1,4,2)

cols<- rainbow(length(samTypes))
cols[4]<- "#FF80FF"


layout(matrix(c(1,1,1,1,2), nrow = 5, ncol = 1, byrow = TRUE))

psimeanBkup<- psimean
psidatBkup<- psidat

psimean<- psimean*100
psidat<- psidat*100

ranges <- list(x = c(min(start(uniqExGr))-100, max(start(uniqExGr))+100), y =c(0,0) )
strTmp<- "+"

xLim<- sort(ranges$x, decreasing=F)
yDist=.14
par(mar=c(yDist/2, 5.1, 4.1, 2.1))


plot(start(uniqExGr), psimean[,levels(groupSam)[1]], lwd=0, col=0, pch=NA,
     xlab="", xaxt='n', ylab=bquote(.("Mean ") ~ Psi ~ .(" (%)")), ylim=c(-5,100), main="OBSCN", xlim=xLim)
for(j in unique(c(firstExInd, lastExInd)))
  abline(v=start(uniqExGr)[j], lty=2, col="grey70", lwd=2)
if (length(dupInd)>0)
  for(j in dupInd)
    abline(v=start(uniqExGr)[j], lty=2, col="grey70", lwd=2)
points(start(uniqExGr), psimean[,(levels(groupSam))[ordIdx[1]]], type="b", pch=16, lwd=2, col=cols[1])
for(n in 2:length(cols))
  points(start(uniqExGr), psimean[,(levels(groupSam))[ordIdx[n]]], type="b", pch=16, lwd=2, col=cols[n])

if(strTmp=="+"){
  points(start(uniqExGr)[firstExInd], rep(-5, length(firstExInd)), col=1, bg=2, pch=24, cex=1.5, lwd=2)
  points(start(uniqExGr)[lastExInd], rep(-5, length(lastExInd)), col=1, bg=2, pch=25, cex=1.5, lwd=2)
} else {
  points(start(uniqExGr)[firstExInd], rep(-5, length(firstExInd)), col=1, bg=2, pch=25, cex=1.5, lwd=2)
  points(start(uniqExGr)[lastExInd], rep(-5, length(lastExInd)), col=1, bg=2, pch=24, cex=1.5, lwd=2)
}
if (length(dupInd)>0)
  points(start(uniqExGr)[dupInd], rep(-5, length(dupInd)), col=1, bg=4, pch=23, cex=1.5, lwd=2)


par(mar=c(5.1, 5.1,yDist/2, 2.1))
tmpPlotVarMed<-as.numeric(apply(psimean,1,var))

if(strTmp=="+"){
  xlabTmp<- "-->\nExon start"
} else{
  xlabTmp<- "<--\nExon start"
}
exonNo<- 1:length(start(uniqExGr))
if(strTmp=="-")
  exonNo<- length(start(uniqExGr))-(1:length(start(uniqExGr)))+1
exonNo<- as.character(exonNo)

xdel<- c(0, start(uniqExGr)[2:length(uniqExGr)]-start(uniqExGr)[1:(length(uniqExGr)-1)])
ydel<- c(0, abs(psimean[,levels(groupSam)[1]][-1]-psimean[,levels(groupSam)[1]][-c(length(uniqExGr))] ))

exonNo[which((xdel<3000) & ydel<10)]<- ""
exonNoTmp<- exonNo


i<-0
while(i<=length(exonNo)){
  i<- i+1
  j<- i
  if(exonNo[i]=="" & i<length(exonNo)){
    while(exonNo[j+1]=="" & (j+1)<length(exonNo)){
      j<- j+1
    }
    if(i!=j){
      exonNo[i+trunc((j-i)/2)]<- paste(i,j,sep="-")
    } else{
      exonNo[i+trunc((j-i)/2)]<- i
    }
    
  }
  i<-j
}

exonNoBkup<- exonNo
exonNoList<- strsplit(exonNo, split="-")

exonNo<- unlist(sapply(exonNoList, function(x){
  if(length(x)==0)
    x<- ""
  out=x
  if(length(x)==2){
    if(x[1]==x[2]){
      out<- x[1]
    } else{
      out<- paste(x[1:2], collapse="-")
    }
  } else{
    out<- x
  }
  return(out)
}))

text(start(uniqExGr), psimean[,levels(groupSam)[1]],  exonNo,
     cex=1, pos=1,col="black", srt=+30) 

strTmp<- "+"
xLim<- sort(ranges$x, decreasing=F)
yDist=.14

par(mar=c(5.1, 5.1,yDist/2, 2.1))
tmpPlotVarMed<-as.numeric(apply(psimean,1,var))

if(strTmp=="+"){
  xlabTmp<- "-->\nExon start"
} else{
  xlabTmp<- "<--\nExon start"
}

# Add legends
samNo<- table(groupSam)[levels(groupSam)]
names(samNo)<- levels(groupSam)
realNames<- paste(levels(groupSam)[ordIdx], " (N = ",as.numeric(table(groupSam)[levels(groupSam)][ordIdx]), ")")
legend(x=228355000,y=35, legend=c(realNames, "first exon", "Last exon"), lty=c(rep(1,length(unique(groupSam))),NA, NA), pch=c(rep(NA,length(unique(groupSam)) ),24,25),col=c(cols,"black","black"),
       pt.bg=c(rep(NA,length(unique(groupSam)) ), "red", "red"), lwd=rep(2, length(unique(groupSam))+2), pt.cex=c(rep(NA,length(unique(groupSam)) ), 1.5,1.5), bg="white")


plot(start(uniqExGr), tmpPlotVarMed, type="l", lwd=0, col=0, pch=NA, ylab="var", xlab=xlabTmp, xlim=xLim)
for(j in unique(c(firstExInd, lastExInd)))
  abline(v=start(uniqExGr)[j], lty=2, col="grey70", lwd=2)
if (length(dupInd)>0){
  for(j in dupInd)
    abline(v=start(uniqExGr)[j], lty=2, col="grey70", lwd=2)
}
points(start(uniqExGr), tmpPlotVarMed, type="l", pch=16, lwd=2, col="purple")

dev.off()

# THEN MODIFY SLIGHTLYT THE EXON NUMBERS !



psimean<- psimeanBkup
psidat<- psidatBkup

#### A3 and A5
load("./obscnUniExGr.rda")
obscnUniExGr
load(".//shinyData/firstExInd.rda")
load(".//shinyData/lastExInd.rda")
firstExInd
lastExInd


obscnUniExMap<- findOverlaps(obscnUniExGr,obscnUniExGr, type="any")
obscnUniExMapSta<- findOverlaps(obscnUniExGr,obscnUniExGr, type="start")
obscnUniExMapEnd<- findOverlaps(obscnUniExGr,obscnUniExGr, type="end")
obscnUniExMap<-obscnUniExMap[which(queryHits(obscnUniExMap)!= subjectHits(obscnUniExMap)),]
# Get rid f one instance when duplicate reverse data exists
rmInd<- c()
for (i in 1: length(obscnUniExMap)){
  tmpInd<- which(queryHits(obscnUniExMap)==subjectHits(obscnUniExMap)[i] & subjectHits(obscnUniExMap)==queryHits(obscnUniExMap)[i])
  if(length(tmpInd)>0 & tmpInd>i){
    rmInd<- c(rmInd,tmpInd)
  }
}
obscnUniExMap<- obscnUniExMap[-c(rmInd)]

obscnUniExMapSta<- obscnUniExMapSta[which(queryHits(obscnUniExMapSta)!= subjectHits(obscnUniExMapSta)),]
obscnUniExMapEnd<- obscnUniExMapEnd[which(queryHits(obscnUniExMapEnd)!= subjectHits(obscnUniExMapEnd)),]


obscnUniExMap<- obscnUniExMap[-c(which((queryHits(obscnUniExMap)%in% c(firstExInd,lastExInd, 126)) &
        (subjectHits(obscnUniExMap)%in% c(firstExInd,lastExInd, 126)))),]

a35Ind<- apply(as.data.frame(obscnUniExMap), 1, c, simplify = F)


## THESE ARE ONLY USEFUL FOR + STRAND. FOR - WRITE SUITABLE SIMILAR CODES
mapStaTmp1<-match(queryHits(obscnUniExMap),queryHits(obscnUniExMapSta))

mapStaTmp2<-match(subjectHits(obscnUniExMap),subjectHits(obscnUniExMapSta))



mapEndTmp1<-lapply(queryHits(obscnUniExMap), function(x) return(which(queryHits(obscnUniExMapEnd)==x)))
mapEndTmp1<- unlist(lapply(1:length(mapEndTmp1), function(x){lapply(mapEndTmp1[[x]], function(y){return(c(x,y))})}), recursive=FALSE)
mapEndTmp1<- matrix(unlist(mapEndTmp1), ncol=2, byrow=T)

mapEndTmp2<-lapply(subjectHits(obscnUniExMap), function(x) return(which(subjectHits(obscnUniExMapEnd)==x)))
mapEndTmp2<- unlist(lapply(1:length(mapEndTmp2), function(x){lapply(mapEndTmp2[[x]], function(y){return(c(x,y))})}), recursive=FALSE)
mapEndTmp2<- matrix(unlist(mapEndTmp2), ncol=2, byrow=T)

mapSta<- which((!is.na(mapStaTmp1))&(!is.na(mapStaTmp2))&(mapStaTmp1==mapStaTmp1))
mapEnd<- mapEndTmp2[which(mapEndTmp2[,1]==mapEndTmp1[,1] & mapEndTmp2[,2]==mapEndTmp1[,2]),1]
mapEnd<- mapEnd[which((!queryHits(obscnUniExMap)%in%firstExInd) & !(subjectHits(obscnUniExMap)%in%firstExInd))]

a3<- obscnUniExMap[mapEnd,]
a3GrList<- apply(as.data.frame(a3),1,function(x){return(obscnUniExGr[x,])})
leavea3In<-c() 


### SATRT OF PLOTTING OF A3s
for(n in 1:length(a3GrList)) {
print(n)
tmpGr<- a3GrList[[n]]
endTmpGr<- GRanges(seqnames=seqnames(tmpGr),
                IRanges(start=start(tmpGr)-1,end=start(tmpGr)-1), strand=strand(tmpGr))

allIsGr<- GRanges(seqnames=rowData(allIsObj)$chr, 
                  IRanges(
                    start=rowData(allIsObj)$begin,
                    end=rowData(allIsObj)$end),
                  strand=rowData(allIsObj)$strand)
allIsGrUniq<- unique(allIsGr)
uniqIsObj<- allIsObj[findOverlaps(allIsGrUniq,allIsGr, type="equal", select="first"),]

mapEnds<-findOverlaps(endTmpGr, allIsGrUniq, type="end",select="all")

leavea3InAdd<- FALSE
if(length(mapEnds)>1){
  leavea3IndAdd<- TRUE
  tmpCounts<-tapply(subjectHits(mapEnds),queryHits(mapEnds),function(x){return(colSums(counts(uniqIsObj)[x,,drop=FALSE]))})
  sumCounts<- matrix(unlist(tmpCounts), nrow=length(tmpCounts), byrow=T)
  (a3sumMat<-t(apply(sumCounts, 1,function(x){return(100*x/(1+colSums(sumCounts)))})))


  colnames(a3sumMat)<- names(tmpCounts[[1]])

  load(".//shinyData/groupSam.rda")
  load(".//shinyData/psidat.rda")

  table(rownames(psidat)==colnames(a3sumMat))

  rownames(a3sumMat)<-as.character(tmpGr)
  groupSamOrd<- factor(as.character(groupSam), levels=c("ADULT HEART", "FETAL HEART", "ADULT MUSCLE", "FETAL MUSCLE"))
  tempPlot<- aggregate(t(a3sumMat), list(groupSamOrd), c)

  listPlot<- unlist(c(tempPlot[,-1]), recursive=FALSE)
  samTypes<- tempPlot[,1]

  listPlotLab<- unlist(lapply(names(listPlot), function(x){
    stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
    scoord<- as.character(tmpGr)
    return(paste(substr(x,1,nchar(x)-1), " ",stype, sep=""))
  }))

  listPlotLabType<- unlist(lapply(names(listPlot), function(x){
    stype= samTypes[as.numeric(substr(x, nchar(x), nchar(x)))]
    return(stype)
  }))

  listPlotLabCoord<- unlist(lapply(names(listPlot), function(x){
    return(substr(x,1,nchar(x)-1))
  }))


  pdf(paste("AllPsiA3_",paste(as.character(tmpGr),collapse="_"),".pdf", sep=""), width=13, height=10, pointsize=14)
  par(mar=c(20.1, 4.1, 4.1, 2.1))
  
  ordIdx<-c(3,1,4,2)
  tmpOrdIdx<-c(2,4,1,3)
  cols<- rainbow(length(samTypes))
  cols[4]<- "#FF80FF"
  

  plotList2<- listPlot
  lowSamInd<- which(unlist(lapply(listPlot,length))<5)
  plotList2[lowSamInd]<- NA
  
  boxplot(plotList2, main="", 
          names=listPlotLabType, border=rep(cols[ordIdx], length(tmpGr)), las=2, ylim=c(0,100),
          ylab=bquote(Psi ~ .("(%)")))
  
  for(k in lowSamInd)
    points(x=rep(k,length(listPlot[[k]])), y=listPlot[[k]], col=rep(cols[ordIdx],length(tmpGr))[k], pch=1)
  
  
  for(k in 1:length(listPlot))
    points(x=k, y=mean(listPlot[[k]]), col= rep(cols[ordIdx],length(tmpGr))[k], pch=18, cex=1.5)
  
  resPlot<-c()
  resPlot$stats<- lapply(listPlot, function(x){return(as.matrix(summary(x)))})
  resPlot$stats<- do.call(cbind.data.frame, resPlot$stats)[-1,]
  
  resPlot$stats
  lapply(seq(from=4.5, by=4, length.out=length(listPlot)/4-1), function(x){abline(v=x, col="grey")})
  mains<- as.character(tmpGr)
  mapTmpObscn<-findOverlaps(tmpGr, obscnUniExGr, type="equal",select="all")
  mains<- paste("Exon ", subjectHits(mapTmpObscn), ": ", mains, sep="")
  for(i in 1:length(mains))
    mtext(text=mains[i], side = 3, line = 1, outer = FALSE, at = (i-1)*4+2+((i-1)*.5), font=2)


  delPsi<- resPlot$stats[3,seq(from=1,by=2,length.out=length(listPlot)/2)] - resPlot$stats[3,seq(from=2,by=2,length.out=length(listPlot)/2)]
  labs<- list()
  labs<- delPsi

  tmpNumb<- format(round(labs, 3), nsmall=3, digits=3)
  tmpNumb[which(round(labs, 4)==0)]<- format(labs[which(round(labs, 4)==0)], nsmall=4, digits=4, scientific=T)
  
  labelsPlot<- 
    paste(rep(c("DeltaPsi(AH/FH)", 
                "DeltaPsi(AM/FM)"),length(labs)/2),
          gsub(" ","",tmpNumb), sep="=")

  labPlot=labelsPlot


  for(i in 1:length(labelsPlot)){

    mtext(text=bquote(Delta*Psi*.(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 9, outer = FALSE, 
          at = seq(from=1.5, by=2, length.out=(2*length(mains)))[i], col="black")
  
  }
  par(xpd=TRUE)
  for(i in seq(from=1, by=2, length.out=4))
    drawBracket(x=c(i,i+1), y=-52)
  
  
  # second level lables
  
  delPsi<- resPlot$stats[3,seq(from=3,by=4,length.out=length(listPlot)/4)] - resPlot$stats[3,seq(from=1,by=4,length.out=length(listPlot)/4)]
  labs<- list()
  labs<- delPsi

  tmpNumb<- format(round(labs, 3), nsmall=3, digits=3)
  tmpNumb[which(round(labs, 4)==0)]<- format(labs[which(round(labs, 4)==0)], nsmall=4, digits=4, scientific=T)
  
  labelsPlot<- 
    paste(rep(c("DeltaPsi(AM/AH)"),length(labs)),
          gsub(" ","",tmpNumb), sep="=")

  labPlot=labelsPlot
  

  
  for(i in 1:length(labelsPlot)){

    mtext(text=bquote(Delta*Psi*.(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 13, outer = FALSE, 
          at = seq(from=2, by=4, length.out=(length(mains)))[i], col="black")
    
  }
  par(xpd=TRUE)
  for(i in seq(from=1, by=4, length.out=2))
    drawBracket(x=c(i,i+2), y=-75)
  
  
  # third level lables
  
  delPsi<- resPlot$stats[3,seq(from=4,by=4,length.out=length(listPlot)/4)] - resPlot$stats[3,seq(from=2,by=4,length.out=length(listPlot)/4)]
  labs<- list()
  labs<- delPsi
  
  tmpNumb<- format(round(labs, 3), nsmall=3, digits=3)
  tmpNumb[which(round(labs, 4)==0)]<- format(labs[which(round(labs, 4)==0)], nsmall=4, digits=4, scientific=T)
  
  labelsPlot<- 
    paste(rep(c("DeltaPsi(FM/FH)"),length(labs)),
          gsub(" ","",tmpNumb), sep="=")
  
  labPlot=labelsPlot
  
  
  
  for(i in 1:length(labelsPlot)){
    
    mtext(text=bquote(Delta*Psi*.(gsub("DeltaPsi","", labelsPlot[i]))), side = 1, line = 17, outer = FALSE, 
          at = seq(from=3, by=4, length.out=(length(mains)))[i], col="black")
    
  }
  par(xpd=TRUE)
  for(i in seq(from=2, by=4, length.out=2))
    drawBracket(x=c(i,i+2), y=-98)

  dev.off()
}
leavea3In<- c(leavea3In, leavea3InAdd)
}
a3GrList<- a3GrList[which(leaveA5In)]

#Rename figure file to FigS2
file.rename(
  "./AllPsiA5_chr1:228377931-228378006:+_chr1:228377937-228378006:+.pdf",
  "./figS2.pdf")

### Make table

obscnUniExGr

# Making index for Meta gene
obscnUniExMap<- findOverlaps(obscnUniExGr,obscnUniExGr, type="any")
obscnUniExMap<-obscnUniExMap[which(queryHits(obscnUniExMap)!= subjectHits(obscnUniExMap)),]
# Get rid f one instance when duplicate reverse data exists
rmInd<- c()
for (i in 1: length(obscnUniExMap)){
  tmpInd<- which(queryHits(obscnUniExMap)==subjectHits(obscnUniExMap)[i] & subjectHits(obscnUniExMap)==queryHits(obscnUniExMap)[i])
  if(length(tmpInd)>0 & tmpInd>i){
    rmInd<- c(rmInd,tmpInd)
  }
}
obscnUniExMap<- obscnUniExMap[-c(rmInd)]
ExClasses<- rep(0,length(obscnUniExGr))
classN<- 1
ExClasses[c(queryHits(obscnUniExMap)[1], subjectHits(obscnUniExMap)[1])]<- classN
for(i in 2:nrow(as.data.frame(obscnUniExMap))){
  elements<- as.numeric(as.data.frame(obscnUniExMap)[i,])
  found<- FALSE
  for(j in 1:classN){
    if(!found & length(which(which(ExClasses==j) %in% elements))>0){
      ExClasses[elements]=j
      found=TRUE
    }
  }
  if(!found){
    classN<- classN+1
    ExClasses[elements]=classN
  }
}   
   

metaId<- 1
for(i in 2:length(ExClasses)){
  if(ExClasses[i]!=0 & ExClasses[i]%in%ExClasses[1:(i-1)]){
    metaId[i]<- unique(metaId[which(ExClasses[1:i-1]==ExClasses[i])])
  } else {
    metaId[i]<- metaId[i-1]+1
  }
    
}

pasteChar<- rep("",length(ExClasses))
for(x in min(ExClasses[ExClasses!=0]):max(ExClasses)){
  indsClass<- which(ExClasses==x)
  indsCool<- indsClass[which(width(obscnUniExGr)[indsClass]==max(width(obscnUniExGr)[indsClass]))]
  if(length(indsClass[indsClass!=indsCool])>0)
    pasteChar[indsClass[indsClass!=indsCool]]<- letters[1:length(indsClass[indsClass!=indsCool])]
} 

metaId<- paste(metaId,pasteChar,sep="")

#

trFile="./OBSCN_ENSID.txt"
selectTr<- scan(trFile, what="character", sep = "\n")
selTr<- gsub("\\..*","",selectTr)

firstexBool<- rep(FALSE, length(obscnUniExGr))
firstexBool[firstExInd]<- TRUE
lastexBool<- rep(FALSE, length(obscnUniExGr))
lastexBool[lastExInd]<- TRUE
lastexBool[length(lastexBool)]<- TRUE
incBool=matrix(0, nrow=length(obscnUniExGr), ncol=length(selTr))

load("./obscnRefUncolExGr.rda")
load("./obscnRefUncol.rda")
obscnRefUncolEx<- obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),]
mapUniqEx<-findOverlaps(obscnUniExGr,obscnRefUncolExGr, type="equal",select="all")
mapUniqExTr<-tapply(queryHits(mapUniqEx),gsub("\\..*","",obscnRefUncolEx$transcript_id)[subjectHits(mapUniqEx)],
       c)
for(i in 1:length(mapUniqExTr))
  incBool[mapUniqExTr[[i]],i]<- 1
colnames(incBool)<- names(mapUniqExTr)
asVec<- rep("",length(obscnUniExGr))

for(j in 1:ncol(incBool)){
    incBool[which(!(incBool[,j]==1)),j]<- NA
    incBool[which(incBool[,j]==1),j]<- 1:length(which(incBool[,j]==1))
}
  
library(openxlsx)
for(i in 1:length(a5GrList)){
  grTmp<- a5GrList[[i]]
  tmpMap<-findOverlaps(grTmp,obscnUniExGr, type="equal")
  for(j in 1:nrow(as.data.frame(tmpMap)))
      asVec[subjectHits(tmpMap)[j]]<- paste(paste(asVec[subjectHits(tmpMap)[j]], ", A5'=", paste(as.character(grTmp)[-c(queryHits(tmpMap)[j])], collapse="/"), sep=""), collapse="")
}

allFirst<- c()
if(length(allFirst)>0){
  allFirstList<- tapply(allFirst, end(allFirst), c)
  allFirstList<- sapply(allFirstList, paste, collapse="/")
  for(i in 1:length(allFirstList)){
    grTmp<- allFirstList[i]
    tmpMap<-match(unlist(strsplit(grTmp, split="/")),as.character(obscnUniExGr))
    for(j in 1:length(tmpMap))
     asVec[tmpMap[j]]<- paste(paste(asVec[tmpMap[j]], ", AF=",as.character(allFirstList)[-c(i)], sep=""), collapse="")
  }
}

allLastList<- tapply(allLast, end(allLast), c)
allLastList<- sapply(allLastList, paste, collapse="/")
for(i in 1:length(allLastList)){
  grTmp<- allLastList[i]
  tmpMap<-match(unlist(strsplit(grTmp, split="/")),as.character(obscnUniExGr))
  for(j in 1:length(tmpMap))
    asVec[tmpMap[j]]<- paste(paste(asVec[tmpMap[j]], ", AL=",as.character(allLastList)[-c(i)], sep=""), collapse="")
}

asVec<- gsub("^, ","",asVec)
exonsOut<- data.frame("Exon number"=1:length(obscnUniExGr),
           "Meta transcript exon number"= metaId, chr=seqnames(obscnUniExGr),
           begin=start(obscnUniExGr), end=end(obscnUniExGr), strand=strand(obscnUniExGr),
           incBool,
           "First exon"=firstexBool,
           "Last exon"=lastexBool,
           "Alternative splicing"= asVec
             )
colnames(exonsOut)<- gsub("\\.", " ", colnames(exonsOut))



wb<- createWorkbook("Unique OBSCN exons")
addWorksheet(wb, "Unique OBSCN exons")
writeData(wb, sheet = 1, exonsOut, rowNames = FALSE, borders="all",
          borderColour="#000000")
rowColInd<-which(as.matrix(!is.na(exonsOut[,7:12])), arr.ind = TRUE)[,1]+1
colColInd<-which(as.matrix(!is.na(exonsOut[,7:12])), arr.ind = TRUE)[,2]+6

rowColInd<-c(rowColInd,which(as.matrix(exonsOut[,13:14]), arr.ind = TRUE)[,1]+1)
colColInd<-c(colColInd,which(as.matrix(exonsOut[,13:14]), arr.ind = TRUE)[,2]+12)

addStyle(wb, sheet = 1, createStyle(border = "TopBottomLeftRight", fgFill = "lightgreen"), 
         rows = rowColInd, cols = colColInd)
addStyle(wb, sheet = 1, createStyle(border = "TopBottomLeftRight", halign = "center",
                                    fgFill = "#999999", borderColour = "#000000"), 
         rows = 1, cols = 1:ncol(exonsOut))


setColWidths(wb, 1, cols = 1, widths = 35)
saveWorkbook(wb, "table2.xlsx", overwrite = TRUE)






###########################
###########################
#### Check exons thoroughly
load(file="./exObj.rda")
load(file="./resExEx.rda")
refEx<- rowData(resExEx)[which(rowData(resExEx)$int_ex=="exon"),]

library(biomaRt)
ensembl<-useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( 
  attributes=c("ensembl_gene_id", "external_gene_name"), 
  filters = "ensembl_gene_id", values =unique(rowData(exObj)$gene_id), 
  mart=ensembl, uniqueRows = TRUE)

genEns<- annoGen[,2]
names(genEns)<-annoGen[,1]

gen<- unique(rowData(exObj)$gene_id)
print(genEns[gen])
indsel<- which(rowData(exObj)$gene_id==gen)
allExGr<- uniqExGr
load(file="./shinyData/uniqExGr.rda")

innerMapUniqEx<- GenomicRanges::findOverlaps(uniqExGr,uniqExGr,
type="any") 
which(table(queryHits(innerMapUniqEx))>1)


mapUniqEx<- GenomicRanges::findOverlaps(uniqExGr, allExGr,
                                        type="any") 
(A53<- sort(table(queryHits(mapUniqEx)), decreasing=T)[which(sort(table(queryHits(mapUniqEx)), decreasing=T)>1)])


load(file="./shinyData/firstExInd.rda")
load(file="./shinyData/lastExInd.rda")

names(A53)[names(A53) %in% c(firstExInd,lastExInd)]

names(A53)[!(names(A53) %in% c(firstExInd,lastExInd))]


mapUniqEx[queryHits(mapUniqEx)==117,]


allExGr[subjectHits(mapUniqEx[queryHits(mapUniqEx)==117,]),]


uniqExGr[117,]


A5<- allExGr[subjectHits(mapUniqEx[queryHits(mapUniqEx)==117,]),]

save(A5, file="A5.rda")

load(file="allIsObj.rda")
allIsGr<- GRanges(as.character(rowData(allIsObj)$chr), 
                  IRanges(start=as.numeric(rowData(allIsObj)$begin),
                          end=as.numeric(rowData(allIsObj)$end)),
                  strand=as.character(rowData(allIsObj)$strand))

allIsUniqGr<- unique(allIsGr)