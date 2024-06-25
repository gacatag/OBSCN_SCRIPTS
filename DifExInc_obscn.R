#### run differential tests with batch adjustment
DifExInc<- function(outPath,
         obscnUniExGr,
         allExObj){

  
# First and Last exons should be zero as they don't have skipped exons.

uniqExObjGr<- GRanges(seqnames= as.character(rowData(uniqExObj)[,"chr"]),
                      IRanges::IRanges(
                        start=as.numeric(rowData(uniqExObj)[,"begin"]) ,
                        end=as.numeric(rowData(uniqExObj)[,"end"])),
                      strand= as.character(rowData(uniqExObj)[,"strand"]))


firstLastExInd<- as.numeric(unlist(tapply(1:length(allExObj), rowData(allExObj)[,"transcript_id"] ,function(x){return(c(min(x),max(x)))})))
firstLastExGr<- GRanges(seqnames= as.character(rowData(allExObj)[firstLastExInd,"chr"]),
                        IRanges::IRanges(
                          start=as.numeric(rowData(allExObj)[firstLastExInd,"begin"]) ,
                          end=as.numeric(rowData(allExObj)[firstLastExInd,"end"])),
                        strand= as.character(rowData(allExObj)[firstLastExInd,"strand"]))

# ENST00000660857 isoform of OBSCN is not completely annotated! It is flagged in 
#ENSEMBL as CDS5' and CDS3' incomplete. It's first and last exons should not be
#taken into account, i.e. their counts must not be set to Zero. We want to have 
#exon inclusion levels for its first and last exons.
# In the same manner, the first exon of ENST00000493977 isoform should also be 
#removed as it is CDS5' incomplete.

rowData(allExObj)[gsub("\\..*","",rowData(allExObj)$transcript_id)=="ENST00000660857",]
# DataFrame with 4 rows and 8 columns
# chr     begin       end      strand int_ex_num      int_ex
# <character> <numeric> <numeric> <character>  <integer> <character>
#   1        chr1 228309484 228309559           +          1        exon
# 2        chr1 228313614 228313886           +          3        exon
# 3        chr1 228314814 228315080           +          5        exon
# 4        chr1 228315847 228316119           +          7        exon
# transcript_id            gene_id
# <character>        <character>
#   1 ENST00000660857.1 ENSG00000154358.23
# 2 ENST00000660857.1 ENSG00000154358.23
# 3 ENST00000660857.1 ENSG00000154358.23
# 4 ENST00000660857.1 ENSG00000154358.23
rowData(allExObj)[gsub("\\..*","",rowData(allExObj)$transcript_id)=="ENST00000493977",]
# DataFrame with 4 rows and 8 columns
# chr     begin       end      strand int_ex_num      int_ex
# <character> <numeric> <numeric> <character>  <integer> <character>
#   1        chr1 228215689 228215829           +          1        exon
# 2        chr1 228216421 228216702           +          3        exon
# 3        chr1 228217013 228217288           +          5        exon
# 4        chr1 228217395 228218222           +          7        exon
# transcript_id            gene_id
# <character>        <character>
#   1 ENST00000493977.2 ENSG00000154358.23
# 2 ENST00000493977.2 ENSG00000154358.23
# 3 ENST00000493977.2 ENSG00000154358.23
# 4 ENST00000493977.2 ENSG00000154358.23


(rmZero<-rowData(allExObj)[c(which(gsub("\\..*","",rowData(allExObj)$transcript_id)=="ENST00000660857"&(rowData(allExObj)$int_ex_num==1 | rowData(allExObj)$int_ex_num==7)),
           which(gsub("\\..*","",rowData(allExObj)$transcript_id)=="ENST00000493977"&rowData(allExObj)$int_ex_num==1)), 1:4])
# DataFrame with 3 rows and 4 columns
# chr     begin       end      strand
# <character> <numeric> <numeric> <character>
#   1        chr1 228309484 228309559           +
#   2        chr1 228315847 228316119           +
#   3        chr1 228215689 228215829           +

rmZeroGr<- GRanges(seqnames= as.character(rmZero[,"chr"]),
                   IRanges::IRanges(
                     start=as.numeric(rmZero[,"begin"]) ,
                     end=as.numeric(rmZero[,"end"])),
                   strand= as.character(rmZero[,"strand"]))
save(rmZeroGr, file="rmZeroGr.rda")
zeroInd<- unique(queryHits(findOverlaps(uniqExObjGr, firstLastExGr, type="equal")))

overlapRm<-findOverlaps(uniqExObjGr[zeroInd,], rmZeroGr, type="equal")

zeroInd<- zeroInd[-c(queryHits(overlapRm))]

assays(uniqExObj)$"counts"[zeroInd, ] <- 0

selUniqExObj<- uniqExObj[,which(gsub("\\.1$","",gsub("_EEES.tsv","",rownames(colData(uniqExObj)))) %in% rownames(psidat))]


  

groupSamOrd<-groupSam[match(gsub("\\.1$","",gsub("_EEES.tsv","",rownames(colData(selUniqExObj)))),rownames(psidat))]
levels(groupSamOrd)
#[1] "ADULT MUSCLE"  "ADULT HEART"   "EMBRYO MUSCLE" "EMBRYO HEART"
cond<- as.character(groupSamOrd)
cond[grep("ADULT", cond)]<- "test"
cond[grep("EMBRYO", cond)]<- "ctrl"

selUniqExObj<- addAnnotation(x=selUniqExObj,
                         sampleAnnotationType="CONDITION",
                         sampleAnnotation=cond)

# Similar to DESeq differential expression analysis, you can also correct for 
# BATCH by using design=~BATCH+CONDITION+CONDITION:IMIS
difAdult<- deseqInterest(
  selUniqExObj,
  design=~batch+CONDITION+CONDITION:EEES, 
  contrast=list("CONDITIONtest.EEESee", "CONDITIONctrl.EEESee"), 
  bpparam = SnowParam(workers=20))


colData(selUniqExObj)<- colData(selUniqExObj)[, -c(ncol(colData(selUniqExObj)))]

cond<- as.character(groupSamOrd)
cond[grep("MUSCLE", cond)]<- "test"
cond[grep("HEART", cond)]<- "ctrl"

selUniqExObj<- addAnnotation(x=selUniqExObj,
                                sampleAnnotationType="CONDITION",
                                sampleAnnotation=cond)

# Similar to DESeq differential expression analysis, you can also correct for 
# BATCH by using design=~BATCH+CONDITION+CONDITION:IMIS
difMuscle<- deseqInterest(
  selUniqExObj,
  design=~batch+CONDITION+CONDITION:EEES, 
  contrast=list("CONDITIONtest.EEESee", "CONDITIONctrl.EEESee"), 
  bpparam = SnowParam(workers=20))



levels(groupSamOrd)
#[1] "FETAL MUSCLE" "ADULT MUSCLE" "FETAL HEART"  "ADULT HEART" 

el1<- c("AM", "AH", "AM", "FM")
el2<-c("FM","FH","AH", "FH")

elVec<- levels(groupSamOrd)
names(elVec)<- c("FM","AM","FH","AH")


# Running differential exon inclusion for comparisons AM/FM, AH/FH,
# AM/AH and FM/FH
colData(selUniqExObj)<- colData(selUniqExObj)[, -c(ncol(colData(selUniqExObj)))]

difList<- list()
batchAdjusted=c()
for(i in 1:length(el1)){


  cond<- rep("no", length(as.character(groupSamOrd)))
  cond[which(as.character(groupSamOrd)==elVec[el2[i]])]<- "ctrl"
  cond[which(as.character(groupSamOrd)==elVec[el1[i]])]<- "test"
  selUniqExObjTmp<- selUniqExObj[,which(cond!="no")]
  #colData(selUniqExObj)<- colData(selUniqExObj)[, -c(ncol(colData(selUniqExObj)))]
  
  selUniqExObjTmp<- addAnnotation(x=selUniqExObjTmp,
                                  sampleAnnotationType="CONDITION",
                                  sampleAnnotation=cond[which(cond!="no")])
  paste(colData(selUniqExObjTmp)$batch,colData(selUniqExObjTmp)$CONDITION, sep="/")
  if(i!=2){
    difTmp<- deseqInterest(
      selUniqExObjTmp,
      design=~batch+CONDITION+CONDITION:EEES, 
      contrast=list("CONDITIONtest.EEESee", "CONDITIONctrl.EEESee"), 
      bpparam = SnowParam(workers=20))
    batchAdjusted=c(batchAdjusted, TRUE)
  } else{
    difTmp<- deseqInterest(
      selUniqExObjTmp,
      design=~CONDITION+CONDITION:EEES, 
      contrast=list("CONDITIONtest.EEESee", "CONDITIONctrl.EEESee"), 
      bpparam = SnowParam(workers=20))
    batchAdjusted=c(batchAdjusted, FALSE)
  }
  difList<- c(difList,list(difTmp))
  
}

names(difList)<- paste(el1,el2, sep="/")



uniqAllExGr<- GRanges(
  seqnames= as.character(rowData(uniqExObj)[,"chr"]),
  IRanges::IRanges(
    start=as.numeric(rowData(uniqExObj)[,"begin"]) ,
    end=as.numeric(rowData(uniqExObj)[,"end"])
  ), strand=rowData(uniqExObj)[,"strand"])

obsInd<- findOverlaps(uniqExGr, uniqAllExGr, type="equal", select= "first")


difMuscleSel<- as.data.frame(difMuscle[obsInd,])
difAdultSel<- as.data.frame(difAdult[obsInd,])
difListSel<- lapply(1:length(difList), function(x){
  y<- as.data.frame(difList[[x]])
  y<- y[obsInd,]
  return(y)})
names(difListSel)<- names(difList)
save(difList, file="difList.rda")
save(difMuscle, file="difMuscle.rda")
save(difAdult, file="difAdult.rda")

save(difListSel, file="shinyData//difListSel.rda")
save(difMuscleSel, file="shinyData//difMuscleSel.rda")
save(difAdultSel, file="shinyData//difAdultSel.rda")
save(obsInd, file="shinyData//obsInd.rda")


}





