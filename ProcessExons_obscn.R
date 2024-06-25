processExons<- function(trFile,
             gtfFile,
             geneName="OBSCN",
             outPath){
setwd(outPath)
# Defining the unique exons
selectTr<- scan(trFile, what="character", sep = "\n")
selectGen<- rep(geneName, length(selectTr))
selectTr<- gsub("\\..*","", selectTr)

library(IntEREst)
allRefUncol<- referencePrepare (sourceBuild="file",
                                filePath=gtfFile,
                                collapseExons=FALSE)
obscnRefUncol<- allRefUncol[which(gsub("\\..*","",allRefUncol[, "transcript_id"])%in%selectTr), ]
limitRange<- GRanges(
  seqnames= as.character(obscnRefUncol$chr)[1],
  IRanges::IRanges(
    start=min(as.numeric(obscnRefUncol$begin))-5 ,
    end=max(as.numeric(obscnRefUncol$end))+5
  ))


save(limitRange, file="limitRange.rda")
save(obscnRefUncol, file="obscnRefUncol.rda")

obscnRefUncolExGr<- GRanges(
  seqnames= as.character(obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),"chr"]),
  IRanges::IRanges(
    start=as.numeric(obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),"begin"]) ,
    end=as.numeric(obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),"end"])
  ), strand=obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),"strand"])

 

obscnUniExGr<- unique(obscnRefUncolExGr)
#length(obscnUniExGr)
#[1] 126
obscnUniExGr<- obscnUniExGr[order(start(obscnUniExGr), decreasing=FALSE),]

allFirstInd<- which(obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),"int_ex_num"]==1)
allLastInd<- as.numeric(tapply(1:length(which(obscnRefUncol$int_ex=="exon")),
                    obscnRefUncol[which(obscnRefUncol$int_ex=="exon"),"transcript_id"],
                    max))
allFirstGr<- unique(obscnRefUncolExGr[allFirstInd,])  
allLastGr<- unique(obscnRefUncolExGr[allLastInd,])  

mapObscnUniEx<- findOverlaps(obscnUniExGr,obscnUniExGr, type="any")

mapObscnUniEx<-mapObscnUniEx[which(queryHits(mapObscnUniEx)!= subjectHits(mapObscnUniEx)),]
# Get rid f one instance when duplicate reverse data exists
rmInd<- c()
for (i in 1: length(mapObscnUniEx)){
  tmpInd<- which(queryHits(mapObscnUniEx)==subjectHits(mapObscnUniEx)[i] & subjectHits(mapObscnUniEx)==queryHits(mapObscnUniEx)[i])
  if(length(tmpInd)>0 & tmpInd>i){
    rmInd<- c(rmInd,tmpInd)
  }
  
}
mapObscnUniEx<- mapObscnUniEx[-c(rmInd)]
mapObscnUniEx
# Hits object with 6 hits and 0 metadata columns:
#   queryHits subjectHits
# <integer>   <integer>
#   [1]         1           2
# [2]         1           3
# [3]         2           3
# [4]         8           9
# [5]       122         123
# [6]       125         126
# -------
#   queryLength: 126 / subjectLength: 126

ObscnUniFirstEx<- findOverlaps(obscnUniExGr,allFirstGr, type="equal")
ObscnUniLastEx<- findOverlaps(obscnUniExGr,allLastGr, type="equal")
allOverlappingInds<- unique(c(queryHits(mapObscnUniEx), subjectHits(mapObscnUniEx)))
(AltExInd<-sapply(list(AFirst=queryHits(ObscnUniFirstEx),
            ALast=queryHits(ObscnUniLastEx),
            A35=(1:length(obscnUniExGr))[which(!((1:length(obscnUniExGr)) %in% c(queryHits(ObscnUniFirstEx),queryHits(ObscnUniLastEx))))]),
       function(x){ return(allOverlappingInds[allOverlappingInds %in% x])}))
# $AFirst
# [1] 1 2 3 9
# 
# $ALast
# [1] 125 126
# 
# $A35
# [1]   8 122 123

sapply(AltExInd, function(x){return(start(obscnUniExGr)[x])})
# $AFirst
# [1] 228208063 228208130 228208044 228215689
# 
# $ALast
# [1] 228378618 228378618
# 
# $A35
# [1] 228215563 228377931 228377937

uniqExGr<- obscnUniExGr
(sameFirst<- table(start(uniqExGr))[table(start(uniqExGr))>1])
# 228378618 
# 2

omittedExon=c()
for(x in as.numeric(names(sameFirst))){
  indsSameFirst<- which(start(uniqExGr)==x)
  omittedExonTmp<- indsSameFirst[indsSameFirst!= min(indsSameFirst)]
  omittedExon<- c(omittedExon, omittedExonTmp)
  uniqExGr<- uniqExGr[which(!((1:length(uniqExGr)) %in% 
                                omittedExonTmp)),]
}


dir.create("shinyData")
save(uniqExGr,file="shinyData/uniqExGr.rda")
save(obscnUniExGr, file="obscnUniExGr.rda")
save(AltExInd, file="AltExInd.rda")
save(ObscnUniFirstEx, file="ObscnUniFirstEx.rda")
save(ObscnUniLastEx, file="ObscnUniLastEx.rda")
save(omittedExon, file="omittedExon.rda")
}
