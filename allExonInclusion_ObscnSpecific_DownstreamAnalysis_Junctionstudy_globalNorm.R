
setwd("/PATH/TO/WORKFOLDER/")
JuncDirPath<-"/PATH/TO/JUNCTION_FILES/"

JuncDirs<- dir(path=JuncDirPath, pattern="_junction.bed", 
               include.dirs = T, full.names = T)


#Loading OBSCN exon coords 
load("./shinyData/uniqExGr.rda")

targCoords<- list(c(min(start(uniqExGr)), max(end(uniqExGr))),
                  c(min(start(uniqExGr)), max(end(uniqExGr))),
                  c(min(start(uniqExGr)), max(end(uniqExGr))))

aimTargCoords<- list(c(min(end(uniqExGr)[11:26]), max(start(uniqExGr)[11:26])),
                  c(min(end(uniqExGr)[45:59]), max(start(uniqExGr)[45:59])),
                  c(min(end(uniqExGr)[96:113]), max(start(uniqExGr)[96:113])))

#targGr<- 
tempDat<- list()
tempCnt<- list()
namesVec<-c()
for(i in 1:length(JuncDirs)){
  for(f in dir(pattern = ".bed",path=JuncDirs[i],full.names = TRUE)){
    print(f)
    namesVec<- c(namesVec, f)
    tmp<- read.table(f,sep="\t", stringsAsFactors = F, header=F)

    tmpStart<- as.numeric(tmp[,2])+as.numeric(sapply(strsplit(tmp[,11], split=","), head, 1))
    tmpEnd<- as.numeric(tmp[,3])- as.numeric(sapply(strsplit(tmp[,11], split=","), function(x){return(x[2])}))
    # Focus on End > Start cases !!!

    tmpGr<- GRanges(  seqnames= as.character(tmp[,1]),
                      IRanges::IRanges(
                        start=tmpStart ,
                        end=tmpEnd+1
                      ), strand=as.character(tmp[,6]))
    nReads<- as.numeric(tmp[,5])
    
    length(which((as.character(seqnames(tmpGr))%in%as.character(seqnames(uniqExGr))) &
          ((start(tmpGr))%in%end(uniqExGr))))
    # 159
    
    length(which((as.character(seqnames(tmpGr))%in%as.character(seqnames(uniqExGr))) &
                   ((end(tmpGr))%in%start(uniqExGr))))
    #  162
    if(i==1 & f==dir(pattern = ".bed",path=JuncDirs[i],full.names = TRUE)[1]){
      inds<- sapply(targCoords, function(x){
        return(which(as.character(seqnames(tmpGr))%in%seqnames(uniqExGr) &
                 start(tmpGr)>x[1] & end(tmpGr)<x[2]))
      }, simplify = F)
      tempDat<- list(tmpGr[inds[[1]]], tmpGr[inds[[2]]], tmpGr[inds[[3]]])
      tempCnt<- list(nReads[inds[[1]]], nReads[inds[[2]]], nReads[inds[[3]]])
    } else{
      inds<- sapply(targCoords, function(x){
        return(which(as.character(seqnames(tmpGr))%in%seqnames(uniqExGr) &
                       start(tmpGr)>x[1] & end(tmpGr)<x[2]))
      }, simplify = F)
      tempDat<- c(tempDat,list(tmpGr[inds[[1]]], tmpGr[inds[[2]]], tmpGr[inds[[3]]]))
      tempCnt<- c(tempCnt,list(nReads[inds[[1]]], nReads[inds[[2]]], nReads[inds[[3]]]))
    }

  }
}
#namesTmpBkup<- namesTmp
namesTmp<- gsub("_junction.bed","",gsub(".*\\/","",namesVec))
namesTmp<- gsub("_FKRN.*|_Musc","",namesTmp)
match(rownames(psidat), namesTmp)
length(which(is.na(match(rownames(psidat), namesTmp))))
#[1] 0
#cat(c(bamF1, bamF2), file=paste(JuncDirPath,"remainingBam.txt",sep="/"), sep = "\n")

exjnc<- list(unique(do.call(c,tempDat[seq(from=1, to=length(tempDat), by=3)])),
             unique(do.call(c,tempDat[seq(from=2, to=length(tempDat), by=3)])),
             unique(do.call(c,tempDat[seq(from=3, to=length(tempDat), by=3)])))
exjnccnt<-list() 
exjnccnt1<- c()
for(i in seq(from=1, to=length(tempDat), by=3)){
  cntTmp<- rep(0, length(exjnc[[1]]))
  cntTmp[findOverlaps(do.call(c,tempDat[i]), exjnc[[1]], type="equal", select="first")]<- unlist(tempCnt[i])
  if(length(exjnccnt1)==0){
    exjnccnt1=cntTmp} else{
      exjnccnt1=cbind(exjnccnt1,cntTmp)
      }
}

exjnccnt2<- c()
for(i in seq(from=2, to=length(tempDat), by=3)){
  cntTmp<- rep(0, length(exjnc[[2]]))
  cntTmp[findOverlaps(do.call(c,tempDat[i]), exjnc[[2]], type="equal", select="first")]<- unlist(tempCnt[i])
  if(length(exjnccnt2)==0){
    exjnccnt2=cntTmp} else{
      exjnccnt2=cbind(exjnccnt2,cntTmp)
    }
}

exjnccnt3<- c()
for(i in seq(from=3, to=length(tempDat), by=3)){
  cntTmp<- rep(0, length(exjnc[[3]]))
  cntTmp[findOverlaps(do.call(c,tempDat[i]), exjnc[[3]], type="equal", select="first")]<- unlist(tempCnt[i])
  if(length(exjnccnt3)==0){
    exjnccnt3=cntTmp} else{
      exjnccnt3=cbind(exjnccnt3,cntTmp)
    }
}

#exjnc<- exjncBkup
exjnccntFil1<- exjnccnt1[,match(rownames(psidat), namesTmp)] 
exjnccntFil2<- exjnccnt2[,match(rownames(psidat), namesTmp)] 
exjnccntFil3<- exjnccnt3[,match(rownames(psidat), namesTmp)]
exjncBkup<-exjnc

if(length(which(rowSums(exjnccntFil1)<=5))>0){
  indRow1<- which((rowSums(exjnccntFil1>5)>=2))
  exjnc[[1]]<- exjnc[[1]][indRow1]
  exjnccntFil1<- exjnccntFil1[indRow1,]
}

if(length(which(rowSums(exjnccntFil2)<=5))>0){
  indRow2<- which((rowSums(exjnccntFil2>5)>=2))
  exjnc[[2]]<- exjnc[[2]][indRow2]
  exjnccntFil2<- exjnccntFil2[indRow2,]
}

if(length(which(rowSums(exjnccntFil3)<=5))>0){
  indRow3<- which((rowSums(exjnccntFil3>5)>=2))
  exjnc[[3]]<- exjnc[[3]][indRow3]
  exjnccntFil3<- exjnccntFil3[indRow3,]
}

aEnd<- GRanges(  seqnames= as.character(seqnames(uniqExGr)),
                  IRanges::IRanges(
                    start=end(uniqExGr),
                    end=end(uniqExGr)
                  ), strand=as.character(strand(uniqExGr)))
aStart<- GRanges(  seqnames= as.character(seqnames(uniqExGr)),
                 IRanges::IRanges(
                   start=start(uniqExGr) ,
                   end=start(uniqExGr)
                 ), strand=as.character(strand(uniqExGr)))

psiFil1<- exjnccntFil1/matrix(rep(as.numeric(colSums(exjnccntFil1)), nrow(exjnccntFil1)), byrow=T, ncol=ncol(exjnccntFil1))
psiFil2<- exjnccntFil2/matrix(rep(as.numeric(colSums(exjnccntFil2)), nrow(exjnccntFil2)), byrow=T, ncol=ncol(exjnccntFil2))
psiFil3<- exjnccntFil3/matrix(rep(as.numeric(colSums(exjnccntFil3)), nrow(exjnccntFil3)), byrow=T, ncol=ncol(exjnccntFil3))

aimTargtInds<- sapply(aimTargCoords, function(y){
  tmpGr<- GRanges(seqnames= as.character(seqnames(uniqExGr)[1]),
  IRanges::IRanges(
    start=y[1] ,
    end=y[2]
  ), strand=as.character(strand(uniqExGr)[1]))
  return(queryHits(findOverlaps(exjnc[[1]], tmpGr, type="any")))
})
psiFil1Bkup<-psiFil1
psiFil2Bkup<-psiFil2
psiFil3Bkup<-psiFil3

psiFil1<- psiFil1[aimTargtInds[[1]],]
psiFil2<- psiFil2[aimTargtInds[[2]],]
psiFil3<- psiFil3[aimTargtInds[[3]],]

exjncBkup<- exjnc

exjnc[[1]]<- exjnc[[1]][aimTargtInds[[1]]]
exjnc[[2]]<- exjnc[[2]][aimTargtInds[[2]]]
exjnc[[3]]<- exjnc[[3]][aimTargtInds[[3]]]

typeA1End<- findOverlaps(exjnc[[1]], aEnd, type="start")
typeA1Sta<- findOverlaps(exjnc[[1]], aStart, type="end")
typeA1<- which((1:length(exjnc[[1]]))%in%queryHits(typeA1End) & (1:length(exjnc[[1]]))%in%queryHits(typeA1Sta))
typeB1<- which(!(1:length(exjnc[[1]]))%in%typeA1)

typeA2End<- findOverlaps(exjnc[[2]], aEnd, type="start")
typeA2Sta<- findOverlaps(exjnc[[2]], aStart, type="end")
typeA2<- which((1:length(exjnc[[2]]))%in%queryHits(typeA2End) & (1:length(exjnc[[2]]))%in%queryHits(typeA2Sta))
typeB2<- which(!(1:length(exjnc[[2]]))%in%typeA2)

typeA3End<- findOverlaps(exjnc[[3]], aEnd, type="start")
typeA3Sta<- findOverlaps(exjnc[[3]], aStart, type="end")
typeA3<- which((1:length(exjnc[[3]]))%in%queryHits(typeA3End) & (1:length(exjnc[[3]]))%in%queryHits(typeA3Sta))
typeB3<- which(!(1:length(exjnc[[3]]))%in%typeA3)

load("./shinyData/groupSam.rda")
load(file="./shinyData/groupSam.rda")
codeGrp<- c("FM", "AM", "FH", "AH")
names(codeGrp)<- levels(groupSam)
psiFilList<- list(psiFil1,psiFil2,psiFil3)

cols<- rainbow(length(unique(groupSam)))

ag1<- aggregate(t(100*psiFilList[[1]]),list(as.character(codeGrp[as.character(groupSam)])),mean)
annoJnc1<- rep("B",nrow(psiFilList[[1]]))
annoJnc1[typeA1]<- "A"
table(annoJnc1)

ag2<- aggregate(t(100*psiFilList[[2]]),list(as.character(codeGrp[as.character(groupSam)])),mean)
annoJnc2<- rep("B",nrow(psiFilList[[2]]))
annoJnc2[typeA2]<- "A"
table(annoJnc2)

ag3<- aggregate(t(100*psiFilList[[3]]),list(as.character(codeGrp[as.character(groupSam)])),mean)
annoJnc3<- rep("B",nrow(psiFilList[[3]]))
annoJnc3[typeA3]<- "A"
table(annoJnc3)

# PSI Filter: Filter to junctions which At least one sample exists with 
# ExJnc higher than 0.1%
quantile(psiFilList[[1]],.9)
# 90% 
# 0.0005270387 
quantile(psiFilList[[2]],.9)
# 90% 
# 0.0004246478
quantile(psiFilList[[3]],.9)
# 90% 
# 0.0007622646

length(which(rowSums(psiFilList[[1]]>.001)>=1))
# [1] 119
length(which(rowSums(psiFilList[[2]]>.001)>=1))
# [1] 49
length(which(rowSums(psiFilList[[3]]>.001)>=1))
# [1] 19

extFil1<- which(rowSums(psiFilList[[1]]>.001)>=1)
psiFilList[[1]]<- psiFilList[[1]][extFil1,]
annoJnc1<-annoJnc1[extFil1]
table(annoJnc1)
# A  B 
# 34 85
extFil2<- which(rowSums(psiFilList[[2]]>.001)>=1)
psiFilList[[2]]<- psiFilList[[2]][extFil2,]
annoJnc2<-annoJnc2[extFil2] 
table(annoJnc2)
# A  B 
# 28 21 
extFil3<- which(rowSums(psiFilList[[3]]>.001)>=1)
psiFilList[[3]]<- psiFilList[[3]][extFil3,]
annoJnc3<-annoJnc3[extFil3] 
table(annoJnc3)
# A  B 
# 17  2 

# ANALYZING FIRST REGION
exjncSt1<- GRanges(  seqnames= as.character(seqnames(exjnc[[1]][extFil1])),
                   IRanges::IRanges(
                     start=start(exjnc[[1]][extFil1]) ,
                     end=start(exjnc[[1]][extFil1])
                   ), strand=as.character(strand(exjnc[[1]][extFil1])))
exjncEn1<- GRanges(  seqnames= as.character(seqnames(exjnc[[1]][extFil1])),
                     IRanges::IRanges(
                       start=end(exjnc[[1]][extFil1]) ,
                       end=end(exjnc[[1]][extFil1])
                     ), strand=as.character(strand(exjnc[[1]][extFil1])))
exjncEn1<- findOverlaps(exjncEn1, uniqExGr, type="start")
exjncSt1<- findOverlaps(exjncSt1, uniqExGr, type="end")

exonsAnnoSt1<- tapply(paste("EX",subjectHits(exjncSt1),sep=""), queryHits(exjncSt1), paste, sep=",")
exonsAnnoEn1<- tapply(paste("EX",subjectHits(exjncEn1),sep=""), queryHits(exjncEn1), paste, sep=",")
exonsAnno1<- rep("", length(annoJnc1))
exonsAnno1[as.numeric(names(exonsAnnoSt1))]<- paste(exonsAnno1[as.numeric(names(exonsAnnoSt1))], exonsAnnoSt1, sep=",")
exonsAnno1[as.numeric(names(exonsAnnoEn1))]<- paste(exonsAnno1[as.numeric(names(exonsAnnoEn1))], exonsAnnoEn1, sep=",")


exjnc1<- GRanges(  seqnames= as.character(seqnames(exjnc[[1]][extFil1])),
                     IRanges::IRanges(
                       start=start(exjnc[[1]][extFil1]) ,
                       end=end(exjnc[[1]][extFil1])
                     ), strand=as.character(strand(exjnc[[1]][extFil1])))
exjncAny1<- findOverlaps(exjnc1, uniqExGr, type="any")
exAnnoAny1<- tapply(subjectHits(exjncAny1), queryHits(exjncAny1), 
  function(x){return(paste(c("",paste("ex",as.numeric(summary(x)[c(1,6)]), sep="")), collapse=","))})

exonsAnno1[which(annoJnc1=="B" & exonsAnno1=="")]<- exAnnoAny1[which(annoJnc1=="B" & exonsAnno1=="")]
exonsAnno1[which(annoJnc1=="B" & exonsAnno1!="")]<- 
  as.character(unlist(lapply(which(annoJnc1=="B" & exonsAnno1!=""), function(x){
    out<-gsub(gsub("EX","ex",exonsAnno1[x]), exonsAnno1[x], exAnnoAny1[x])
    return(out)
  })))




# ANALYZING SECOND REGION

exjncSt2<- GRanges(  seqnames= as.character(seqnames(exjnc[[2]][extFil2])),
                     IRanges::IRanges(
                       start=start(exjnc[[2]][extFil2]) ,
                       end=start(exjnc[[2]][extFil2])
                     ), strand=as.character(strand(exjnc[[2]][extFil2])))
exjncEn2<- GRanges(  seqnames= as.character(seqnames(exjnc[[2]][extFil2])),
                     IRanges::IRanges(
                       start=end(exjnc[[2]][extFil2]) ,
                       end=end(exjnc[[2]][extFil2])
                     ), strand=as.character(strand(exjnc[[2]][extFil2])))
exjncEn2<- findOverlaps(exjncEn2, uniqExGr, type="start")
exjncSt2<- findOverlaps(exjncSt2, uniqExGr, type="end")

exonsAnnoSt2<- tapply(paste("EX",subjectHits(exjncSt2),sep=""), queryHits(exjncSt2), paste, sep=",")
exonsAnnoEn2<- tapply(paste("EX",subjectHits(exjncEn2),sep=""), queryHits(exjncEn2), paste, sep=",")
exonsAnno2<- rep("", length(annoJnc2))
exonsAnno2[as.numeric(names(exonsAnnoSt2))]<- paste(exonsAnno2[as.numeric(names(exonsAnnoSt2))], exonsAnnoSt2, sep=",")
exonsAnno2[as.numeric(names(exonsAnnoEn2))]<- paste(exonsAnno2[as.numeric(names(exonsAnnoEn2))], exonsAnnoEn2, sep=",")



exjnc2<- GRanges(  seqnames= as.character(seqnames(exjnc[[2]][extFil2])),
                   IRanges::IRanges(
                     start=start(exjnc[[2]][extFil2]) ,
                     end=end(exjnc[[2]][extFil2])
                   ), strand=as.character(strand(exjnc[[2]][extFil2])))
exjncAny2<- findOverlaps(exjnc2, uniqExGr, type="any")
exAnnoAny2<- tapply(subjectHits(exjncAny2), queryHits(exjncAny2), 
                    function(x){return(paste(c("",paste("ex",as.numeric(summary(x)[c(1,6)]), sep="")), collapse=","))})

exonsAnno2[which(annoJnc2=="B" & exonsAnno2=="")]<- exAnnoAny2[which(annoJnc2=="B" & exonsAnno2=="")]
exonsAnno2[which(annoJnc2=="B" & exonsAnno2!="")]<- 
  as.character(unlist(lapply(which(annoJnc2=="B" & exonsAnno2!=""), function(x){
    out<-gsub(gsub("EX","ex",exonsAnno2[x]), exonsAnno2[x], exAnnoAny2[x])
    return(out)
  })))





#### ANALYZING THIRD REGION, It turned out to be an alternative last exon!


exjncSt3<- GRanges(  seqnames= as.character(seqnames(exjnc[[3]][extFil3])),
                     IRanges::IRanges(
                       start=start(exjnc[[3]][extFil3]) ,
                       end=start(exjnc[[3]][extFil3])
                     ), strand=as.character(strand(exjnc[[3]][extFil3])))
exjncEn3<- GRanges(  seqnames= as.character(seqnames(exjnc[[3]][extFil3])),
                     IRanges::IRanges(
                       start=end(exjnc[[3]][extFil3]) ,
                       end=end(exjnc[[3]][extFil3])
                     ), strand=as.character(strand(exjnc[[3]][extFil3])))
exjncEn3<- findOverlaps(exjncEn3, uniqExGr, type="start")
exjncSt3<- findOverlaps(exjncSt3, uniqExGr, type="end")

exonsAnnoSt3<- tapply(paste("EX",subjectHits(exjncSt3),sep=""), queryHits(exjncSt3), paste, sep=",")
exonsAnnoEn3<- tapply(paste("EX",subjectHits(exjncEn3),sep=""), queryHits(exjncEn3), paste, sep=",")
exonsAnno3<- rep("", length(annoJnc3))
exonsAnno3[as.numeric(names(exonsAnnoSt3))]<- paste(exonsAnno3[as.numeric(names(exonsAnnoSt3))], exonsAnnoSt3, sep=",")
exonsAnno3[as.numeric(names(exonsAnnoEn3))]<- paste(exonsAnno3[as.numeric(names(exonsAnnoEn3))], exonsAnnoEn3, sep=",")


exjnc3<- GRanges(  seqnames= as.character(seqnames(exjnc[[3]][extFil3])),
                   IRanges::IRanges(
                     start=start(exjnc[[3]][extFil3]) ,
                     end=end(exjnc[[3]][extFil3])
                   ), strand=as.character(strand(exjnc[[3]][extFil3])))
exjncAny3<- findOverlaps(exjnc3, uniqExGr, type="any")
exAnnoAny3<- tapply(subjectHits(exjncAny3), queryHits(exjncAny3), 
                    function(x){return(paste(c("",paste("ex",as.numeric(summary(x)[c(1,6)]), sep="")), collapse=","))})

exonsAnno3[which(annoJnc3=="B" & exonsAnno3=="")]<- exAnnoAny3[which(annoJnc3=="B" & exonsAnno3=="")]
exonsAnno3[which(annoJnc3=="B" & exonsAnno3!="")]<- 
  as.character(unlist(lapply(which(annoJnc3=="B" & exonsAnno3!=""), function(x){
    out<-gsub(gsub("EX","ex",exonsAnno3[x]), exonsAnno3[x], exAnnoAny3[x])
    return(out)
  })))







#####

save(bfinalPlotList, file="bfinalPlotList.rda")
save(bfinalPm, file="bfinalPm.rda")
save(bfinalPh, file="bfinalPh.rda")
save(bfinalMain, file="bfinalMain.rda")
save(bfinalDelM, file="bfinalDelM.rda")
save(bfinalDelH, file="bfinalDelH.rda")
save(bfinalTitleList, file="bfinalTitleList.rda")


save(finalPlotList, file="finalPlotList.rda")
save(finalPm, file="finalPm.rda")
save(finalPh, file="finalPh.rda")
save(finalMain, file="finalMain.rda")
save(finalDelM, file="finalDelM.rda")
save(finalDelH, file="finalDelH.rda")
save(finalTitleList, file="finalTitleList.rda")

subFig<- LETTERS 

# PLOTTING INTERESTING EVENTS FROM FIRST REGION (ON 5' END), 
# SECOND REGION (IN THE MIDDLE FO GENE) AND THRIRD REGION (ON THE 3' END).

k<- 0
pdf("./Fig6.pdf", height=8, width=10, pointsize = 9)
cols<- rainbow(length(unique(groupSam)))
cols[4]<- "#FF80FF"
par(mfrow=c(3,3))
layout(mat = matrix(c(1,2,3,4,5,6,8,7,8), 
                    ncol = 3, byrow = TRUE)
)
par(mar=c(10.1, 4.1, 4.1, 2.1))
tmpOrdIdx<-c(2,4,1,3)
library(clinfun)

chooseInd<- unlist(lapply(gsub("−","-", gsub(",","",gsub(" \\(.*","",c("chr1:228,256,926−228,273,272 (connecting exons 11 and 18)",
  "chr1:228,244,571−228,256,651 (connecting exons 15 and 19)",
  "chr1:228,288,888−228,292,523 (connecting exons 47 and 50)",
  "chr1:228,288,888−228,293,353 (connecting exons 47 and 59)", 
  "chr1:228,293,518−228,294,318 (connecting exons 51 and 56)", 
  "chr1:228,294,415−228,294,780 (connecting exons 52 and 59)",
  "chr1:228353062-228362576 (Ex96, 99)")))), function(x){return(grep(x,unlist(finalMain), perl=T))}))
colTitles<- rep(1, length(chooseInd))
for(i in chooseInd){
  k<- k+1
  
  
  
  
  plotList2<- finalPlotList[[i]]
  lowSamInd<- as.numeric(which(unlist(lapply(finalPlotList[[i]],length))<5))
  plotList2[lowSamInd]<- NA
  #names(plotList2)<- names(plotList)
  boxplot(plotList2,
          border=cols[tmpOrdIdx], main=finalMain[[i]], xlab="", ylab="EJ(%)")
  
  for(z in lowSamInd)
    points(x=rep(z,length(finalPlotList[[i]][[z]])), y=finalPlotList[[i]][[z]], col=cols[tmpOrdIdx][z], pch=1)
  
  
  for(z in 1:length(finalPlotList[[i]]))
    points(x=z, y=mean(finalPlotList[[i]][[z]]), col= cols[tmpOrdIdx][z], pch=18, cex=1.5)
  
  
  
  
  for(l in 1:length(finalTitleList[[i]])){
    if(l %% 2 == 0){
      if(length(grep("^\\*",finalTitleList[[i]][[l-1]]))==0){
        title(xlab=bquote(Delta*E*J*.(finalTitleList[[i]][[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
      } else {
        title(xlab=bquote(.("*")*Delta*E*J*.(finalTitleList[[i]][[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
      }
    } else {
      title(xlab=finalTitleList[[i]][[l]], line = 2+(l*1.5), col.lab = colTitles[l])
    }
    
  }
  mtext(subFig[k],  2, adj=4, padj = -10, las=1, line=0, font=2)
}
dev.off()



pdf("./Fig7.pdf", height=6, width=10, pointsize = 9)

cols<- rainbow(length(unique(groupSam)))
cols[4]<- "#FF80FF"
par(mfrow=c(2,3))
layout(mat = matrix(c(1,2,3,4,5,6), 
                    ncol = 3, byrow = TRUE)
)
par(mar=c(10.1, 4.1, 4.1, 2.1))
tmpOrdIdx<-c(2,4,1,3)
library(clinfun)

chooseInd<- unlist(lapply(gsub("−","-", gsub(",","",gsub(" \\(.*","",c("chr1:228,243,458−228,246,528 (overlapping exons16-18)",
  "chr1:228,305,386−228,306,372 (overlapping exons 54-56)", 
  "chr1:228,292,161−228,294,152 (overlapping exons 56 and 57)", 
  "chr1:228,295,020−228,300,123 (overlapping exons 53-56)", 
  "chr1:228,300,020−228,303,814 (overlapping exons 45-106)",
  "chr1:228,369,162−228,369,942 (overlapping exons 45-106)")))), function(x){return(grep(x,unlist(bfinalMain), perl=T))}))
colTitles<- rep(1, length(chooseInd))
k<- 0
for(i in chooseInd){
  k<- k+1
  # boxplot(bfinalPlotList[[i]],
  #         border=cols[tmpOrdIdx], main=bfinalMain[[i]], xlab="", ylab="EJ(%)")
  
  
  
  plotList2<- bfinalPlotList[[i]]
  lowSamInd<- as.numeric(which(unlist(lapply(bfinalPlotList[[i]],length))<5))
  plotList2[lowSamInd]<- NA
  #names(plotList2)<- names(plotList)
  boxplot(plotList2,
          border=cols[tmpOrdIdx], main=bfinalMain[[i]], xlab="", ylab="EJ(%)")
  
  for(z in lowSamInd)
    points(x=rep(z,length(bfinalPlotList[[i]][[z]])), y=bfinalPlotList[[i]][[z]], col=cols[tmpOrdIdx][z], pch=1)
  
  
  for(z in 1:length(bfinalPlotList[[i]]))
    points(x=z, y=mean(bfinalPlotList[[i]][[z]]), col= cols[tmpOrdIdx][z], pch=18, cex=1.5)
  
  
  
  for(l in 1:length(bfinalTitleList[[i]])){
    if(l %% 2 == 0){
      if(length(grep("^\\*",bfinalTitleList[[i]][[l-1]]))==0){
        title(xlab=bquote(Delta*E*J*.(bfinalTitleList[[i]][[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
      } else {
        title(xlab=bquote(.("*")*Delta*E*J*.(bfinalTitleList[[i]][[l]])), line = 2+(l*1.5), col.lab = colTitles[l])
      }
    } else {
      title(xlab=bfinalTitleList[[i]][[l]], line = 2+(l*1.5), col.lab = colTitles[l])
    }
    
  }
  
  mtext(subFig[k],  2, adj=4, padj = -11.5, las=1, line=0, font=2)
  
}
dev.off()