# It is assumed that the working directory is where the files of the scripts are
workDir= "/PATH/TO/WORKFOLDER/AND/SCRIPTS/FILES"

setwd(workDir)
# The object files necessary for the R/shiny based interactive 
# visualization (OBSCN_PSI) is stored in the 'shinyData' folder
dir.create("./shinyData/")

resFilesPath="./res/"


######### Make the data for Shiny FOR UNADJUSTED PSI levels


# Defining the unique exons
ProcessExons(trFile="./OBSCN_ENSID.txt",
             gtfFile="PATH/TO/GENE/ANNOTATION/GTF/FILE",
             geneName="OBSCN",
             outPath="./")
  

#Load all saved objects
for(f in dir("./","\\.rda$"))
  load(f)

#  Make result objects and measure PSI Exon Inclusion
MeasurePSIUnadj(datapath="./",
                outPath="./",
                uniqExGr=uniqExGr)
#Load all saved objects
for(f in dir("./","\\.rda$"))
  load(f)
# 
DifExInc(outPath="./",
         uniqExGr=uniqExGr,
         uniqExObj=uniqExObj,
         allExObj)

# finding first and last exons
trFile<- "./OBSCN_ENSID.txt"
selectTr<- scan(trFile, what="character", sep = "\n")
selectTr<- gsub("\\..*","", selectTr)

obscnExObj<- allExObj[which(gsub("\\..*","",rowData(allExObj)$gene_id)=="ENSG00000154358"),]
obscnFirstIndTmp<- tapply(1:length(obscnExObj), rowData(obscnExObj)$transcript_id, function(x){return(x[which(rowData(obscnExObj)$int_ex_num[x]==1)])})
obscnLastIndTmp<- tapply(1:length(obscnExObj), rowData(obscnExObj)$transcript_id, function(x){return(x[which(rowData(obscnExObj)$int_ex_num[x]==max(rowData(obscnExObj)$int_ex_num[x]))])})

obscnFirstInd<- as.numeric(obscnFirstIndTmp[match(selectTr,gsub("\\..*","",names(obscnFirstIndTmp)))])
obscnLastInd<- as.numeric(obscnLastIndTmp[match(selectTr,gsub("\\..*","",names(obscnLastIndTmp)))])

obscnFirstGr<- unique(GRanges(  seqnames= as.character(rowData(obscnExObj)[,"chr"]),
                           IRanges::IRanges(
                             start=as.numeric(rowData(obscnExObj)[,"begin"]) ,
                             end=as.numeric(rowData(obscnExObj)[,"end"])
                           ), strand=as.character(rowData(obscnExObj)[,"strand"]))[obscnFirstInd,])

obscnLastGr<- unique(GRanges(  seqnames= as.character(rowData(obscnExObj)[,"chr"]),
                                IRanges::IRanges(
                                  start=as.numeric(rowData(obscnExObj)[,"begin"]) ,
                                  end=as.numeric(rowData(obscnExObj)[,"end"])
                                ), strand=as.character(rowData(obscnExObj)[,"strand"]))[obscnLastInd,])

### Deallig with the two isoforms with CDS5' and CDS3' incomplete 
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

obscnFirstGr<- obscnFirstGr[-c(queryHits(findOverlaps(obscnFirstGr, rmZeroGr, type="equal"))),]

obscnLastGr<- obscnLastGr[-c(queryHits(findOverlaps(obscnLastGr, rmZeroGr, type="equal"))),]

load(file="shinyData/uniqExGr.rda")

firstExInd<- unique(queryHits(findOverlaps(uniqExGr,obscnFirstGr, type="equal")))
lastExInd<- unique(queryHits(findOverlaps(uniqExGr,obscnLastGr, type="equal")))

save(firstExInd, file="shinyData//firstExInd.rda")
save(lastExInd, file="shinyData//lastExInd.rda")

dupInd<- c()
save(dupInd, file="shinyData//dupInd.rda")


exObj<- obscnExObj[gsub("\\..*","",rowData(obscnExObj)$transcript_id) %in% selectTr,]
save(exObj, file="exObj.rda")
