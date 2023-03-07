x<-read.table("/projects/bionano/tmp/test")
bin_segment<-function(x)
{
library(gtools)
library(tidyr)

colnames(x)<-c("mapID","chr","start_q","end_q","start_r","end_r","strand")
x<-x[order(x$mapID,x$chr,x$start_q),]
x$end_r<-as.numeric(x$end_r)
tmp<-x       
tmp_start<-min(tmp$start_r)
tmp_end<-max(tmp$end_r) 
#window<-ceiling(tmp_end/tmp_start) ##how to set the window size
window=10000
tmp_set<-seq(tmp_start,tmp_end,window)  

data_tmp <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = length(tmp_set),
                          ncol = 3))   
           
data_tmp$X1<-unique(tmp$chr)
data_tmp$X2<-tmp_set+1
data_tmp$X3<-tmp_set+window
data_tmp$X4<-paste0("A",(1:length(tmp_set)))
colnames(data_tmp)<-c("chr","start","end","note")
data_tmp$chr<-paste0("chr",data_tmp$chr)

##for each contig
mapid<-unique(tmp$mapID)
data_all<-c()
#one mapID have one contig
for (i in 1:length(mapid)){
	contig<-tmp[which(tmp$mapID==mapid[i]),c("chr","start_r","end_r","strand")]
	
	if(nrow(contig)==1){
	   contig_start<-min(contig$start_r)
       contig_end<-max(contig$end_r) 
       #window<-ceiling(tmp_end/tmp_start) ##how to set the window size??
       window=10000
       contig_set<-seq(contig_start,contig_end,window)
       
       ##make a dataframe with a window gap for each contig
       contig_tmp <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = length(contig_set),
                          ncol = 4))       
       contig_tmp$X1<-unique(contig$chr)
       contig_tmp$X2<-contig_set+1
       contig_tmp$X3<-contig_set+window
       contig_tmp$X4<-as.character(contig$strand)
       colnames(contig_tmp)<-c("chr","start","end","strand")
       contig_tmp$chr<-paste0("chr",contig_tmp$chr)
       contig_tmp$block_rank<-1 
       }      
 ##one mapID have multiple contigs      
   if(nrow(contig)>1){
  	   contig_tmp<-c()
      for(j in 1:nrow(contig)){
       	contig_start<-min(contig[j,]$start_r)
       contig_end<-max(contig[j,]$end_r) 
       #window<-ceiling(tmp_end/tmp_start) ##how to set the window size
       window=10000
       contig_set<-seq(contig_start,contig_end,window)
  
       ##make a dataframe with a window gap for each contig
       contig_tmp1 <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = length(contig_set),
                          ncol = 4))       
       contig_tmp1$X1<-unique(contig$chr)
       contig_tmp1$X2<-contig_set+1
       contig_tmp1$X3<-contig_set+window
       contig_tmp1$X4<-as.character(contig[j,]$strand)
       colnames(contig_tmp1)<-c("chr","start","end","strand")
       contig_tmp1$chr<-paste0("chr",contig_tmp1$chr)
       contig_tmp1$block_rank<-j                 ##add contig block order
       contig_tmp<-rbind(contig_tmp,contig_tmp1)
          }
       }
##make overlap
       refGR <- makeGRangesFromDataFrame(data_tmp,keep.extra.columns=TRUE,  ignore.strand=FALSE)
       testGR <- makeGRangesFromDataFrame(contig_tmp,keep.extra.columns=TRUE,  ignore.strand=FALSE)
       hits <- findOverlaps(refGR, testGR)
##Compute percent overlap and filter the hits
       overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
       percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
       hits <- hits[percentOverlap > 0.5]
       data_final<-contig_tmp[subjectHits(hits),]
       data_final$mapID<-mapid[i]
       data_final$note<-data_tmp[queryHits(hits),]$note
##reverse block order if strand is "-"
       blocks<-unique(data_final$block_rank)
          all<-c()
          for (k in 1:length(blocks)){
             	    block_k<-data_final[which(data_final$block_rank==k),]
             	    if(block_k$strand=="+"){
             	    block_k$note<-block_k$note
             	     } else {
             	    block_k$note<-rev(block_k$note)
             	     }
             	     all<-rbind(all,block_k)
         }
     data_all<-rbind(data_all,all)
     }  
