process_bins<-function(data_all)
{
un_blocks<-unique(data_all$mapID)
all_final<-c()
##set gap between bins
gap_blocks<-5 
for (i in 1:length(un_blocks)){
	un_contig<-data_all[which(data_all$mapID==un_blocks[i]),c("mapID","block_rank","note")]
	un_contig$block_rank<-as.factor(as.character(un_contig$block_rank))
	##split each block with a list
	un_contig_list<-split(un_contig,un_contig$block_rank)
	##mapID has one contig
	if(length(un_contig_list)==1){
		df_tmp<-un_contig_list[[1]]
        df_tmp$num<-as.numeric(gsub("A","",df_tmp$note))
        df_tmp$diff<-ave(df_tmp$num, factor(df_tmp$block_rank), FUN=function(x) c(NA,diff(x)))
        block_gap<-which(df_tmp$diff>gap_blocks)
           if(length(block_gap)>0){
                df_block<-df_tmp$note[block_gap] 
                df_block_forward<-df_tmp$note[block_gap-1]
                df_block_forward<-paste0("-",df_block_forward,"//") ##add connect symbol 
                df_gap_merge<-paste0(c(rbind(df_block_forward,df_block)),collapse="")               ##gather gap blocks
                df_all<-paste0(df_tmp$note[1],df_gap_merge,"-",df_tmp$note[nrow(df_tmp)])           ##add head and tail bins
		        all<-data.frame(ID=un_blocks[i],seq=df_all)
		        }else{
		        	 h1<-df_tmp$note[1] ##extract first row for one block
			         h2<-df_tmp$note[nrow(df_tmp)] ##extract last row for one block
			         block_h<-paste(h1,h2,sep="-") #combine first and last block
		        	 all<-data.frame(ID=un_blocks[i],seq=block_h)}    ##no gap blocks
                     }
        else{  
       ##mapID has multiple contigs
		combine<-c()
		for(h in 1:length(un_contig_list)){
            n_row<-nrow(un_contig_list[[h]])
            df_tmp<-un_contig_list[[h]]
            df_tmp$num<-as.numeric(gsub("A","",df_tmp$note))
            df_tmp$diff<-ave(df_tmp$num, factor(df_tmp$block_rank), FUN=function(x) c(NA,diff(x)))
            block_gap<-which(df_tmp$diff>gap_blocks)
            if(length(block_gap)>0){
            df_block<-df_tmp$note[block_gap] 
            df_block_forward<-df_tmp$note[block_gap-1]
            df_block_forward<-paste0("-",df_block_forward,"//")                                 ##add connect symbol 
            df_gap_merge<-paste0(c(rbind(df_block_forward,df_block)),collapse="")               ##gather gap blocks
            df_all<-paste0(df_tmp$note[1],df_gap_merge,"-",df_tmp$note[nrow(df_tmp)])           ##add head and tail bins
			combine<-paste(c(combine,df_all),collapse="|")                                            ##combine all blocks
	        } else {
	        	 n_row<-nrow(un_contig_list[[h]])
			     h1<-un_contig_list[[h]]$note[1] ##extract first row for one block
			     h2<-un_contig_list[[h]]$note[n_row] ##extract last row for one block
			     block_h<-paste(h1,h2,sep="-") #combine first and last block
			     combine<-paste(c(combine,block_h),collapse="|") #combine all blocks
	             }  
	        all<-data.frame(ID=un_blocks[i],seq=combine)
	            }		   
}	
 all_final<-rbind(all_final,all)    
}
}
