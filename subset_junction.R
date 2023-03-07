subset_junction<-function(all_final)
{
library(stringr)
tmp_seq<-sapply(all_final$seq,function(x) str_extract_all(x, "(A[0-9]+\\|A[0-9]+)"))
tmp_seq<-sapply(tmp_seq,function(x) paste0(x,collapse=","))
all_final$uncontinued_seq<-tmp_seq

all_final1<-all_final %>% 
    mutate(note = strsplit(as.character(uncontinued_seq), "\\,")) %>% 
    unnest(note) %>%
    as.data.frame %>%
    subset(select=c(ID,note))
colnames(all_final1)<-c("mapID","note")
}
