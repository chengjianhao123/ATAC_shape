if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("GO.db")
install.packages("car")
library(GO.db)
library(clusterProfiler)
library(org.Mm.eg.db)
library(car)
library(dplyr)
install.packages("reticulate")
library(reticulate)
library(reticulate)
source_python("pickle_reader.py")
library(ggplot2)
library(plyr)


#read json result
json_tmp<-list.files(path="uropa",pattern="finalhits.txt")
json_result<-paste0("uropa/",json_tmp)



##name the plot name and ego result
plot1name<-paste0("count",c(1:31),".png")
plot2name<-paste0("Peak_cluster",c(1:31),".png")
plot3name<-paste0("GT_cluster",c(1:31),".png")
ego<-list()
length(ego)<-31
sample_result<-list()
length(sample_result)<-31
pkl_name<-paste0("pkl",1:31)
for (k in 1:31){
  final<-read_pickle_file(pkl_name[k])
  final$`Unnamed: 0`<-NULL
  colnames(final)[14]<-c("cluster_labels")
  #colnames(final)[15]<-c("count")
  
  #calculate reads in each plot
  for (s in 1:nrow(final)){
    chr_count=final$group[s]
    gene_name=final$gene_rank2[s]
    final$count[s]=sum(myViews[[k]][[chr_count]][[gene_name]])
    print(s)
  }
  
  #plot count information
  plot0<-ggplot(final, aes(x=as.factor(cluster_labels), y=count,fill=as.factor(cluster_labels))) + 
    geom_boxplot() +
    labs(fill = "cluster") +
    xlab("cluster") +
    ggtitle("Count and Cluster relationship")+
    theme(plot.title = element_text(hjust = 0.5))
  plot0<-plot0 + scale_fill_manual(values=c("red", "blue", "green"))
  ggsave(paste0("../count/",plot1name[k]),plot0,width =5, height =5, dpi = 600)
  
  #connect to json_result
  upora_res<-read.table(json_result[[k]],header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
  tmp_res<-dplyr::select(upora_res,peak_id,gene_id)
  # remove duplicated result
  tmp_res=tmp_res[!duplicated(tmp_res$gene_id),]
  final<-dplyr::left_join(final,tmp_res,by="gene_id")
  #change into categorical data
  final$peak_id[!is.na(final$peak_id)]="peak"
  final$peak_id[is.na(final$peak_id)]="non_peak"
  final$peak_id=as.factor(final$peak_id)
  #build chi-square test
  ass_res1<-data.frame(matrix(0,nrow=6,ncol=4))
  ass_res1[,1]<-c(rep("1",2),rep("2",2),rep("3",2))
  ass_res1[,2]<-rep(c("non_peak","peak"),3)
  for (i in 1:3){
    a<-final%>%filter(cluster_labels==i)%>%filter(peak_id=="non_peak")%>%nrow()
    b<-final%>%filter(cluster_labels==i)%>%filter(peak_id=="peak")%>%nrow()
    ass_res1[(2*i-1),3]<-a
    ass_res1[(2*i),3]<-b
    ass_res1[(2*i-1),4]<-round(a/(a+b),2)
    ass_res1[(2*i),4]<-round(b/(a+b),2)
  }
  colnames(ass_res1)<-c("cluster","peak","number","percentile")
  plot1<-ggplot(ass_res1,aes(fill=peak,y=percentile,x=cluster)) + 
  geom_bar(position="fill",stat="identity")+ geom_text(aes(label=ifelse(percentile!=0,paste0(percentile*100,"%"),NA)),position=position_stack(vjust=0.5),size=5)+
    ggtitle("Peaks and Cluster association")+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("../peak_count/",plot2name[k]),plot1,width =5, height = 5, dpi = 600)
  
  
  #association with coding  genes
  ass_res2<-data.frame(matrix(0,nrow=6,ncol=4))
  ass_res2[,2]<-rep(c("cluster1","cluster2","cluster3"),2)
  ass_res2[,1]<-c(rep("coding_gene",3),rep("non_coding",3))
  for (i in 1:3){
    a<-final%>%filter(cluster_labels==i)%>%filter(label=="coding_gene")%>%nrow()
    b<-final%>%filter(cluster_labels==i)%>%filter(label=="noncoding_genes")%>%nrow()
    ass_res2[i,3]<-a
    ass_res2[(i+3),3]<-b
  }
  #calculate the precentage
  gene_sum<-aggregate(ass_res2$X3,by=list(ass_res2$X1),sum)
  ass_res2$X4=round(ass_res2$X3/c(rep(gene_sum[1,2],3),rep(gene_sum[2,2],3)),2)
  colnames(ass_res2)<-c("gene_type","cluster","number","percentile")
  plot2<-ggplot(ass_res2,aes(fill=cluster,y=percentile,x=gene_type)) + 
  geom_bar(position="fill",stat="identity")+geom_text(aes(label=paste0(percentile*100,"%")),position = position_stack(vjust = 0.5),size=5)+
  scale_fill_manual(values=c("red", "blue", "green"))+
  ggtitle("Cluster and Gene Type association")+
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("../gene_type/",plot2name[k]),plot2,width =5, height =5, dpi=600)
  
  
  
  #geo analysis for one_sample
  str_string <- function(x) {
    x<-unlist(strsplit(x, "\\."))
    y<-x[1]
    return(y)
  }
  gene_id2<-sapply(final$gene_id,str_string)
  names(gene_id2)<-NULL
  final$gene_id2<-gene_id2
  for (i in 1:3){
    cluster1gene=as.vector(unlist(filter(final,cluster_labels==i)%>%dplyr::select(gene_id2)))
    ego[[k]][[i]]<- enrichGO(gene  =cluster1gene,
                             universe = gene_table_final$gene_id2,
                             OrgDb         = org.Mm.eg.db,
                             keyType       = 'ENSEMBL',
                             ont           = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)}
  sample_result[[k]]<-final
  print(k)
}




