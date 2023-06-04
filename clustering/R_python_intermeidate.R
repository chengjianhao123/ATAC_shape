library(rtracklayer)
library(IRanges)


#load the data to Myviews
q<-list.files("test/",pattern=".bigwig",recursive = T)
sample_table[sample_table$strain=="C57BL/6J",]


#my Views
myViews<-NULL
for (i in 1:31){
  cov.bw=rtracklayer::import(q[i], which=Ipromoter.gr,as = "RleList")
  myViews[[i]]=Views(cov.bw,Ipromoter.gr)
}


#export normalize reads to python
norm_reads_matrix<-paste0("norm_reads_matrix",c(1:31),".csv")
gene_name_var<-paste0("gene_name_text",c(1:31),".csv")
for (i in 1:31){
  reads<-matrix(0,ncol=2001,nrow=36706)
  for (j in 1:36706){
    gene_name=gene_table_final$gene_rank2[j]
    chr_count=gene_table_final$group[j]
    reads[j,]<-as.numeric(myViews[[i]][[chr_count]][[gene_name]])
    print(j)
  }
  reads_sum<-apply(reads,1,sum)
  NA_row<-which(reads_sum==0)
  gene_name_text<-gene_table_final[-NA_row,]
  reads<-reads[-NA_row,]
  write.csv(gene_name_text,gene_name_var[i])
  write.csv(reads,norm_reads_matrix[i])
}
