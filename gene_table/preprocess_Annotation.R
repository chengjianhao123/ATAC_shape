install.packages("tidyverse")
library(tidyverse)
library(rtracklayer)
library(IRanges)

#GT version
GT<-as.data.frame(rtracklayer::import("/Users/cheng/Desktop/bioinfo/gencode.mouse.v1.annotation.gtf"))
GT2<-GT[GT$type=="gene",]
GT2$gene_type<-as.factor(GT2$gene_type)
summary(GT2$gene_type)


#add new label
GT2$label<-NULL
#gtf3 file new_label
GT2$label[GT2$gene_type!="protein_coding"]<-c("noncoding_genes")
GT2$label[GT2$gene_type=="protein_coding"]<-c("coding_gene")
gene_table<-select(GT2,seqnames,start,gene_id,gene_name,gene_type,label)%>%mutate(tss_s=start-5000,tss_e=start+5000)




#make the gene_table
gene_expression$name<-rownames(gene_expression)
gene_table<-gene_table %>%inner_join(gene_expression,by=c("gene_name"="name"))%>%filter(seqnames!="chrM")
str(gene_table)
gene_table$gene_type<-as.factor(gene_table$gene_type)
summary(gene_table$gene_type)


# check the integrity of chrom.size
chrom.size<-read.table("/Users/cheng/Desktop/bioinfo/norm/mm9.chrom.sizes.txt",header = F)
colnames(chrom.size)<-c("seqnames","base_length")
chrom.size<-chrom.size[1:21,]
#check tss_s
tmp0<-which(gene_table$tss_s<0)

#check tss_e
tmp3<-NULL
chrom_length<-NULL
genet_length<-NULL
length(gene_table$seqnames)
gene_table$seqnames=as.character(gene_table$seqnames)
for (i in 1:length(gene_table$seqnames)){
  chrom_length<-chrom.size$base_length[which(gene_table$seqnames[i]==chrom.size$seqnames)]
  genet_length<-gene_table$tss_e[i]
   if (genet_length>chrom_length){
    tmp3<-append(tmp3,i)
 }
  else{
    print(i)
  }
}

final_filter<-c(tmp0,tmp3)
final_filter
gene_table2<-gene_table
str(gene_table2)

#promoter.gr
range<-IRanges(start=gene_table2$tss_s,end=gene_table2$tss_e)
promoter.gr<-GRanges(seqnames=gene_table2$seqnames,ranges=range,strand=rep("*",20789),names=gene_table2$gene_name)
names(promoter.gr)<-gene_table2$gene_id
length(promoter.gr)

# filter out blacklist regions
library(rtracklayer)
blacklistregion<-rtracklayer::import("/Users/cheng/Desktop/bioinfo/mm9-blacklist.bed",format="bed")
hits <- findOverlaps(promoter.gr, blacklistregion)
promoter.gr.filtered<-promoter.gr[setdiff(1:length(promoter.gr),queryHits(hits))]
promoter.gr.filtered
Ipromoter.gr<-as(promoter.gr.filtered,"IRangesList")
Ipromoter.gr$chrM<-NULL
#after fliter gene_table
gene_table3<-gene_table2[match(names(promoter.gr.filtered),gene_table2$gene_id,),]
gene_table4<-gene_table3[gene_table3$seqnames!="chrM",]
coding_num<-which(gene_table4$label=="coding_gene")
non_coding_num<-which(gene_table4$label=="noncoding_genes")
#final_table
str(gene_table4)
length(coding_num)
length(non_coding_num)

#load the file
myViews<-NULL
q<-list.files(recursive = T)
for (i in 1:11){
  cov.bw=import(q[i], which=promoter.gr.filtered,as = "RleList")
  myViews[[i]]=Views(cov.bw,Ipromoter.gr)
}

