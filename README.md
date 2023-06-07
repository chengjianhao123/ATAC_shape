# ATAC_shape
This ATAC_shape project contains three different part.we list them in order
1.the gene table part.This part contains the how to calcuate the starting and ending position of a gene based on GENEONTOLOGY V1 annotation
file downloaded from the GENECODE project.the result were listed in csv file.
2.The clustering procedure.This part contains how to load the normalized bigwig file to R,then passing the each gene's read to python 
to do the dimension reduction and clustering procedure.we do not upload normalized bigwig file due to file size limitation.
3 combining clustering to other gene charcteristic (Total count and gene type and Peak association and RNA analysis).this task also need
the peak annotation result and rna-seq result,which were also omit due to file size limitation.
