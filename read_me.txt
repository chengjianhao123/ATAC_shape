The source code part contain 4 parts
1.the gene_table file.Based on the gene ontolgy V1 file download from Genecode,we calculated each gene's start poistion and end position.The result were written in csv files
2the raw data part,after download the raw data from GEO database, we normalized the data using normalized.bigwig.py module,and we then upload the normalized ATAC-seq data.Due to space limitation we only upload one normalized file
3.Shape analysis part were written in python.First we load the raw data to R.then passing each gene's read to python.then use python to do the helliger distance calculation and K-means clustering procedures.
3.Other characteristic.after getting the K-means result, we pass the result back to R and then calcuate each cluster's other characteristics.we also provide the the Peak annotation result in the uropa files.
