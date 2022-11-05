####DESEQ2 data############
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
##################################

######## select packages from tab manualy
######### the condition row should have the same name for each in a set each time, ignoring the identifier number
######### change file names,  res needed (paiwise comparisons)


library("DESeq2")
library("ggplot2")
library('reshape')
#setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Transcriptomic_Data/Euc-Mut_Raw_Data/pre_processed_raw_point_data')
#setwd('E:/WSU Work/Transcriptomic_Data/Euc-Mut_Raw_Data/pre_processed_raw_point_data')
setwd('E:/WSU Work/Transcriptomic_Data/Fungi-PIS-Mut_Raw_Data-Revisited')

#data=read.table("SI14_euc-pis_preprocessed_raw_point_timecourse_data.txt", header = TRUE, row.names = 1)
data=read.table("Piso_SI_14_Updated_Count_Timecourse_Pre-Proccessed.txt", header = TRUE, row.names = 1)

setwd('E:/WSU Work/Transcriptomic_Data/Fungi-PIS-Mut_Raw_Data-Revisited/Comparison_Construction_Files')


#data<-data.frame(format(round(data))) # needs interger values
#write.table(data,"SI14_euc-pis_preprocessed_raw_point_timecourse_data.txt")
colnames(data)
sample=colnames(data)

condition =c("Ctrl_plant","Ctrl_plant","Ctrl_plant",
             "precontact_24hr","precontact_24hr","precontact_24hr",
             "postcontact_24hr","postcontact_24hr","postcontact_24hr",
             "postcontact_48hr","postcontact_48hr","postcontact_48hr",
             "postcontact_1wk","postcontact_1wk","postcontact_1wk",
             "postcontact_2wk","postcontact_2wk","postcontact_2wk")

#"GeneID",
colData=cbind(colnames(data), sample, condition)

######  verification step
colData=data.frame(colData)
write.table(colData, file="Condition_verification.txt", sep="\t")

######boxplot of data transformation
epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(data + epsilon)), breaks=100, col="blue", border="white",
     main="Expressed_Log2-transformed counts per gene_sans20", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
dev.copy(png,'Log2-transformed_counts_expressed_sans20.png')
dev.off()
condition_1 <- as.data.frame(t(condition))
boxplot(log2(data + epsilon), condition_1$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(counts +1)")
dev.copy(png,'boxplot_expressed_transformed-data_sans20.png')
dev.off()

######dds=data deseq
dds= DESeqDataSetFromMatrix(data, colData, ~ condition)
dds= DESeq(dds)
results(dds)
str(dds)

###### normalized data output ############
norm_counts=counts(dds, normalized=TRUE)
write.table(norm_counts, file="Normalized_expressed_sans20.txt", sep="\t")
##rlog tranformation
rld<- rlogTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)

library("gplots")
library('RColorBrewer')

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(10, 10))
dev.copy(png,'Expressed_sample-to-sample_heatmap_sans20.png')
dev.off()
print(plotPCA(rld))
dev.copy(png, 'Expressed_seq_PCA')
dev.off()
##### know the size of the file, col number and lines) ###
dim(norm_counts)
str(dds)

#####Filtering 10 reads per sample repetition->excel########
##############Hclust##########
## log2(x+1) ##
norm_counts_log=log2(norm_counts+1)

### calculation of Pearson correlation ##

r = cor(norm_counts_log, method = "pearson")
write.table(r, file="Correlation_values_sans20.txt", sep="\t")
## dissimilarity calculation ##
d=1-r
## Clustering ##
h = hclust(as.dist(d), method = "complete", members = NULL)
plot(h)
dev.copy(png, 'Expressed_seq_hclust_sans20')
dev.off()

####comparaison####
####Pairwise comparison#####

res1 <- results(dds, contrast=c("condition","precontact_24hr","Ctrl_plant"))     ### test group first, control group second
res2 <- results(dds, contrast=c("condition","postcontact_24hr","Ctrl_plant"))
res3 <- results(dds, contrast=c("condition","postcontact_48hr","Ctrl_plant"))
res4 <- results(dds, contrast=c("condition","postcontact_1wk","Ctrl_plant"))
res5 <- results(dds, contrast=c("condition","postcontact_2wk","Ctrl_plant"))


###### Formatting #########
res1=res1[,c(2,6)]
res2=res2[,c(2,6)]
res3=res3[,c(2,6)]
res4=res4[,c(2,6)]
res5=res5[,c(2,6)]


#res12=res12[,c(2,6)]

#Final
colnames(res1)=c("Precontact_vs_Control_log2FC","Precontact_vs_Control_Padj")
colnames(res2)=c("Contact_24hr__vs_Control_log2FC","Contact_24hr__vs_Control_Padj")
colnames(res3)=c("Contact_48hr__vs_Control_log2FC","Contact_48hr__vs_Control_Padj")
colnames(res4)=c("Contact_1week__vs_Control_log2FC","Contact_1week__vs_Control_Padj")
colnames(res5)=c("Contact_2week__vs_Control_log2FC","Contact_2week__vs_Control_Padj")


Diff_Expressed_Genes=cbind(res1, res2, res3, res4, res5)


### Create output file in desegnated output folder

setwd('E:/WSU Work/Transcriptomic_Data/Fungi-PIS-Mut_Comparison_Data-Revisited')


#write.csv(Diff_Expressed_Genes, file="Euc-Pis_Compared_Relitive-Timecourse_DGE-Checked.csv", sep="\t")
write.csv(Diff_Expressed_Genes, file="Pisolithus-Genes_Euc-Pis_Compared_Relitive-Timecourse_DGE-Checked.csv", sep="\t")


