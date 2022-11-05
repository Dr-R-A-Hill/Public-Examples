#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

browseVignettes("clusterProfiler")


### Setting work directory
setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Transcriptomic_Data/Euc-Mut_Comparison_Data')

reddeer<-read.csv("SI14_Euc_Comb_Timecourse_DGE(1).csv",fill=TRUE, header = TRUE)#xlsx


### Orginise and retrieve a list of genes 

reddeer5<-subset(reddeer,Week_2_Contact_vs_Control_Sig == "TRUE")
antilope5<-reddeer5[order(reddeer5$Week_2_Contact_vs_Control_Padj),]

mspacman<-antilope5[0:3073,"pacId"]                     ### gives pacId
pacman<-as.character(mspacman)

panthera<-antilope5[0:3073,"Panther"]                     ### gives Panther
panther<-as.character(panthera)
felidae<-data.frame(panther)
grrrrr<-c(felidae)

gogogadget<-antilope5[0:3073,"GO"]                     ### gives GO
gogadgetgo<-as.character(gogogadget)
gogetem<-data.frame(gogadgetgo)

itsako<-antilope5[0:3073,"KO"]                     ### gives KO
absoluteko<-as.character(itsako)
onehitko<-data.frame(absoluteko)



bestguess<-antilope5[0:3073,"Best_hit_arabi_name"]                     ### gives ?
goodguess<-as.character(bestguess)




### Export list as file

##setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Transcriptomic_Data/Euc-Mut_Raw_Data/Processed_Data')
setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Networks/Eucalyptus - Mutulist Networks/2 Week Post-Contact Condition')


#write.table(pacman,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_pacId.txt")
#write.table(grrrrr,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_Panther3.txt")
#write.csv(panther,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_Panther.csv")
#write.table(gogadgetgo,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_GO.txt")
#write.csv(gogadgetgo,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_GO.csv")
#write.table(absoluteko,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_KO.txt")
write.csv(goodguess,"3073_sig_euc-SI14_genes_mut_2week_rallconditions_Best_hit_arabi_name.csv")

## ----style, echo=FALSE, results="asis", message=FALSE--------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE)



if (!requireNamespace("goana", quietly = TRUE))
  install.packages("goana")

## Linear model usage:

fit <- lmFit(gogadgetgo, design)
fit <- eBayes(fit)

## Standard GO analysis

go.fisher <- goana(fit, species="Hs")
topGO(go.fisher, sort = "up")
topGO(go.fisher, sort = "down")






### Trying a more direct way  # doesnt really work due to "org...."package
### INstall and lead stuff

#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("biomaRt")

library("org.Hs.eg.db") # remember to install it if you don't have it already
symbols <- mapIds(org.Hs.eg, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")

library("biomaRt")

ensembl = useMart("ensembl",dataset=pacman)
getBM(attributes='hgnc_symbol', 
      filters = 'ensembl_gene_id', 
      values = ensemblsIDS, 
      mart = ensembl)

#library("org.Hs.eg.db") # remember to install it if you don't have it already
#symbols <- mapIds(org.Hs.eg, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
#library("biomaRt")
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#getBM(attributes='hgnc_symbol', 
#      filters = 'ensembl_gene_id', 
#      values = ensemblsIDS, 
#      mart = ensembl)
