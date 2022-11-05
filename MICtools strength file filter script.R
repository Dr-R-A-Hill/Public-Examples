### Setting work directory
setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes')
setwd('//ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Pisolithus - Eucalyptus Networks/6687 euc genes and 141 pis exudes')


### Call files
rodent<-read.table("all_141-Pis-Exudes_and_all_6687_sig_FPKM_ForMICtools_strength.txt",fill=TRUE, header = TRUE)


###Filter
rodent_f1<-subset(rodent, MICe >= 0.8)
rodent_f1$PosNeg<- ifelse(rodent_f1$SpearmanRho >=0, "Pos", "Neg")   ### only do once

#write.csv(rodent_f1,"all_141-Pis-Exudes_and_all_6687_sig_strength_FILTERED.csv")

### STAGE 2

### Call files
rodent_f1<-read.csv("all_141-Pis-Exudes_and_all_6687_sig_strength_FILTERED.csv",fill=TRUE, header = TRUE)





#rodent_f2 <- factor(c("Prot"))
#rodent_f2<-rodent_f1[rodent_f1$Var1 %like% "P"]

#Months[Name %like% "mb"]
#x <- x[grepl("Aisle", x[["column1"]])


#rodent_f2 <- rodent_f1[grepl("Prot", rodent_f1[["Var1"]])]

##rodent_f2<-subset(rodent_f1, Var1 %in% "Prot")

###add genoral ID columns
rodent_f2<-substr(rodent_f1$Var1, 1,2)
rodent_f1$var1_id<- rodent_f2
rodent_f2<-substr(rodent_f1$Var2, 1,2)
rodent_f1$var2_id<- rodent_f2
write.csv(rodent_f1,"all_141-Pis-Exudes_and_all_6687_sig_strength_FILTERED_IDs.csv")


rodent_f2_1<-subset(rodent_f1, var1_id == "Pr")
rodent_f2_2<-subset(rodent_f1, var2_id == "Pr")

rodent_f3<-rbind(rodent_f2_1,rodent_f2_2)
write.csv(rodent_f3,"all_141-Pis-Exudes_and_all_6687_sig_strength_FILTERED_IDs_PisInteractionsOnly.csv")


###Stage 3


### Setting work directory
setwd('E:/WSU Work/Transcriptomic_Data/Euc-Mut_Comparison_Data')
### Call files
euc_genes<-read.csv("6687_all_sig_euc-SI14_genes_mut_all_rallcon.csv",fill=TRUE, header = TRUE)

### Setting work directory
setwd('E:/WSU Work/Transcriptomic_Data/Fungi-PIS-Mut_Comparison_Data-Revisited')
### Call files
pis_genes<-read.csv("Pisolithus-Genes_Euc-Pis_Compared_Relitive-Timecourse_DGE-WithDescriptions.csv",fill=TRUE, header = TRUE)


### Setting work directory
setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes')
setwd('//ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Pisolithus - Eucalyptus Networks/6687 euc genes and 141 pis exudes')


### Call files
canis<-read.csv("all_141-Pis-Exudes_and_all_6687_sig_strength_FILTERED_IDs_PisInteractionsOnly.csv",fill=TRUE, header = TRUE)

canis_v1<-as.character(canis$Var1)
canis_v2<-as.character(canis$Var2)
canis_canis<-union(canis_v1,canis_v2)   ## all genes/SSPs - no repeats - 796


canis_euc<-subset(euc_genes,LocusName %in% canis_canis)
canis_pis<-subset(pis_genes,ProtID_Labelled %in% canis_canis)

setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes/Descriptions')
write.csv(canis_euc,"all_141-Pis-Exudes_and_all_6687___655_Euc_Descriptions.csv")
write.csv(canis_pis,"all_141-Pis-Exudes_and_all_6687___141_Pis_Descriptions.csv")




###Stage 4

setwd('C:/Users/Richard/OneDrive - Western Sydney University/Auxin Paper/Auxin Genes')
auxin_genes<-read.csv("Selected All Auxin Enzyme Genes.csv",fill=TRUE, header = TRUE)

setwd('C:/Users/Richard/OneDrive - Western Sydney University/Carotenoid and ABA Paper/Carotenoid Heatmaps')
carotenoid_genes<-read.csv("15_KEGG_significant_carotenoid_genes_only.csv",fill=TRUE, header = TRUE)


setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes/Descriptions')
canis_desc<-read.csv("Euc-Pis_PisOnlyInteractions_Network_Descriptions.csv",fill=TRUE, header = TRUE)

canis_desc_auxin<-subset(auxin_genes, LocusName %in% canis_desc$Name1)
canis_desc_caro<-subset(carotenoid_genes, LocusName %in% canis_desc$Name1)

setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes/Descriptions')
write.csv(canis_desc_auxin,"all_141-Pis-Exudes_and_all_6687___1_auxin_gene.csv")
write.csv(canis_desc_caro,"all_141-Pis-Exudes_and_all_6687___4_carotenoid_genes.csv")


