### Stage 1
###Filter strength file
### Setting work directory
#setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes')
#setwd('//ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Pisolithus - Eucalyptus Networks/6687 euc genes and 141 pis exudes')
## Call files
#rodent<-read.table("all_141-Pis-Exudes_and_all_6687_sig_FPKM_ForMICtools_strength.txt",fill=TRUE, header = TRUE)


###Filter
#rodent_f1<-subset(rodent, MICe >= 0.8)
#rodent_f1$PosNeg<- ifelse(rodent_f1$SpearmanRho >=0, "Pos", "Neg")   ### only do once

#strength_f1_prot1<-filter(strength_f1, Var1 %in% prot)
#strength_f1_prot2<-filter(strength_f1, Var2 %in% prot)

#strength_f1$Var1_Prot<- strength_f1$Var1["P"]
#strength_f1$Var2_Prot<- substr(strength_f1, Var2, 1, 4)
#strength_f1_prot1<-subset(strength_f1, Var1_Prot == "Prot")
#strength_f1_prot2<-subset(strength_f1, Var2_Prot == "Prot")


## Stage 1 - call biosynthesis files

setwd('E:/WSU Work/Networks/Auxin-ABA Networks')

biosyn<-read.csv("IAA-ABA Biosynthesis Genes Only.csv",fill=TRUE, header = TRUE)

indu<-read.csv("IAA-ABA Biosynthesis and Induced Genes.csv",fill=TRUE, header = TRUE)

induonly<-read.csv("IAA-ABA Induced Genes Only.csv",fill=TRUE, header = TRUE)

resp<-read.csv("IAA-ABA Biosynthesis and Additional Responsive Genes.csv",fill=TRUE, header = TRUE)


### stage 2 - call strength file


setwd('E:/WSU Work/Networks/Pisolithus Networks/6687 euc genes and 141 pis exudes')

mice<-read.csv("all_141-Pis-Exudes_and_all_6687_sig_strength_FILTERED_IDs.csv",fill=TRUE, header = TRUE)


### stage 3 compare and extract


biosyn_mice_var1<-subset(mice,Var1 %in% biosyn$LocusName)
biosyn_mice_var2<-subset(mice,Var2 %in% biosyn$LocusName)



indu_mice_var1<-subset(mice,Var1 %in% indu$LocusName)
indu_mice_var2<-subset(mice,Var2 %in% indu$LocusName)


induonly_mice_var1<-subset(mice,Var1 %in% induonly$LocusName)
induonly_mice_var2<-subset(mice,Var2 %in% induonly$LocusName)

resp_mice_var1<-subset(mice,Var1 %in% resp$LocusName)
resp_mice_var2<-subset(mice,Var2 %in% resp$LocusName)



### stage 4 bind columns

biosyn_mice<-rbind(biosyn_mice_var1,biosyn_mice_var2)

indu_mice<-rbind(indu_mice_var1,indu_mice_var2)


induonly_mice<-rbind(induonly_mice_var1,induonly_mice_var2)

resp_mice<-rbind(resp_mice_var1,resp_mice_var2)


### stage 5 write to file 

setwd('E:/WSU Work/Networks/Auxin-ABA Networks')


write.csv(biosyn_mice,"IAA-ABA Biosynthesis Genes Only - strength.csv")

write.csv(indu_mice,"IAA-ABA Biosynthesis and Induced Genes - strength.csv")


write.csv(induonly_mice,"IAA-ABA Induced Genes Only - strength.csv")

write.csv(resp_mice,"IAA-ABA Biosynthesis and Additional Responsive Genes - strength.csv")



###Stage 6 - timecourse extraction

setwd('E:/WSU Work/Transcriptomic_Data/Euc-Mut_Comparison_Data')


reddeer<-read.csv("Euc-Pis_Compared_Relitive-Timecourse_DGE-Checked_Sig0.05_WithDescriptions.csv",fill=TRUE, header = TRUE)

biosyn_timecourse<-subset(reddeer,LocusName %in% biosyn$LocusName)

indu_timecourse<-subset(reddeer,LocusName %in% indu$LocusName)

induonly_timecourse<-subset(reddeer,LocusName %in% induonly$LocusName)

resp_timecourse<-subset(reddeer,LocusName %in% resp$LocusName)


setwd('E:/WSU Work/Networks/Auxin-ABA Networks')


write.csv(biosyn_timecourse,"IAA-ABA Biosynthesis Genes Only - timecourse.csv")

write.csv(indu_timecourse,"IAA-ABA Biosynthesis and Induced Genes - timecourse.csv")

write.csv(induonly_timecourse,"IAA-ABA Induced Genes Only - timecourse.csv")

write.csv(resp_timecourse,"IAA-ABA Biosynthesis and Additional Responsive Genes - timecourse.csv")
