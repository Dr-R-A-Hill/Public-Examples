### Setting work directory
#setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Transcriptomic_Data/Euc-Mut_Comparison_Data')
setwd('E:/WSU Work/Transcriptomic_Data/Euc-Mut_Comparison_Data')
setwd('D:/WSU Work/Transcriptomic_Data/Euc-Mut_Comparison_Data')


### Call files
#reddeer<-read.csv("SI14_Euc-Pis_Compared_Timecourse_DGE_sig_named.csv",fill=TRUE, header = TRUE)
reddeer<-read.csv("Euc-Pis_Compared_Relitive-Timecourse_DGE-Checked_Sig0.05_WithDescriptions.csv",fill=TRUE, header = TRUE)

#yellowdeer<-read.csv("top_2517_sig_euc-SI14_genes_mut_all_rallcon.csv",fill=TRUE, header = TRUE)  ## sig genes only
#yellowdeer<-read.table("ABA Responsive Gene IDs.txt",fill=TRUE, header = TRUE)  ## sig genes only
yellowdeer<-read.csv("euc-pis_sig_pos_precontact_genes.csv",fill=TRUE, header = TRUE)  ## sig genes only

#greendeer<-as.character(yellowdeer$LocusName)

bluedeer<-subset(reddeer,LocusName %in% yellowdeer$x)

write.csv(bluedeer,"euc-pis_sig_pos_precontact_genes_AllData.csv")


#str(reddeer$LocusName)


### Mean condition groups and then add as a new column to the table

#elk<-rowMeans(cbind(reddeer$cond_1v2_log2FC,reddeer$cond_1v3_log2FC,reddeer$cond_1v4_log2FC,reddeer$cond_1v5_log2FC,reddeer$cond_1v6_log2FC),na.rm=FALSE, dims=+1)
#str(elk)
#reddeer$cond_1_mean<-elk



### Arange the two lists in order and take a subset of the 500 most differenceated

##attach(reddeer)
reddeer1<-subset(reddeer, Precontact_vs_Control_Sig == "TRUE")
reddeer1a<-subset(reddeer1, Precontact_vs_Control_log2FC >= 1 | Precontact_vs_Control_log2FC <= -1)
reddeer1a_neg<-subset(reddeer1, Precontact_vs_Control_log2FC <= -1)
reddeer1a_pos<-subset(reddeer1, Precontact_vs_Control_log2FC >= 1)
antilope1<-reddeer1a[order(reddeer1a$Precontact_vs_Control_Padj),]    ### orders the data, giving the highest first
impalas<-antilope1[0:312,"LocusName"]         # formally 2314,378   now 312  ### gets the first xxx objects 
impala<-as.character(impalas)
#impala[2314]
impala1<-antilope1[0:312,"arabi_defline"] 
impala1a<-as.character(impala1)
impala1b<-antilope1[0:312,"KO"] 
impala1c<-as.character(impala1b)
impala1d<-antilope1[0:312,"GeneName"]     ### For gene enrichment
impala1e<-as.character(impala1d)
impala1f<-antilope1[0:312,"Precontact_vs_Control_Padj"] 
impala1g<-as.character(impala1f)
#impala2<-mapply(c, impala,impala1a, SIMPLIFY = FALSE)
impala2<-data.frame("LocusName"=impala,"arabi_defline"=impala1a,"KO"=impala1c,"GeneName"=impala1e, "Treatment_vs_Control_Padj"=impala1g)#)

reddeer2<-subset(reddeer, Contact_24hr_vs_Control_Sig == "TRUE")
reddeer2a<-subset(reddeer2, Contact_24hr_vs_Control_log2FC >= 1 | Contact_24hr_vs_Control_log2FC <= -1)
reddeer2a_pos<-subset(reddeer2, Contact_24hr_vs_Control_log2FC >= 1)
reddeer2a_neg<-subset(reddeer2, Contact_24hr_vs_Control_log2FC <= -1)
antilope2<-reddeer2a[order(reddeer2a$Contact_24hr_vs_Control_Padj),]
pronghorns<-antilope2[0:3480,"LocusName"]   # formally 3074,5590 now 3480
pronghorn<-as.character(pronghorns)
#pronghorn[450]
##detach(reddear)
pronghorn1<-antilope2[0:3480,"arabi_defline"] 
pronghorn1a<-as.character(pronghorn1)
pronghorn1b<-antilope2[0:3480,"KO"] 
pronghorn1c<-as.character(pronghorn1b)
pronghorn1d<-antilope2[0:3480,"GeneName"] 
pronghorn1e<-as.character(pronghorn1d)
pronghorn1f<-antilope2[0:3480,"Contact_24hr_vs_Control_Padj"] 
pronghorn1g<-as.character(pronghorn1f)
#pronghorn2<-mapply(c, pronghorn,pronghorn1a, SIMPLIFY = FALSE)
pronghorn2<-data.frame("LocusName"=pronghorn,"arabi_defline"=pronghorn1a,"KO"=pronghorn1c,"GeneName"=pronghorn1e,"Treatment_vs_Control_Padj"=pronghorn1g)#)


reddeer3<-subset(reddeer, Contact_48hr_vs_Control_Sig == "TRUE")
reddeer3a<-subset(reddeer3, Contact_48hr_vs_Control_log2FC >= 1 | Contact_48hr_vs_Control_log2FC <= -1)
reddeer3a_pos<-subset(reddeer3, Contact_48hr_vs_Control_log2FC >= 1)
reddeer3a_neg<-subset(reddeer3, Contact_48hr_vs_Control_log2FC <= -1)
antilope3<-reddeer3a[order(reddeer3a$Contact_48hr_vs_Control_Padj),]
prongs<-antilope3[0:3300,"LocusName"]   # formally 2844,5267 now 3300
prong<-as.character(prongs) 
#prong
prong1<-antilope3[0:3300,"arabi_defline"] 
prong1a<-as.character(prong1)
prong1b<-antilope3[0:3300,"KO"] 
prong1c<-as.character(prong1b)
prong1d<-antilope3[0:3300,"GeneName"] 
prong1e<-as.character(prong1d)
prong1f<-antilope3[0:3300,"Contact_48hr_vs_Control_Padj"] 
prong1g<-as.character(prong1f)
#prong2<-mapply(c, prong,prong1a, SIMPLIFY = FALSE)
prong2<-data.frame("LocusName"=prong,"arabi_defline"=prong1a,"KO"=prong1c,"GeneName"=prong1e,"Treatment_vs_Control_Padj"=prong1g)


reddeer4<-subset(reddeer, Contact_1week_vs_Control_Sig == "TRUE")
reddeer4a<-subset(reddeer4, Contact_1week_vs_Control_log2FC >= 1 | Contact_1week_vs_Control_log2FC <= -1)
reddeer4a_pos<-subset(reddeer4, Contact_1week_vs_Control_log2FC >= 1)
reddeer4a_neg<-subset(reddeer4,  Contact_1week_vs_Control_log2FC <= -1)
antilope4<-reddeer4a[order(reddeer4a$Contact_1week_vs_Control_Padj),]
horns<-antilope4[0:3111,"LocusName"]   # Formally 467,5902 now 3111
horn<-as.character(horns) 
#horn[450]
horn1<-antilope4[0:3111,"arabi_defline"] 
horn1a<-as.character(horn1)
horn1b<-antilope4[0:3111,"KO"] 
horn1c<-as.character(horn1b)
horn1d<-antilope4[0:3111,"GeneName"] 
horn1e<-as.character(horn1d)
horn1f<-antilope4[0:3111,"Contact_1week_vs_Control_Padj"] 
horn1g<-as.character(horn1f)
#horn2<-mapply(c, horn,horn1a, SIMPLIFY = FALSE)
horn2<-data.frame("LocusName"=horn,"arabi_defline"=horn1a,"KO"=horn1c,"GeneName"= horn1e,"Treatment_vs_Control_Padj"=horn1g)


reddeer5<-subset(reddeer,Contact_2week_vs_Control_Sig == "TRUE")
reddeer5a<-subset(reddeer5, Contact_2week_vs_Control_log2FC >= 1 | Contact_2week_vs_Control_log2FC <= -1)
reddeer5a_pos<-subset(reddeer5, Contact_2week_vs_Control_log2FC >= 1)
reddeer5a_neg<-subset(reddeer5, Contact_2week_vs_Control_log2FC <= -1)

antilope5<-reddeer5a[order(reddeer5a$Contact_2week_vs_Control_Padj),]
horses<-antilope5[0:922,"LocusName"]   # Formally 3073 1217     now 922            ###horse is a factor that gives me the row number when I try to use it
horse<-as.character(horses)                           ### gets the titles only and puts them in a seporate object
#horses
horse1<-antilope5[0:922,"arabi_defline"] 
horse1a<-as.character(horse1)
horse1b<-antilope5[0:922,"KO"] 
horse1c<-as.character(horse1b)
horse1d<-antilope5[0:922,"GeneName"] 
horse1e<-as.character(horse1d)
horse1f<-antilope5[0:922,"Contact_2week_vs_Control_Padj"] 
horse1g<-as.character(horse1f)
#horse2<-mapply(c,horse,horse1a, SIMPLIFY = FALSE)
horse2<-data.frame("LocusName"=horse,"arabi_defline"=horse1a,"KO"=horse1c,"GeneName"=horse1e,"Treatment_vs_Control_Padj"=horse1g)
##horse1b<-horse1a[""=="n^1"]
#horse1a[horse1a==""] <- 4
#horse1a[1:10]



carotenoid_genes<-c("Eucgr.H04384.v2.0", "Eucgr.L02103.v2.0","Eucgr.F02914.v2.0","Eucgr.C03415.v2.0","Eucgr.A02546.v2.0","Eucgr.A01599.v2.0","Eucgr.I00923.v2.0","Eucgr.B02180.v2.0","Eucgr.D00291.v2.0","Eucgr.D02555.v2.0","Eucgr.F01409.v2.0", "Eucgr.F03199.v2.0","Eucgr.G03255.v2.0","Eucgr.D01823.v2.0", "Eucgr.K00294.v2.0")
precontact_carotenoid_genes<-intersect(impala,carotenoid_genes)
x24hr_carotenoid_genes<-intersect(pronghorn,carotenoid_genes)
x48hr_carotenoid_genes<-intersect(prong,carotenoid_genes)
x1week_carotenoid_genes<-intersect(horn,carotenoid_genes)
x2weeks_carotenoid_genes<-intersect(horse,carotenoid_genes)




#write.csv(reddeer1a_pos$LocusName,"euc-pis_sig_pos_precontact_genes.csv")
#write.csv(reddeer2a_pos$LocusName,"euc-pis_sig_pos_24hr_genes.csv")
#write.csv(reddeer3a_pos$LocusName,"euc-pis_sig_pos_48hr_genes.csv")
#write.csv(reddeer4a_pos$LocusName,"euc-pis_sig_pos_1week_genes.csv")
#write.csv(reddeer5a_pos$LocusName,"euc-pis_sig_pos_2week_genes.csv")


#write.csv(reddeer1a_neg$LocusName,"euc-pis_sig_neg_precontact_genes.csv")
#write.csv(reddeer2a_neg$LocusName,"euc-pis_sig_neg_24hr_genes.csv")
#write.csv(reddeer3a_neg$LocusName,"euc-pis_sig_neg_48hr_genes.csv")
#write.csv(reddeer4a_neg$LocusName,"euc-pis_sig_neg_1week_genes.csv")
#write.csv(reddeer5a_neg$LocusName,"euc-pis_sig_neg_2week_genes.csv")


### combine all the lists of sig genes by LocusName

herd1<-union(impala,pronghorn)     ### Joins 2 lists with no repeats
herd2<-union(herd1,prong)
herd3<-union(herd2,horn)
herd4<-union(herd3,horse)
herds<-subset(reddeer,LocusName %in% herd4)
herding<-data.frame("LocusName"=herds$LocusName,"arabi_defline"=herds$arabi_defline,"KO"=herds$KO,"GeneName"=herds$GeneName)

herd1_pos<-union(reddeer1a_pos$LocusName,reddeer2a_pos$LocusName)     ### Joins 2 lists with no repeats
herd2_pos<-union(herd1_pos,reddeer3a_pos$LocusName)
herd3_pos<-union(herd2_pos,reddeer4a_pos$LocusName)
herd4_pos<-union(herd3_pos,reddeer5a_pos$LocusName)
herds_pos<-subset(reddeer,LocusName %in% herd4_pos)
herding_pos<-data.frame("LocusName"=herds_pos$LocusName,"arabi_defline"=herds_pos$arabi_defline,"KO"=herds_pos$KO,"GeneName"=herds_pos$GeneName)


write.csv(herds_pos,"3435_all-pos_sig_euc-SI14_genes_mut_all_rallcon.csv")


herd1_neg<-union(reddeer1a_neg$LocusName,reddeer2a_neg$LocusName)     ### Joins 2 lists with no repeats
herd2_neg<-union(herd1_neg,reddeer3a_neg$LocusName)
herd3_neg<-union(herd2_neg,reddeer4a_neg$LocusName)
herd4_neg<-union(herd3_neg,reddeer5a_neg$LocusName)
herds_neg<-subset(reddeer,LocusName %in% herd4_neg)
herding_neg<-data.frame("LocusName"=herds_neg$LocusName,"arabi_defline"=herds_neg$arabi_defline,"KO"=herds_neg$KO,"GeneName"=herds_neg$GeneName)


write.csv(herds_neg,"3372_all-neg_sig_euc-SI14_genes_mut_all_rallcon.csv")



### Find the sig genes at one timepoint that are in common with all the previously sig expressed genes in the lists
### basically, we've seen these genes before, this is what they are doing now

sharedherd1<-intersect(impala,pronghorn)     
sharedherd2<-intersect(herd1,prong)   # all sig 48hr genes in lists found in all the genes common to precontact and 24hr
sharedherd3<-intersect(herd2,horn)
sharedherd4<-intersect(herd3,horse)



### combine all the lists of sig genes by LocusName - only including those appearing in all lists - be careful of object names
### Find the persistant genes and orginise the data 

iherd1<-intersect(impala,pronghorn)     
iherd2<-intersect(iherd1,prong)
iherd3<-intersect(iherd2,horn)
iherd4<-intersect(iherd3,horse)

#intersected_gene_info<-subset(reddeer,LocusName %in% iherd4)
#write.csv(intersected_gene_info,"109-intersected_gene_info_euc-SI14_genes_mut_persistent_rallcon_AllData.csv")


gooddeer<-intersect(herd4,yellowdeer$LocusName)
## gooddeer<-subset(olddeera,LocusName %in% newdeer$LocusName)

### Sourcing ggene lists from elsewhere
### Setting work directory
#setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Metabolism Driven Data Analysis Pipeline/Euc-Pis/All Significant Genes')

### Call files
#meta_genes<-read.csv("all-6687_sig_euc-SI14_genes_mut_callcon_rallconditions_KEGG_50KO_Photosynthisis-Genes.csv",fill=TRUE, header = TRUE)
#meta_gene_names<-as.character(meta_genes$LocusName) # might not be neccesary, 4 genes missing due to repeats
#meta_gene_info<-subset(meta_genes,LocusName %in% meta_gene_names)



### Get the annotations for our selected genes  --- You only need to change the list and the output file name at this point

iherd_detail<-subset(reddeer,LocusName %in% iherd4) # all the information regarding persistant genes
#iherd_detail<-subset(reddeer,LocusName %in% impala2$LocusName)
#iherd_detail<-subset(reddeer,LocusName %in% meta_genes$LocusName)


#iherd_point<-data.frame(iherd_detail$LocusName,iherd_detail$arabi_defline,iherd_detail$KEGG_ec,iherd_detail$KO,iherd_detail$GO) # the abridged version
iherd_point<-data.frame(iherd_detail$LocusName,iherd_detail$arabi_defline,iherd_detail$KEGG_ec,iherd_detail$KO,iherd_detail$GO,iherd_detail$GeneName)

#iherd_point_plus1<-rbind(iherd_point,precon_unique_point[1:5,])
#iherd_point_plus2<-rbind(iherd_point_plus1,X24hr_unique_point[1:5,])
#iherd_point_plus3<-rbind(iherd_point_plus2,X48hr_unique_point[1:5,])
#iherd_point_plus4<-rbind(iherd_point_plus3,X1week_unique_point[1:5,])
#iherd_point_plus5<-rbind(iherd_point_plus4,X2week_unique_point[1:5,])


#peasants<-Reduce(intersect, list(impala,pronghorn,prong,horn,horse))

#commondata<-subset(reddeer, LocusName %in% iherd4)   ### change iherd number to investigate a different length of persistance
#commondata<-subset(reddeer, LocusName %in% sharedherd4)
#filteredcommondata<-data.frame(as.character(commondata$LocusName),as.character(commondata$Precontact_vs_Control_log2FC),as.character(commondata$Contact_24hr_vs_Control_log2FC),as.character(commondata$Contact_48hr_vs_Control_log2FC),as.character(commondata$Week_1_Contact_vs_Control_log2FC),as.character(commondata$Week_2_Contact_vs_Control_log2FC),as.character(commondata$LocusName),as.character(commondata$arabi_defline),as.character(commondata$KO))
#filteredcommondata<-data.frame(commondata$LocusName,commondata$Precontact_vs_Control_log2FC,commondata$Contact_24hr_vs_Control_log2FC,commondata$Contact_48hr_vs_Control_log2FC,commondata$Week_1_Contact_vs_Control_log2FC,commondata$Week_2_Contact_vs_Control_log2FC,commondata$LocusName,commondata$arabi_defline)
#filteredcommondata<-cbind(as.character(commondata$LocusName),as.character(commondata$Precontact_vs_Control_log2FC))             
#setnames(filteredcommondata, old=c("old_name","another_old_name"), new=c(
#names(filteredcommondata) <- c("LocusName","Precontact_vs_Control_log2FC","Contact_24hr_vs_Control_log2FC","Contact_48hr_vs_Control_log2FC","Week_1_Contact_vs_Control_log2FC","Week_2_Contact_vs_Control_log2FC","LocusName","arabi_defline","KO")

#fcd_key<-data.frame(as.character(filteredcommondata$LocusName),as.character(filteredcommondata$arabi_defline),as.character(filteredcommondata$KO))

#herd1a<-union(impala2,pronghorn2)

### combine all the lists of sig genes by arabi_defline - seems to take each row as a single entry?cbind
#herd1a<-rbind(impala2,pronghorn2)     ### Joins 2 lists with repeats
#herd2aa<-rbind(herd1a,prong2)
#herd2a<-unique(herd2aa) 
#herd3aa<-rbind(herd2a,horn2)
#herd3a<-unique(herd3aa) 
#herd4aa<-rbind(herd3a,horse2)
#herd4a<-unique(herd4aa)                ### Gets rid of repeat rows



#herda_names<-as.character(herd4a$LocusName)
#herda_descriptions<-as.character(herd4a$arabi_defline)
#herda_uni<-data.frame("LocusName"=herda_names,"arabi_defline"=herda_descriptions)    ### This is nessesary because the table will otherwise treat the data from different tables seporate 


#usaid<-subset(herda_uni, LocusName %in% iherd4) ### Creates a subset of the full significant gene list and extracts only those genes appearing in all conditions

###
####### Extracting the raw data
###
### Setting work directory

#setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Transcriptomic_Data/Euc-Mut_Raw_Data')
setwd('E:/WSU Work/Transcriptomic_Data/Euc-Mut_Raw_Data')

### Call files
hornbill<-read.csv("SI14_Timecourse_Euc_FPKM_JMPcalc_07092019.csv",fill=TRUE, header = TRUE)#xlsx


###   Extract the rows of data that are the same as those from the significant combined datasets
#parrot<-subset(hornbill,GeneID %in% (filteredcommondata$LocusName))
#parrot<-subset(hornbill,GeneID %in% (iherd_point$iherd_detail.LocusName))
#parrot<-subset(hornbill,GeneID %in% (herd4))  ### precontact only

parrot<-subset(hornbill,GeneID %in% (herds$LocusName))

### Sort the data


sorted_parrot <- with(parrot, parrot[order(parrot$GeneID) , ])


#sorted_names <- herding[order(herding$LocusName),] 
#sorted_names <- fcd_key[order(fcd_key$as.character.filteredcommondata.LocusName),] 
sorted_names <- iherd_point[order(iherd_point$iherd_detail.LocusName),] 


### Add the gene descriptions to the significant raw data 


puffin<-cbind(sorted_parrot,sorted_names)


### Write to file

#setwd('\\\\ad.uws.edu.au/dfshare/HomesBLK$/90940729/Desktop/Transcriptomic_Data/Euc-Mut_Raw_Data/Processed_Data')

#write.csv(bluedeer,"Carotenoid_Pathway_Genes_v2_AllData.csv")

setwd('E:/WSU Work/Transcriptomic_Data/Euc-Mut_Raw_Data/Processed_Data')

write.csv(parrot,"all_6687_sig_euc-SI14_FPKM_mut_all_rallcon.csv")


