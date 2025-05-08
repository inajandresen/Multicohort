rm(list=ls())
library(gage)
library(PADOG)
library(org.Hs.eg.db)
library(foreach)
library(doRNG)
library(tidyverse)
library(reshape2)
library(ggforce)
library(ggpubr)
library(writexl)



# Format dataset ----------------------------------------------------------
mygslist=readList("../Data/c2.cp.v2023.2.Hs.entrez.gmt")

load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
prot <- fData(Merged_cohorts_wins_Detrended_exprset)
prot$TargetFullName_EntrezGeneSymbol <- paste(prot$EntrezGeneID, prot$EntrezGeneSymbol, prot$TargetFullName, prot$SomaId)
prot2 <- prot %>%
  dplyr::select(TargetFullName_EntrezGeneSymbol, EntrezGeneID)
dup <- prot2$EntrezGeneID[duplicated(prot2$EntrezGeneID)]
dup <- unique(dup)


# Time interval 1 ---------------------------------------------------------
T2_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2")]
T2_MoM <- exprs(T2_Merged_cohorts_wins_Detrended_exprset)
T2_MoM <- as.data.frame(T2_MoM)
T2_MoM <- merge(prot2, T2_MoM, by = "row.names")
T2_MoM <- T2_MoM %>%
  rename("SomaId" = Row.names) 

#Filter on duplicates and make into long format
T2_MoM_dub <- filter(T2_MoM, grepl(paste(dup, collapse = "|"), EntrezGeneID))
T2_MoM_dub_long <- melt(T2_MoM_dub) 
T2_MoM_dub_long <- T2_MoM_dub_long %>%
 rename("SampleID" = variable) %>%
  rename("MoM" = value)

T2_samp <- pData(T2_Merged_cohorts_wins_Detrended_exprset)
T2_samp <- T2_samp %>%
  dplyr::select(Group) %>%
  rownames_to_column(var = "SampleID")

T2_MoM_dub_long <- merge(T2_samp, T2_MoM_dub_long, by = "SampleID")


#Plot boxplots of duplicates
T2_plot <- ggplot(T2_MoM_dub_long, aes(x = Group, y=MoM, fill = Group)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~TargetFullName_EntrezGeneSymbol, scales = "free", ncol = 6) +
  #stat_compare_means(method = "t.test") +
  stat_compare_means(method = "wilcox.test") +
  theme_classic()  

ggsave(filename= "T2_Boxplots_duplicates_wilcoxon.png",
       plot = T2_plot,
       device = "png",
       path = "../Enrichment/Padog/",
       width = 60,
       height = 55,
       units = "cm")

#The duplicated proteins I want to remove
#"718"    SL000313 SL000312 SL003362       
#"720 721"   SL000316   
#"727"     SL000319
#"2243 2244 2266" SL000424
#"4049 4050"  SL000508  
#"5340"    SL000540 SL000541       
#"6368"   SL003302       
#"7422"   SL000002      
#"2159"  SL003324         
#"5443"   SL000300       
#"5624"  SL000048      
#"6320"    SL004363   
#"2158" SL004400         
#"5473"   SL003191      
#"2253" SL014308     

remove <- c("SL000313", "SL000312", "SL003362", "SL000316", "SL000319",
            "SL000424", "SL000508", "SL000540", "SL000541", "SL003302",
            "SL000002", "SL003324", "SL000300", "SL000048", "SL004363",
            "SL004400", "SL003191", "SL014308")

T2_MoM <- T2_MoM[!T2_MoM$SomaId %in% remove ,]

T2_MoM <- T2_MoM %>%
  remove_rownames() %>%
  dplyr::select(-SomaId, -TargetFullName_EntrezGeneSymbol) %>%
  column_to_rownames(var = "EntrezGeneID") %>%
  as.matrix()


T2_samp$Group=(ifelse(T2_samp$Group=="PE","d","c"))
T2_samp <- T2_samp %>%
  mutate("Sample2" = SampleID) %>%
  column_to_rownames(var = "Sample2")

all(rownames(T2_samp)== colnames(T2_MoM))


#############PADOG
T2_myr=padog(
  esetm=T2_MoM,
  group=T2_samp$Group,
  paired=FALSE,
  annotation=NULL,
  gslist=mygslist,
  organism="hsa",
  verbose=TRUE,
  Nmin=5,
  NI=2000,
  plots=FALSE,
  parallel = T,
  ncr = 8)

# use Ppadog_FDR to determine significance
T2_myr$Ppadog_FDR<-p.adjust(T2_myr$Ppadog,"fdr")
T2_myr$PmeanAbsT_FDR<-p.adjust(T2_myr$PmeanAbsT,"fdr")

T2_dotplot_padog <- dotplot(T2_myr, x = "richFactor", showCategory = 20, title = "Disease ontology") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) #Adjusts the length of the labels on y-axis


write_xlsx(T2_myr, path = "../Enrichment/Padog/Padog_res_T2.xlsx",
           col_names = TRUE)

# Time interval 2 ---------------------------------------------------------
T3_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T3")]
T3_MoM <- exprs(T3_Merged_cohorts_wins_Detrended_exprset)
T3_MoM <- as.data.frame(T3_MoM)
T3_MoM <- merge(prot2, T3_MoM, by = "row.names")
T3_MoM <- T3_MoM %>%
  rename("SomaId" = Row.names) 

#Filter on duplicates and make into long format
T3_MoM_dub <- filter(T3_MoM, grepl(paste(dup, collapse = "|"), EntrezGeneID))
T3_MoM_dub_long <- melt(T3_MoM_dub) 
T3_MoM_dub_long <- T3_MoM_dub_long %>%
  rename("SampleID" = variable) %>%
  rename("MoM" = value)

T3_samp <- pData(T3_Merged_cohorts_wins_Detrended_exprset)
T3_samp <- T3_samp %>%
  dplyr::select(Group) %>%
  rownames_to_column(var = "SampleID")

T3_MoM_dub_long <- merge(T3_samp, T3_MoM_dub_long, by = "SampleID")


#Plot boxplots of duplicates
T3_plot <- ggplot(T3_MoM_dub_long, aes(x = Group, y=MoM, fill = Group)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~TargetFullName_EntrezGeneSymbol, scales = "free", ncol = 6) +
  #stat_compare_means(method = "t.test") +
  stat_compare_means(method = "wilcox.test") +
  theme_classic()  

ggsave(filename= "T3_Boxplots_duplicates_wilcoxon.png",
       plot = T3_plot,
       device = "png",
       path = "../Padog/",
       width = 60,
       height = 55,
       units = "cm")

#The duplicated proteins I want to remove
#"718"    SL000313 SL000312 SL000456       
#"720 721"   SL000316   
#"727"     SL000319
#"2243 2244 2266" SL000424
#"4049 4050"  SL000508  
#"5340"    SL000540 SL000541       
#"6368"   SL003301       
#"7422"   SL000002      
#"2159"  SL000360         
#"5443"   SL000300       
#"5624"  SL000048      
#"6320"    SL004363   
#"2158" SL000357         
#"5473"   SL003191      
#"2253" SL014308     

remove <- c("SL000313", "SL000312", "SL000456", "SL000316", "SL000319", "SL000424", "SL000508", "SL000540", "SL000541",
            "SL003301", "SL000002", "SL000360", "SL000300", "SL000048", "SL004363", "SL000357", "SL003191", "SL014308")

T3_MoM <- T3_MoM[!T3_MoM$SomaId %in% remove ,]

T3_MoM <- T3_MoM %>%
  remove_rownames() %>%
  dplyr::select(-SomaId, -TargetFullName_EntrezGeneSymbol) %>%
  column_to_rownames(var = "EntrezGeneID") %>%
  as.matrix()


T3_samp$Group=(ifelse(T3_samp$Group=="PE","d","c"))
T3_samp <- T3_samp %>%
  mutate("Sample2" = SampleID) %>%
  column_to_rownames(var = "Sample2")

all(rownames(T3_samp)== colnames(T3_MoM))


#############PADOG
T3_myr=padog(
  esetm=T3_MoM,
  group=T3_samp$Group,
  paired=FALSE,
  annotation=NULL,
  gslist=mygslist,
  organism="hsa",
  verbose=TRUE,
  Nmin=5,
  NI=2000,
  plots=FALSE,
  parallel = T,
  ncr = 8)

# use Ppadog_FDR to determine significance
T3_myr$Ppadog_FDR<-p.adjust(T3_myr$Ppadog,"fdr")
T3_myr$PmeanAbsT_FDR<-p.adjust(T3_myr$PmeanAbsT,"fdr")

write_xlsx(T3_myr, path = "../Enrichment/Padog/Padog_res_T3.xlsx",
           col_names = TRUE)

# Time interval 2 ---------------------------------------------------------
T4_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T4")]
T4_MoM <- exprs(T4_Merged_cohorts_wins_Detrended_exprset)
T4_MoM <- as.data.frame(T4_MoM)
T4_MoM <- merge(prot2, T4_MoM, by = "row.names")
T4_MoM <- T4_MoM %>%
  rename("SomaId" = Row.names) 

#Filter on duplicates and make into long format
T4_MoM_dub <- filter(T4_MoM, grepl(paste(dup, collapse = "|"), EntrezGeneID))
T4_MoM_dub_long <- melt(T4_MoM_dub) 
T4_MoM_dub_long <- T4_MoM_dub_long %>%
  rename("SampleID" = variable) %>%
  rename("MoM" = value)

T4_samp <- pData(T4_Merged_cohorts_wins_Detrended_exprset)
T4_samp <- T4_samp %>%
  dplyr::select(Group) %>%
  rownames_to_column(var = "SampleID")

T4_MoM_dub_long <- merge(T4_samp, T4_MoM_dub_long, by = "SampleID")


#Plot boxplots of duplicates
T4_plot <- ggplot(T4_MoM_dub_long, aes(x = Group, y=MoM, fill = Group)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~TargetFullName_EntrezGeneSymbol, scales = "free", ncol = 6) +
  #stat_compare_means(method = "t.test") +
  stat_compare_means(method = "wilcox.test") +
  theme_classic()  

ggsave(filename= "T4_Boxplots_duplicates_wilcoxon.png",
       plot = T4_plot,
       device = "png",
       path = "../Padog/",
       width = 60,
       height = 55,
       units = "cm")

#The duplicated proteins I want to remove
#"718"    SL000313 SL000312 SL000456       
#"720 721"   SL000316   
#"727"     SL000319
#"2243 2244 2266" SL000424
#"4049 4050"  SL000507  
#"5340"    SL000268 SL000540       
#"6368"   SL003302       
#"7422"   SL000002      
#"2159"  SL000360         
#"5443"   SL003461       
#"5624"  SL003974      
#"6320"    SL004363   
#"2158" SL000357         
#"5473"   SL004708      
#"2253" SL004342     

remove <- c("SL000313", "SL000312", "SL000456", "SL000316", "SL000319", "SL000424", "SL000507", 
            "SL000268", "SL000540", "SL003302", "SL000002", "SL000360", "SL003461", "SL003974",
            "SL004363", "SL000357", "SL004708", "SL004342")

T4_MoM <- T4_MoM[!T4_MoM$SomaId %in% remove ,]

T4_MoM <- T4_MoM %>%
  remove_rownames() %>%
  dplyr::select(-SomaId, -TargetFullName_EntrezGeneSymbol) %>%
  column_to_rownames(var = "EntrezGeneID") %>%
  as.matrix()


T4_samp$Group=(ifelse(T4_samp$Group=="PE","d","c"))
T4_samp <- T4_samp %>%
  mutate("Sample2" = SampleID) %>%
  column_to_rownames(var = "Sample2")

all(rownames(T4_samp)== colnames(T4_MoM))


#############PADOG
T4_myr=padog(
  esetm=T4_MoM,
  group=T4_samp$Group,
  paired=FALSE,
  annotation=NULL,
  gslist=mygslist,
  organism="hsa",
  verbose=TRUE,
  Nmin=5,
  NI=2000,
  plots=FALSE,
  parallel = T,
  ncr = 8)

# use Ppadog_FDR to determine significance
T4_myr$Ppadog_FDR<-p.adjust(T4_myr$Ppadog,"fdr")
T4_myr$PmeanAbsT_FDR<-p.adjust(T4_myr$PmeanAbsT,"fdr")

write_xlsx(T4_myr, path = "../Enrichment/Padog/Padog_res_T4.xlsx",
           col_names = TRUE)
