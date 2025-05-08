library(org.Hs.eg.db)
library(readxl)
library(tidyverse)
library(writexl)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(DOSE)

# Load the significantly differential proteins ----------------------------
SigProteins_T2 <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
SigProteins_T2 <- SigProteins_T2 %>%
  filter(adj.P.Val < 0.1)

SigProteins_T3 <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
SigProteins_T3 <- SigProteins_T3 %>%
  filter(adj.P.Val < 0.1)

SigProteins_T4 <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")
SigProteins_T4 <- SigProteins_T4 %>%
  filter(adj.P.Val < 0.1)

#Load the background matrix of the 1006 proteins
background <- read_xlsx("../Enrichment/ORA/Background_list_1006.xlsx",
                        col_names = TRUE)


# ORA of Biological processes, time interval 2 --------

#Biological process
T3_enrichgo_BP <- enrichGO(SigProteins_T3$EntrezGeneID, OrgDb = "org.Hs.eg.db", ont = "BP", universe = background$EntrezGeneID, 
                           pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

#GeneRatio: the ratio of input genes that are annotated in a term
#BgRatio: the ratio of all genes that are annotated in this term
#RichFactor: the ratio of input genes (DEGs) that are annotated in a term to all genes that are annotated in this term in the backgrond
#Calculate RichFactor:
T3_enrichgo_BP <- mutate(T3_enrichgo_BP, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 

#Extract the results 
T3_enrichgo_BP_res <- T3_enrichgo_BP@result
write_xlsx(T3_enrichgo_BP_res, path = "../ORA/T3_BP_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T3_enrichgo_BP_res_sig <- T3_enrichgo_BP_res %>%
  filter(p.adjust < 0.2)
dim(T3_enrichgo_BP_res_sig)
#2 sig terms (<0.1)
#27 sig terms (<0.2)

#Map geneID to gene Symbol
T3_enrichgo_BP_clust <- pairwise_termsim(T3_enrichgo_BP)
T3_treeplot_BP <- treeplot(T3_enrichgo_BP_clust, showCategory = 27)

ggsave(filename= "T3_Biological_process_treeplot.png",
       plot = T3_treeplot_BP,
       device = "png",
       path = "../ORA/",
       width = 37,
       height = 20,
       units = "cm")


#Simplify
T3_enrichgo_BP_simp <- simplify(T3_enrichgo_BP)
T3_enrichgo_BP_res_simp <- T3_enrichgo_BP_simp@result
write_xlsx(T3_enrichgo_BP_res_simp, path = "../ORA/T3_BP_res_simp_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T3_enrichgo_BP_res_simp_sig <- T3_enrichgo_BP_res_simp %>%
  filter(p.adjust < 0.2)
dim(T3_enrichgo_BP_res_simp_sig)
#2 sig terms (<0.1)
#18 sig terms (<0.2)


#Map geneID to gene Symbol
T3_enrichgo_BP_clust_simp <- pairwise_termsim(T3_enrichgo_BP_simp)
T3_treeplot_BP_simp <- treeplot(T3_enrichgo_BP_clust_simp, showCategory = 18)

ggsave(filename= "T3_Biological_process_treeplot_simp.png",
       plot = T3_treeplot_BP_simp,
       device = "png",
       path = "../ORA/",
       width = 35,
       height = 15,
       units = "cm")

# ORA of Biological processes, time interval 3 --------

#Biological process
T4_enrichgo_BP <- enrichGO(SigProteins_T4$EntrezGeneID, OrgDb = "org.Hs.eg.db", ont = "BP", universe = background$EntrezGeneID, 
                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

#GeneRatio: the ratio of input genes that are annotated in a term
#BgRatio: the ratio of all genes that are annotated in this term
#RichFactor: the ratio of input genes (DEGs) that are annotated in a term to all genes that are annotated in this term in the backgrond
#Calculate RichFactor:
T4_enrichgo_BP <- mutate(T4_enrichgo_BP, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 

#Extract the results 
T4_enrichgo_BP_res <- T4_enrichgo_BP@result
write_xlsx(T4_enrichgo_BP_res, path = "../ORA/T4_BP_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T4_enrichgo_BP_res_sig <- T4_enrichgo_BP_res %>%
  filter(p.adjust < 0.1)
dim(T4_enrichgo_BP_res_sig)
#0 sig terms (<0.1)
#29 sig terms (<0.2)

#Map geneID to gene Symbol
T4_enrichgo_BP_clust <- pairwise_termsim(T4_enrichgo_BP)
T4_treeplot_BP <- treeplot(T4_enrichgo_BP_clust, showCategory = 29)

ggsave(filename= "T4_Biological_process_treeplot.png",
       plot = T4_treeplot_BP,
       device = "png",
       path = "../ORA/",
       width = 37,
       height = 20,
       units = "cm")


#Simplify
T4_enrichgo_BP_simp <- simplify(T4_enrichgo_BP)
T4_enrichgo_BP_res_simp <- T4_enrichgo_BP_simp@result
write_xlsx(T4_enrichgo_BP_res_simp, path = "../ORA/T4_BP_res_simp_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T4_enrichgo_BP_res_simp_sig <- T4_enrichgo_BP_res_simp %>%
  filter(p.adjust < 0.1)
dim(T4_enrichgo_BP_res_simp_sig)
#0 sig terms (<0.1)
#26 sig terms (<0.2)


#Map geneID to gene Symbol
T4_enrichgo_BP_clust_simp <- pairwise_termsim(T4_enrichgo_BP_simp)
T4_treeplot_BP_simp <- treeplot(T4_enrichgo_BP_clust_simp, showCategory = 26)

ggsave(filename= "T4_Biological_process_treeplot_simp.png",
       plot = T4_treeplot_BP_simp,
       device = "png",
       path = "../ORA/",
       width = 35,
       height = 20,
       units = "cm")


# Compare time biological processes interval 2 and 3 -------------------------------------------
T3_sig_proteins <- data.frame("Time" = "Time interval 2", "EntrezGeneID" = SigProteins_T3$EntrezGeneID)
T4_sig_proteins <- data.frame("Time" = "Time interval 3", "EntrezGeneID" = SigProteins_T4$EntrezGeneID)

T3T4_sig_proteins <- rbind(T3_sig_proteins, T4_sig_proteins)

T3T4_compare_BP <- compareCluster(EntrezGeneID~Time, data=T3T4_sig_proteins, 
                                  fun="enrichGO", ont = "BP", 
                                  OrgDb = "org.Hs.eg.db", universe = background$EntrezGeneID,
                                  pvalueCutoff = 1, qvalueCutoff = 0.15, readable = TRUE)

Compare_BP <- dotplot(T3T4_compare_BP, x = "Time", 
                      showCategory = 1000,
                      color = "qvalue") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + #Adjusts the length of the labels on y-axis)
  xlab(" ")

ggsave(filename= "CompareT2T3_BP_dotplot_q0.15.png",
       plot = Compare_BP,
       device = "png",
       path = "../ORA/",
       width = 25,
       height = 30,
       units = "cm")

#Simplyfy
T3T4_compare_BP_simp <- simplify(T3T4_compare_BP)
Compare_BP_simp <- dotplot(T3T4_compare_BP_simp, x = "Time", 
                           showCategory = 1000,
                           color = "qvalue") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + #Adjusts the length of the labels on y-axis)
  xlab(" ")

ggsave(filename= "CompareT2T3_BP_dotplot_q0.15_simp.png",
       plot = Compare_BP_simp,
       device = "png",
       path = "../ORA/",
       width = 25,
       height = 25,
       units = "cm")


# ORA of Mulecular functions, time interval 2 --------

#Biological process
T3_enrichgo_MF <- enrichGO(SigProteins_T3$EntrezGeneID, OrgDb = "org.Hs.eg.db", ont = "MF", universe = background$EntrezGeneID, 
                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

#GeneRatio: the ratio of input genes that are annotated in a term
#BgRatio: the ratio of all genes that are annotated in this term
#RichFactor: the ratio of input genes (DEGs) that are annotated in a term to all genes that are annotated in this term in the backgrond
#Calculate RichFactor:
T3_enrichgo_MF <- mutate(T3_enrichgo_MF, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 

#Extract the results 
T3_enrichgo_MF_res <- T3_enrichgo_MF@result
write_xlsx(T3_enrichgo_MF_res, path = "../ORA/T3_MF_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T3_enrichgo_MF_res_sig <- T3_enrichgo_MF_res %>%
  filter(p.adjust < 0.1)
dim(T3_enrichgo_MF_res_sig)
#6 sig terms (<0.1)
#37 sig terms (<0.2)

#Map geneID to gene Symbol
T3_enrichgo_MF_clust <- pairwise_termsim(T3_enrichgo_MF)
T3_treeplot_MF <- treeplot(T3_enrichgo_MF_clust, showCategory = 37)

ggsave(filename= "T3_Molecular_functions_treeplot.png",
       plot = T3_treeplot_MF,
       device = "png",
       path = "../ORA/",
       width = 37,
       height = 30,
       units = "cm")


#Simplify
T3_enrichgo_MF_simp <- simplify(T3_enrichgo_MF)
T3_enrichgo_MF_res_simp <- T3_enrichgo_MF_simp@result
write_xlsx(T3_enrichgo_MF_res_simp, path = "../ORA/T3_MF_res_simp_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T3_enrichgo_MF_res_simp_sig <- T3_enrichgo_MF_res_simp %>%
  filter(p.adjust < 0.1)
dim(T3_enrichgo_MF_res_simp_sig)
#5 sig terms (<0.1)
#22 sig terms (<0.2)


#Map geneID to gene Symbol
T3_enrichgo_MF_clust_simp <- pairwise_termsim(T3_enrichgo_MF_simp)
T3_treeplot_MF_simp <- treeplot(T3_enrichgo_MF_clust_simp, showCategory = 22)

ggsave(filename= "T3_Molecular_functions_treeplot_simp.png",
       plot = T3_treeplot_MF_simp,
       device = "png",
       path = "../ORA/",
       width = 35,
       height = 15,
       units = "cm")

# ORA of Mulecular functions, time interval 3 --------

T4_enrichgo_MF <- enrichGO(SigProteins_T4$EntrezGeneID, OrgDb = "org.Hs.eg.db", ont = "MF", universe = background$EntrezGeneID, 
                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

#GeneRatio: the ratio of input genes that are annotated in a term
#BgRatio: the ratio of all genes that are annotated in this term
#RichFactor: the ratio of input genes (DEGs) that are annotated in a term to all genes that are annotated in this term in the backgrond
#Calculate RichFactor:
T4_enrichgo_MF <- mutate(T4_enrichgo_MF, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 

#Extract the results 
T4_enrichgo_MF_res <- T4_enrichgo_MF@result
write_xlsx(T4_enrichgo_MF_res, path = "../ORA/T4_MF_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T4_enrichgo_MF_res_sig <- T4_enrichgo_MF_res %>%
  filter(p.adjust < 0.2)
dim(T4_enrichgo_MF_res_sig)
#24 sig terms (<0.1)
#24 sig terms (<0.15)
#27 sig terms (<0.2)

#Map geneID to gene Symbol
T4_enrichgo_MF_clust <- pairwise_termsim(T4_enrichgo_MF)
T4_treeplot_MF <- treeplot(T4_enrichgo_MF_clust, showCategory = 27)

ggsave(filename= "T4_Molecular_functions_treeplot.png",
       plot = T4_treeplot_MF,
       device = "png",
       path = "../ORA/",
       width = 37,
       height = 25,
       units = "cm")


#Simplify
T4_enrichgo_MF_simp <- simplify(T4_enrichgo_MF)
T4_enrichgo_MF_res_simp <- T4_enrichgo_MF_simp@result
write_xlsx(T4_enrichgo_MF_res_simp, path = "../ORA/T4_MF_res_simp_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T4_enrichgo_MF_res_simp_sig <- T4_enrichgo_MF_res_simp %>%
  filter(p.adjust < 0.2)
dim(T4_enrichgo_MF_res_simp_sig)
#13 sig terms (<0.1)
#13 sig terms (<0.15)
#14 sig terms (<0.2)


#Map geneID to gene Symbol
T4_enrichgo_MF_clust_simp <- pairwise_termsim(T4_enrichgo_MF_simp)
T4_treeplot_MF_simp <- treeplot(T4_enrichgo_MF_clust_simp, showCategory = 14)

ggsave(filename= "T4_Molecular_functions_treeplot_simp.png",
       plot = T4_treeplot_MF_simp,
       device = "png",
       path = "../ORA/",
       width = 35,
       height = 15,
       units = "cm")
# Compare MF time interval 2 and 3 -------------------------------------------
T3_sig_proteins <- data.frame("Time" = "Time interval 2", "EntrezGeneID" = SigProteins_T3$EntrezGeneID)
T4_sig_proteins <- data.frame("Time" = "Time interval 3", "EntrezGeneID" = SigProteins_T4$EntrezGeneID)

T3T4_sig_proteins <- rbind(T3_sig_proteins, T4_sig_proteins)

T3T4_compare_MF <- compareCluster(EntrezGeneID~Time, data=T3T4_sig_proteins, 
                                  fun="enrichGO", ont = "MF", 
                                  OrgDb = "org.Hs.eg.db", universe = background$EntrezGeneID,
                                  pvalueCutoff = 1, qvalueCutoff = 0.15, readable = TRUE)

Compare_MF <- dotplot(T3T4_compare_MF, x = "Time", 
                      showCategory = 1000,
                      color = "qvalue") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + #Adjusts the length of the labels on y-axis)
  xlab(" ")

ggsave(filename= "CompareT2T3_MF_dotplot_q0.15.png",
       plot = Compare_MF,
       device = "png",
       path = "../ORA/",
       width = 25,
       height = 30,
       #height = 25,
       units = "cm")

#Simplyfy
T3T4_compare_MF_simp <- simplify(T3T4_compare_MF)
Compare_MF_simp <- dotplot(T3T4_compare_MF_simp, x = "Time", 
                           showCategory = 1000,
                           color = "qvalue") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + #Adjusts the length of the labels on y-axis)
  xlab(" ")

ggsave(filename= "CompareT2T3_MF_dotplot_q0.15_simp.png",
       plot = Compare_MF_simp,
       device = "png",
       path = "../ORA/",
       width = 25,
       #height = 25,
       height = 20,
       units = "cm")




# ORA of KEGG time interval 2 ---------------------------------------------

#KEGG pathway does not work. There is something wrong with the database connection
#I found this: https://github.com/YuLab-SMU/clusterProfiler/issues/305
#And tried: 
getOption("clusterProfiler.download.method")
install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method","auto")
#It worked!

#Do the ORA
enrich_KEGG <- enrichKEGG(SigProteins_T3$EntrezGeneID, organism = "hsa", keyType = "kegg", universe = background$EntrezGeneID, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
#Change to gene name
enrich_KEGG <- setReadable(enrich_KEGG, "org.Hs.eg.db", "ENTREZID")
#Calculate rich factor
enrich_KEGG <- mutate(enrich_KEGG, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 

#Extract and save the resutls
enrich_KEGG_results <- enrich_KEGG@result
write_xlsx(enrich_KEGG_results, path = "../ORA/KEGG_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)

enrich_KEGG_sig <- enrich_KEGG %>%
  filter(p.adjust < 0.1)
enrich_KEGG_sig_res <- enrich_KEGG_sig@result
#3 significant terms

#Plot the results in a dotplot
dotplot_KEGG <- dotplot(enrich_KEGG_sig, x = "richFactor", showCategory = 3,
                        title = "KEGG pathway") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) #Adjusts the length of the labels on y-axis

ggsave(filename= "KEGG_dotplot_ORA_adj_cohort_sig.png",
       plot = dotplot_KEGG,
       device = "png",
       path = "../Enrichment/ORA",
       width = 20,
       height = 10,
       units = "cm")

#Make treeplot
enrichgo_KEGG_sig_clust <- pairwise_termsim(enrich_KEGG_sig)
treeplot_enrichgo_KEGG_sig <- treeplot(enrichgo_KEGG_sig_clust, showCategory = 3,
                                       title = "KEGG pathway",
                                       by = "richFactor")

ggsave(filename= "KEGG_treeplot_ORA_adj_cohort_sig.png",
       plot = treeplot_enrichgo_KEGG_sig,
       device = "png",
       path = "../Enrichment/ORA/",
       width = 35,
       height = 10,
       units = "cm")


# ORA of reactome pathway time interval 2 ---------------------------------
T3_enrichRA <- enrichPathway(SigProteins_T3$EntrezGeneID, organism = "human", 
                             pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                             universe = background$EntrezGeneID, readable = TRUE) 

#Calculate richFactor
T3_enrichRA <- mutate(T3_enrichRA, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 


#Extract the results 
T3_enrichRA_result <- T3_enrichRA@result
write_xlsx(T3_enrichRA_result, path = "../ORA/T3_RA_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T3_enrichRA_res_sig <- T3_enrichRA %>%
  filter(p.adjust < 0.2)
dim(T3_enrichRA_res_sig)
#9 sig terms (<0.1)
#12 sig terms (<0.15)
#21 sig terms (<0.2)


#Plot the results in a dotplot
T3_dotplot_enrichRA <- dotplot(T3_enrichRA, x = "richFactor", showCategory = 12,
                               title = "Reactome pathway") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) #Adjusts the length of the labels on y-axis

ggsave(filename= "T3_RA_dotplot_ORA_adj_cohort_q0.15.png",
       plot = T3_dotplot_enrichRA,
       device = "png",
       path = "../ORA/",
       width = 25,
       height = 10,
       units = "cm")

#Make treeplot
T3_enrichgo_RA_clust <- pairwise_termsim(T3_enrichRA)
T3_treeplot_enrichgo_RA <- treeplot(T3_enrichgo_RA_clust, showCategory = 12,
                                    title = "Reactome pathway",
                                    by = "richFactor")

ggsave(filename= "T3_RA_treeplot_ORA_adj_cohort_q0.15.png",
       plot = T3_treeplot_enrichgo_RA,
       device = "png",
       path = "../ORA/",
       width = 40,
       height = 15,
       units = "cm")

# ORA of reactome pathway time interval 3 ---------------------------------
T4_enrichRA <- enrichPathway(SigProteins_T4$EntrezGeneID, organism = "human", 
                             pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                             universe = background$EntrezGeneID, readable = TRUE) 

#Calculate richFactor
T4_enrichRA <- mutate(T4_enrichRA, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 


#Extract the results 
T4_enrichRA_result <- T4_enrichRA@result
write_xlsx(T4_enrichRA_result, path = "../ORA/T4_RA_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T4_enrichRA_res_sig <- T4_enrichRA %>%
  filter(p.adjust < 0.1)
dim(T4_enrichRA_res_sig)
#6 sig terms (<0.1)
#19 sig terms (<0.15)
#35 sig terms (<0.2)


#Plot the results in a dotplot
T4_dotplot_enrichRA <- dotplot(T4_enrichRA, x = "richFactor", showCategory = 19,
                               title = "Reactome pathway") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) #Adjusts the length of the labels on y-axis

ggsave(filename= "T4_RA_dotplot_ORA_adj_cohort_q0.15.png",
       plot = T4_dotplot_enrichRA,
       device = "png",
       path = "../ORA/",
       width = 25,
       height = 16,
       units = "cm")

#Make treeplot
T4_enrichgo_RA_clust <- pairwise_termsim(T4_enrichRA)
T4_treeplot_enrichgo_RA <- treeplot(T4_enrichgo_RA_clust, showCategory = 19,
                                    title = "Reactome pathway",
                                    by = "richFactor")

ggsave(filename= "T4_RA_treeplot_ORA_adj_cohort_q0.15.png",
       plot = T4_treeplot_enrichgo_RA,
       device = "png",
       path = "../ORA/",
       width = 43,
       height = 15,
       units = "cm")





# Compare RA time interval 2 and 3 ----------------------------------------
T3_sig_proteins <- data.frame("Time" = "Time interval 2", "EntrezGeneID" = SigProteins_T3$EntrezGeneID)
T4_sig_proteins <- data.frame("Time" = "Time interval 3", "EntrezGeneID" = SigProteins_T4$EntrezGeneID)

T3T4_sig_proteins <- rbind(T3_sig_proteins, T4_sig_proteins)

T3T4_compare_RA <- compareCluster(EntrezGeneID~Time, data=T3T4_sig_proteins, 
                                  fun="enrichPathway", universe = background$EntrezGeneID,
                                  pvalueCutoff = 1, qvalueCutoff = 0.15, readable = TRUE, pAdjustMethod = "BH")

T3T4_compare_RA_result <- T3T4_compare_RA@result

Compare_RA <- dotplot(T3T4_compare_RA, x = "Time", 
                      showCategory = 1000,
                      color = "qvalue") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + #Adjusts the length of the labels on y-axis)
  xlab(" ")

ggsave(filename= "CompareT2T3_RA_dotplot_q0.15.png",
       plot = Compare_RA,
       device = "png",
       path = "../ORA/",
       width = 25,
       #height = 30,
       height = 28,
       units = "cm")

# ORA of wikipathways time interval 2 ---------------------------------
getOption("clusterProfiler.download.method")
install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method","auto")

T3_enrichWP <- enrichWP(gene = SigProteins_T3$EntrezGeneID, organism = "Homo sapiens",  
                        pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                        universe = background$EntrezGeneID) 

#Calculate richFactor
T3_enrichRA <- mutate(T3_enrichRA, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 


#Extract the results 
T3_enrichRA_result <- T3_enrichRA@result
write_xlsx(T3_enrichRA_result, path = "../ORA/T3_RA_res_ORA_adj_cohort.xlsx",
           col_names = TRUE)
T3_enrichRA_res_sig <- T3_enrichRA %>%
  filter(p.adjust < 0.2)
dim(T3_enrichRA_res_sig)
#9 sig terms (<0.1)
#12 sig terms (<0.15)
#21 sig terms (<0.2)


#Plot the results in a dotplot
T3_dotplot_enrichRA <- dotplot(T3_enrichRA, x = "richFactor", showCategory = 12,
                               title = "Reactome pathway") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) #Adjusts the length of the labels on y-axis

ggsave(filename= "T3_RA_dotplot_ORA_adj_cohort_q0.15.png",
       plot = T3_dotplot_enrichRA,
       device = "png",
       path = "../ORA/",
       width = 25,
       height = 10,
       units = "cm")

#Make treeplot
T3_enrichgo_RA_clust <- pairwise_termsim(T3_enrichRA)
T3_treeplot_enrichgo_RA <- treeplot(T3_enrichgo_RA_clust, showCategory = 12,
                                    title = "Reactome pathway",
                                    by = "richFactor")

ggsave(filename= "T3_RA_treeplot_ORA_adj_cohort_q0.15.png",
       plot = T3_treeplot_enrichgo_RA,
       device = "png",
       path = "../ORA/",
       width = 40,
       height = 15,
       units = "cm")


# Compare DO time interval 2 and 3 ----------------------------------------
T3_sig_proteins <- data.frame("Time" = "Time interval 2", "EntrezGeneID" = SigProteins_T3$EntrezGeneID)
T4_sig_proteins <- data.frame("Time" = "Time interval 3", "EntrezGeneID" = SigProteins_T4$EntrezGeneID)

T3T4_sig_proteins <- rbind(T3_sig_proteins, T4_sig_proteins)

T3T4_compare_DO <- compareCluster(EntrezGeneID~Time, data=T3T4_sig_proteins, 
                                  fun="enrichDO", universe = background$EntrezGeneID,
                                  pvalueCutoff = 1, qvalueCutoff = 0.15, readable = TRUE, pAdjustMethod = "BH")

Compare_DO <- dotplot(T3T4_compare_DO, x = "Time", 
                      showCategory = 1000,
                      color = "qvalue") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) + #Adjusts the length of the labels on y-axis)
  xlab(" ")

ggsave(filename= "CompareT2T3_DO_dotplot_q0.15.png",
       plot = Compare_DO,
       device = "png",
       path = "../ORA/",
       width = 25,
       #height = 30,
       height = 28,
       units = "cm")


