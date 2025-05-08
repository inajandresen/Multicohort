library(tidyverse)
library(pheatmap)
library(readxl)
library(gridExtra)
library(ggrepel)
library(writexl)
library(grid)
library(ggvenn)
library(ggpubr)
library(ggplotify)


# Histogram of p-values ---------------------------------------------------

T2 <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T3 <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T4 <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")

T2_hist <- ggplot(T2, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Time interval 1", subtitle = "Weeks 12-19") +
  theme_classic()

T3_hist <- ggplot(T3, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Time interval 2", subtitle = "Weeks 19-27") +
  theme_classic()

T4_hist <- ggplot(T4, aes(P.Value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "black", colour = "white") +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Time interval 3", subtitle = "Weeks 27-34") +
  theme_classic()



histogram <- grid.arrange(T2_hist, T3_hist, T4_hist, nrow = 1)

ggsave(filename= "Histogram_p-values_adj_cohort.png",
       plot = histogram,
       device = "png",
       path = "../Limma_adj_cohort/",
       width = 20,
       height = 8,
       units = "cm")


# Venn diagram ------------------------------------------------------------
T2 <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T2 <- T2 %>%
  filter(adj.P.Val < 0.1)

T3 <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T3 <- T3 %>%
  filter(adj.P.Val < 0.1)

T4 <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")
T4 <- T4 %>%
  filter(adj.P.Val < 0.1)



T2T3T4_venn <- list("Time interval 1" = T2$AptName, "Time interval 2" = T3$AptName, "Time interval 3" = T4$AptName)

vennplot <- ggvenn(T2T3T4_venn, 
                   fill_color = c("#A3A500", "#00B0F6", "#E76BF3"),
                   stroke_size = 0.5, 
                   set_name_size = 5,
                   show_percentage = FALSE,
                   text_size = 4.3)

ggsave(filename= "Venn_limma_adj_cohort_pval_0.1.png",
       plot = vennplot,
       device = "png",
       path = "../Limma_adj_cohort/",
       width = 10,
       height = 10,
       units = "cm",
       bg = "white")

overlap_T3_T4 <- intersect(T3$EntrezGeneSymbol, T4$EntrezGeneSymbol)


# Volcanoplot, qvalue < 0.1 -----------------
load("../Limma_adj_cohort/fit_adj_cohort.RData")
rm(fit)
rm(design)

#Extract all proteins
T2_all <- topTable(efit, coef = "T2", number = "all", adjust.method="BH")
T3_all <- topTable(efit, coef = "T3", number = "all", adjust.method="BH")
T4_all <- topTable(efit, coef = "T4", number = "all", adjust.method="BH")


#Add a column displaying the differential expressed proteins to the data frame
T2_all$DE_0.05 <- "NO"
T3_all$DE_0.05 <- "NO"
T4_all$DE_0.05 <- "NO"


#Set the threshould for upregulated proteins
T2_all$DE_0.05[T2_all$adj.P.Val < 0.1 & T2_all$logFC > 0] <- "UP"
T3_all$DE_0.05[T3_all$adj.P.Val < 0.1 & T3_all$logFC > 0] <- "UP"
T4_all$DE_0.05[T4_all$adj.P.Val < 0.1 & T4_all$logFC > 0] <- "UP"

#Set the threshould for downregulated proteins
T2_all$DE_0.05[T2_all$adj.P.Val < 0.1 & T2_all$logFC < 0] <- "DOWN"
T3_all$DE_0.05[T3_all$adj.P.Val < 0.1 & T3_all$logFC < 0] <- "DOWN"
T4_all$DE_0.05[T4_all$adj.P.Val < 0.1 & T4_all$logFC < 0] <- "DOWN"

#Add labels
T2_all$Label <- NA
T2_all$Label[T2_all$adj.P.Val < 0.1] <- T2_all$EntrezGeneSymbol[T2_all$adj.P.Val < 0.1]
#T2_all$Label[T2_all$logFC > 0.25] <- T2_all$EntrezGeneSymbol[T2_all$logFC > 0.25]
#T2_all$Label[T2_all$logFC < -0.4] <- T2_all$EntrezGeneSymbol[T2_all$logFC < -0.4]

T3_all$Label <- NA
T3_all$Label[T3_all$adj.P.Val < 0.1] <- T3_all$EntrezGeneSymbol[T3_all$adj.P.Val < 0.1]
#T3_all$Label[T3_all$logFC < -0.4] <- T3_all$EntrezGeneSymbol[T3_all$logFC < -0.4]

T4_all$Label <- NA
T4_all$Label[T4_all$adj.P.Val < 0.1] <- T4_all$EntrezGeneSymbol[T4_all$adj.P.Val < 0.1]
#T4_all$Label[T4_all$logFC < -0.4] <- T4_all$EntrezGeneSymbol[T4_all$logFC < -0.4]

#Plot the volcanoplot
T2_0.1 <- ggplot(data = T2_all, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  xlim(-0.8, 1) +
  ylim(0, 6.5) +
  #labs(tag = "A)") +
  xlab("log2(Fold change)") +
  ylab("-log10(q-value)") +
  geom_text_repel(aes(label = Label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  theme_classic() +
  scale_color_manual(values = c("grey", "brown2")) +
  theme(legend.position = "none") +
  labs(title = "Time interval 1", face = "bold", subtitle = "Weeks 12-19") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


T3_0.1 <- ggplot(data = T3_all, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  xlim(-0.8, 1) +
  ylim(0, 6.5) +
  #labs(tag = "B)") +
  geom_text_repel(aes(label = Label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  theme_classic() +
  scale_color_manual(values = c("cornflowerblue", "grey", "brown2")) +
  #geom_hline(yintercept = -log10(0.05), col = "red") +
  theme(legend.position = "none") +
  labs(title = "Time interval 2", face = "bold", subtitle ="Weeks 19-27") +
  xlab("log2(Fold change)") +
  ylab("-log10(q-value)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


T4_0.1 <- ggplot(data = T4_all, aes(x = logFC, y = -log10(adj.P.Val), col = DE_0.05)) +
  geom_point() +
  xlim(-0.8, 1) +
  ylim(0, 6.5) +
  #labs(tag = "C)") +
  geom_text_repel(aes(label = Label),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  xlab("log2(Fold change)") +
  ylab("-log10(q-value)") +
  theme_classic() +
  scale_color_manual(values = c("cornflowerblue", "grey", "brown2")) +
  #geom_hline(yintercept = -log10(0.05), col = "red") +
  theme(legend.position = "none") +
  labs(title = "Time interval 3", face = "bold", subtitle ="Weeks 27-34") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))



volcanoplot <- grid.arrange(T2_0.1, T3_0.1, T4_0.1, nrow = 2)

ggsave(filename= "Volcanoplot_adj_cohort_pval_0.1.png",
       plot = volcanoplot,
       device = "png",
       path = "../Limma_adj_cohort/",
       width = 23,
       height = 23,
       units = "cm")


# Heat map of DE proteins, pval >0.1  ---------------------------------------------------
T2_all <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T3_all <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T4_all <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")

T2 <- T2_all %>%
  dplyr::select(c(logFC, SomaId, EntrezGeneSymbol, TargetFullName)) %>%
  dplyr::rename("T2" = logFC)

T3 <- T3_all %>%
  dplyr::select(c(logFC, SomaId)) %>%
  dplyr::rename("T3" = logFC)

T4 <- T4_all %>%
  dplyr::select(c(logFC, SomaId)) %>%
  dplyr::rename("T4" = logFC)

LogFC <- merge(T2, T3, by = "SomaId", all = "TRUE") %>%
  merge(T4, by = "SomaId", all = "TRUE") 


#Filter on the significant proteins 
SigProteins_T2 <- T2_all %>%
  filter(adj.P.Val < 0.1)
SigProteins_T3 <- T3_all %>%
  filter(adj.P.Val < 0.1)
SigProteins_T4 <- T4_all %>%
  filter(adj.P.Val < 0.1)

SigProteins <- rbind(SigProteins_T2, SigProteins_T3)
SigProteins <- rbind(SigProteins, SigProteins_T4)
SigProteins <- SigProteins$SomaId
SigProteins <- unique(SigProteins)

LogFC2 <- filter(LogFC, grepl(paste(SigProteins, collapse = "|"), SomaId)) #Filter on common SomaIds
LogFC2 <- LogFC2 %>%
  dplyr::mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
LogFC2$TargetFullName_EntrezGeneSymbol <- paste(LogFC2$TargetFullName, LogFC2$EntrezGeneSymbol2)

#Check duplicates, Targetfullname
TargetFullName_EntrezGeneSymbol <- LogFC2$TargetFullName_EntrezGeneSymbol #Extract to check for duplicate names
dup <- TargetFullName_EntrezGeneSymbol[duplicated(TargetFullName_EntrezGeneSymbol)] #Extract the duplicates


LogFC2 <- LogFC2 %>% 
  remove_rownames() %>% 
  column_to_rownames(var="TargetFullName_EntrezGeneSymbol") %>%
  #column_to_rownames(var="EntrezGeneSymbol2") %>%
  select(T2, T3, T4)
colnames(LogFC2) <- c("Time interval 1", "Time interval 2", "Time interval 3")

#Set value zero to be white color
Breaks <- c(seq(min(LogFC2), 0, length.out=ceiling(100/2) + 1), 
            seq(max(LogFC2)/100, max(LogFC2), length.out=floor(100/2)))


heatmap <- LogFC2 %>%
  pheatmap(border_color = NA,
           breaks = Breaks,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           color = colorRampPalette(c("cornflowerblue", "white", "brown2"))(100),
           fontsize_row = 16,
           fontsize_col = 15,
           treeheight_row = 0,
           cellwidth = 20,
           main = "")

ggsave(filename= "Heatmap_DE_proteins_adj_cohort_pval_0.1.png",
       plot = heatmap,
       device = "png",
       path = "../Limma_adj_cohort/",
       width = 35,
       height = 42,
       units = "cm")


# List of results ---------------------------------------------------------
T2 <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T2_sig <- T2 %>%
  filter(adj.P.Val < 0.1)
T2_sig <- T2_sig$AptName
T2_list <- T2 %>%
  dplyr::select(SomaId, AptName, TargetFullName, EntrezGeneSymbol, logFC, adj.P.Val) %>%
  dplyr::rename("log FC, T2" = logFC) %>%
  dplyr::rename("q-value, T2" = adj.P.Val)

T3 <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T3_sig <- T3 %>%
  filter(adj.P.Val < 0.1)
T3_sig <- T3_sig$AptName
T3_list <- T3 %>%
  dplyr::select(SomaId, AptName, TargetFullName, EntrezGeneSymbol, logFC, adj.P.Val) %>%
  dplyr::rename("log FC, T3" = logFC) %>%
  dplyr::rename("q-value, T3" = adj.P.Val)

T4 <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")
T4_sig <- T4 %>%
  filter(adj.P.Val < 0.1)
T4_sig <- T4_sig$AptName
T4_list <- T4 %>%
  dplyr::select(SomaId, AptName, TargetFullName, EntrezGeneSymbol, logFC, adj.P.Val) %>%
  dplyr::rename("log FC, T4" = logFC) %>%
  dplyr::rename("q-value, T4" = adj.P.Val)

Limma_results <- merge(T2_list, T3_list, by = c("AptName", "TargetFullName", "EntrezGeneSymbol"), all = TRUE)
Limma_results <- merge(Limma_results, T4_list, by = c("AptName", "TargetFullName", "EntrezGeneSymbol"), all = TRUE)
Limma_results <- Limma_results %>%
  dplyr::select(-SomaId.x, -SomaId.y) %>%
  dplyr::select(SomaId, everything()) 

write_xlsx(Limma_results, path = "../Limma_adj_cohort/Limma_results.xlsx",
           col_names = TRUE)


SigProteins <- c(T2_sig, T3_sig, T4_sig)
SigProteins <- unique(SigProteins)

Limma_results_sigproteins <- filter(Limma_results, grepl(paste(SigProteins, collapse = "|"), AptName)) #Filter on common SomaIds

write_xlsx(Limma_results_sigproteins, path = "../Limma_adj_cohort/Limma_results_sigproteins.xlsx",
           col_names = TRUE)



