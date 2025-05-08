library(tidyverse)
library(Biobase)
library(readxl)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(grid)


# Load and format data ----------------------------------------------------

#Load the significant proteins from limma test of the multicohort
T2 <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T2 <- T2 %>%
  filter(adj.P.Val < 0.1)
T2 <- T2$SomaId
T3 <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T3 <- T3 %>%
  filter(adj.P.Val < 0.1)
T3 <- T3$SomaId
T4 <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")
T4 <- T4 %>%
  filter(adj.P.Val < 0.1)
T4 <- T4$SomaId

sig_proteins_multi <- c(T3, T4) %>%
  unique()

rm(T2, T3, T4)

#Load the oslo data and filter on the significant proteins 
T2_oslo <- read_xlsx("../Limma_sep_cohorts/T2_oslo.xlsx",
                        col_names = TRUE)
T2_oslo <- filter(T2_oslo, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T2_oslo <- T2_oslo %>% 
  dplyr::select(AptName, SomaId, TargetFullName, EntrezGeneSymbol, logFC, adj.P.Val) %>%
  rename("Oslo_logFC" = logFC) %>%
  rename("Oslo_adj.P.Val" = adj.P.Val)


T3_oslo <- read_xlsx("../Limma_sep_cohorts/T3_oslo.xlsx",
                        col_names = TRUE)
T3_oslo <- filter(T3_oslo, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T3_oslo <- T3_oslo %>% 
  dplyr::select(AptName, SomaId, TargetFullName, EntrezGeneSymbol, logFC, adj.P.Val) %>%
  rename("Oslo_logFC" = logFC) %>%
  rename("Oslo_adj.P.Val" = adj.P.Val)

T4_oslo <- read_xlsx("../Limma_sep_cohorts/T4_oslo.xlsx",
                        col_names = TRUE)
T4_oslo <- filter(T4_oslo, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T4_oslo <- T4_oslo %>% 
  dplyr::select(AptName, SomaId, TargetFullName, EntrezGeneSymbol, logFC, adj.P.Val) %>%
  rename("Oslo_logFC" = logFC) %>%
  rename("Oslo_adj.P.Val" = adj.P.Val)


#Load the stanford results and filter on significant proteins from the multicohort
T2_stanford <- read_xlsx("../Limma_sep_cohorts/T2_stanford.xlsx",
                        col_names = TRUE)
T2_stanford <- filter(T2_stanford, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T2_stanford <- T2_stanford %>% 
  dplyr::select(SomaId, logFC, adj.P.Val) %>%
  rename("Stanford_logFC" = logFC) %>%
  rename("Stanford_adj.P.Val" = adj.P.Val)

T3_stanford <- read_xlsx("../Limma_sep_cohorts/T3_stanford.xlsx",
                        col_names = TRUE)
T3_stanford <- filter(T3_stanford, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T3_stanford <- T3_stanford %>% 
  dplyr::select(SomaId, logFC, adj.P.Val) %>%
  rename("Stanford_logFC" = logFC) %>%
  rename("Stanford_adj.P.Val" = adj.P.Val)

T4_stanford <- read_xlsx("../Limma_sep_cohorts/T4_stanford.xlsx",
                        col_names = TRUE)
T4_stanford <- filter(T4_stanford, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T4_stanford <- T4_stanford %>% 
  dplyr::select(SomaId, logFC, adj.P.Val) %>%
  rename("Stanford_logFC" = logFC) %>%
  rename("Stanford_adj.P.Val" = adj.P.Val)

#Load the detroit results and filter on the significant proteins from the multicohort
T2_detroit <- read_xlsx("../Limma_sep_cohorts/T2_detroit.xlsx",
           col_names = TRUE)
T2_detroit <- filter(T2_detroit, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T2_detroit <- T2_detroit %>% 
  dplyr::select(SomaId, logFC, adj.P.Val) %>%
  rename("Detroit_logFC" = logFC) %>%
  rename("Detroit_adj.P.Val" = adj.P.Val)

T3_detroit <- read_xlsx("../Limma_sep_cohorts/T3_detroit.xlsx",
           col_names = TRUE)
T3_detroit <- filter(T3_detroit, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T3_detroit <- T3_detroit %>% 
  dplyr::select(SomaId, logFC, adj.P.Val) %>%
  rename("Detroit_logFC" = logFC) %>%
  rename("Detroit_adj.P.Val" = adj.P.Val)

T4_detroit <- read_xlsx("../Limma_sep_cohorts/T4_detroit.xlsx",
           col_names = TRUE)
T4_detroit <- filter(T4_detroit, grepl(paste(sig_proteins_multi, collapse = "|"), SomaId)) #Filter on common SomaIds
T4_detroit <- T4_detroit %>% 
  dplyr::select(SomaId, logFC, adj.P.Val) %>%
  rename("Detroit_logFC" = logFC) %>%
  rename("Detroit_adj.P.Val" = adj.P.Val)

#Mere the data
T2 <- merge(T2_oslo, T2_detroit, by = "SomaId")
T2 <- merge(T2, T2_stanford, by = "SomaId")

T3 <- merge(T3_oslo, T3_detroit, by = "SomaId")
T3 <- merge(T3, T3_stanford, by = "SomaId")

T4 <- merge(T4_oslo, T4_detroit, by = "SomaId")
T4 <- merge(T4, T4_stanford, by = "SomaId")



# Plot correlations -------------------------------------------------------


T2_oslo_detroit <- ggplot(data = T2, aes(x = Oslo_logFC, y = Detroit_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Detroit log Fold Change") +
  xlab("Oslo log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  #scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.25, label.y = -0.62) +
  ggtitle("Oslo vs Detroit")

T2_oslo_stanford <- ggplot(data = T2, aes(x = Oslo_logFC, y = Stanford_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Stanford log Fold Change") +
  xlab("Oslo log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.2, label.y = -1.5) +
  ggtitle("Oslo vs Stanford")

T2_detroit_stanford <- ggplot(data = T2, aes(x = Detroit_logFC, y = Stanford_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Stanford log Fold Change") +
  xlab("Detroit log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.1, label.y = -1.5) +
  ggtitle("Stanford vs Detroit")


#Correlation plots, time interval two 
T3_oslo_detroit <- ggplot(data = T3, aes(x = Oslo_logFC, y = Detroit_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Detroit log Fold Change") +
  xlab("Oslo log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.1, label.y = -0.7) +
  ggtitle("Oslo vs Detroit")

T3_oslo_stanford <- ggplot(data = T3, aes(x = Oslo_logFC, y = Stanford_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Stanford log Fold Change") +
  xlab("Oslo log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.1, label.y = -1.5) +
  ggtitle("Oslo vs Stanford") 

T3_detroit_stanford <- ggplot(data = T3, aes(x = Detroit_logFC, y = Stanford_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Stanford log Fold Change") +
  xlab("Detroit log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.4, label.y = -1.5) +
  ggtitle("Stanford vs Detroit")


#Correlation plots, time interval three 
T4_oslo_detroit <- ggplot(data = T4, aes(x = Oslo_logFC, y = Detroit_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Detroit log Fold Change") +
  xlab("Oslo log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
   stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.4, label.y = -0.7) +
  ggtitle("Oslo vs Detroit")

T4_oslo_stanford <- ggplot(data = T4, aes(x = Oslo_logFC, y = Stanford_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Stanford log Fold Change") +
  xlab("Oslo log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
   stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.5, label.y = -1.2) +
  ggtitle("Oslo vs Stanford")

T4_detroit_stanford <- ggplot(data = T4, aes(x = Detroit_logFC, y = Stanford_logFC), aplha = 0.5) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "lightgrey", size=0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_text_repel(aes(label = EntrezGeneSymbol),
                  #label.padding = 1, 
                  point.padding = 0.5,
                  max.overlaps = 12,
                  segment.color = 'grey50',
                  size = 3.3) +
  ylab("Stanford log Fold Change") +
  xlab("Detroit log Fold Change") +
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=10)) +
  stat_cor(method = "pearson", r.digits = 2, p.digits = 2, color = "blue",
           label.x = 0.4, label.y = -1.3) +
  ggtitle("Stanford vs Detroit")


# Save plots --------------------------------------------------------------

title_T2 <- textGrob("Time interval 1\nweeks 12-19", gp = gpar(fontsize = 20))
title_T3 <- textGrob("Time interval 2\nweeks 19-27", gp = gpar(fontsize = 20))
title_T4 <- textGrob("Time interval 3\nweeks 27-34", gp = gpar(fontsize = 20))

Row1 <- grid.arrange(title_T2, T2_oslo_detroit, T2_oslo_stanford, T2_detroit_stanford, 
                     ncol = 1, nrow = 4, heights = c(0.2, 1, 1, 1))

Row2 <- grid.arrange(title_T3, T3_oslo_detroit, T3_oslo_stanford, T3_detroit_stanford, 
                     ncol = 1, nrow = 4, heights = c(0.2, 1, 1, 1))

Row3 <- grid.arrange(title_T4, T4_oslo_detroit, T4_oslo_stanford, T4_detroit_stanford, 
                     ncol = 1, nrow = 4, heights = c(0.2, 1, 1, 1))


Plots <- grid.arrange(Row1, Row2, Row3, ncol = 3)

ggsave(filename= "Correlations_foldchange_pearson.png",
       plot = Plots,
       device = "png",
       path = "../Limma_sep_cohorts/",
       width = 35,
       height = 35,
       units = "cm")


#Plot to main figure 
Plot_main <- grid.arrange(T3_oslo_detroit, T3_oslo_stanford, T3_detroit_stanford, 
                     ncol = 3, nrow = 1)


ggsave(filename= "Correlations_foldchange_pearson_T3.png",
       plot = Plot_main,
       device = "png",
       path = "../Limma_sep_cohorts/",
       width = 35,
       height = 13,
       units = "cm")

#Plot to main figure 
Plot_main <- grid.arrange(T4_oslo_detroit, T4_oslo_stanford, T4_detroit_stanford, 
                          ncol = 3, nrow = 1)


ggsave(filename= "Correlations_foldchange_pearson_T4.png",
       plot = Plot_main,
       device = "png",
       path = "../Limma_sep_cohorts/",
       width = 35,
       height = 13,
       units = "cm")


