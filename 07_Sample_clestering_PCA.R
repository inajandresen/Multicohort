library(tidyverse)
library(Biobase)
library(pcaMethods)
library(scales)
library(gridExtra)
library(ggpubr)
library(grid)


# Before detrending -------------------------------------------------------

#Time interval 1
load("../Data/Merged_cohorts_wins_lgRFU_exprset.RData")
T1_Merged_cohorts_wins_lgRFU_exprset <- Merged_cohorts_wins_lgRFU_exprset[, Merged_cohorts_wins_lgRFU_exprset$Timepoint %in% c("T2")]
T1_Merged_cohorts_wins_lgRFU_exprset$Group <- gsub("PE", "LOPE", T1_Merged_cohorts_wins_lgRFU_exprset$Group)

pc_mc_lg_t1 <- pca(T1_Merged_cohorts_wins_lgRFU_exprset)
df_mc_lg_t1 <- merge(scores(pc_mc_lg_t1), pData(T1_Merged_cohorts_wins_lgRFU_exprset), by = 0)
pca_mc_lg_t1 <- ggplot(df_mc_lg_t1, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 1", subtitle = "Weeks 12-19") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t1@R2[1], digits = 1),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t1@R2[2], digits = 1),"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")

#Time interval 2
load("../Data/Merged_cohorts_wins_lgRFU_exprset.RData")
T2_Merged_cohorts_wins_lgRFU_exprset <- Merged_cohorts_wins_lgRFU_exprset[, Merged_cohorts_wins_lgRFU_exprset$Timepoint %in% c("T3")]
T2_Merged_cohorts_wins_lgRFU_exprset$Group <- gsub("PE", "LOPE", T2_Merged_cohorts_wins_lgRFU_exprset$Group)

pc_mc_lg_t2 <- pca(T2_Merged_cohorts_wins_lgRFU_exprset)
df_mc_lg_t2 <- merge(scores(pc_mc_lg_t2), pData(T2_Merged_cohorts_wins_lgRFU_exprset), by = 0)
pca_mc_lg_t2 <- ggplot(df_mc_lg_t2, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 2", subtitle = "Weeks 19-27") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t2@R2[1], digits = 1),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t2@R2[2], digits = 2),"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")

#Time interval 3
load("../Data/Merged_cohorts_wins_lgRFU_exprset.RData")
T3_Merged_cohorts_wins_lgRFU_exprset <- Merged_cohorts_wins_lgRFU_exprset[, Merged_cohorts_wins_lgRFU_exprset$Timepoint %in% c("T4")]
T3_Merged_cohorts_wins_lgRFU_exprset$Group <- gsub("PE", "LOPE", T3_Merged_cohorts_wins_lgRFU_exprset$Group)

pc_mc_lg_t3 <- pca(T3_Merged_cohorts_wins_lgRFU_exprset)
df_mc_lg_t3 <- merge(scores(pc_mc_lg_t3), pData(T3_Merged_cohorts_wins_lgRFU_exprset), by = 0)
pca_mc_lg_t3 <- ggplot(df_mc_lg_t3, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 3", subtitle = "Weeks 27-34") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t3@R2[1], digits = 1),"%)")) + 
  #ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t3@R2[2], digits = 3),"%)")) +
  ylab(paste0("PC2 (R2 = 7.00%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")


# After detrending --------------------------------------------------------

#Time interval 1
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
T1_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2")]
Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", Merged_cohorts_wins_Detrended_exprset$Group)

pc_mc_lg_t1_dt <- pca(Merged_cohorts_wins_Detrended_exprset)
df_mc_lg_t1_dt <- merge(scores(pc_mc_lg_t1_dt), pData(Merged_cohorts_wins_Detrended_exprset), by = 0)
pca_mc_lg_t1_dt <- ggplot(df_mc_lg_t1_dt, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 1", subtitle = "Weeks 12-19") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t1_dt@R2[1], digits = 1),"%)")) + 
  #ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t1_dt@R2[2], digits = 2),"%)")) +
  ylab(paste0("PC2 (R2 = 5.90%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")

#Time interval 2
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
T2_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T3")]
T2_Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", T2_Merged_cohorts_wins_Detrended_exprset$Group)

pc_mc_lg_t2_dt <- pca(T2_Merged_cohorts_wins_Detrended_exprset)
df_mc_lg_t2_dt <- merge(scores(pc_mc_lg_t2_dt), pData(T2_Merged_cohorts_wins_Detrended_exprset), by = 0)
pca_mc_lg_t2_dt <- ggplot(df_mc_lg_t2_dt, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 2", subtitle = "Weeks 19-27") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t2_dt@R2[1], digits = 1),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t2_dt@R2[2], digits = 2),"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")

#Time interval 3
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
T3_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T4")]
T3_Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", T3_Merged_cohorts_wins_Detrended_exprset$Group)

pc_mc_lg_t3_dt <- pca(T3_Merged_cohorts_wins_Detrended_exprset)
df_mc_lg_t3_dt <- merge(scores(pc_mc_lg_t3_dt), pData(T3_Merged_cohorts_wins_Detrended_exprset), by = 0)
pca_mc_lg_t3_dt <- ggplot(df_mc_lg_t3_dt, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 3", subtitle = "Weeks 27-34") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t3_dt@R2[1], digits = 1),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t3_dt@R2[2], digits = 2),"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")


# Legend for cohorts ------------------------------------------------------

#Function for extracting legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

for_legend_cohorts <- ggplot(df_mc_lg_t3_dt, aes(PC1, PC2, color = Cohort)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF" )) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 3", subtitle = "Weeks 27-34") +
  xlab(paste0("PC1 (R2 = ",100*pc_mc_lg_t3_dt@R2[1],"%)")) + 
  ylab(paste0("PC2 (R2 = ",100*pc_mc_lg_t3_dt@R2[2],"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

legend_cohorts <- get_legend(for_legend_cohorts)


# After detrending, PE --------------------------------------------------------

#Time interval 1
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
T1_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2")]
Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", Merged_cohorts_wins_Detrended_exprset$Group)

pc_mc_lg_t1_dt_pe <- pca(Merged_cohorts_wins_Detrended_exprset)
df_mc_lg_t1_dt_pe <- merge(scores(pc_mc_lg_t1_dt_pe), pData(Merged_cohorts_wins_Detrended_exprset), by = 0)
pca_mc_lg_t1_dt_pe <- ggplot(df_mc_lg_t1_dt_pe, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#3399CC", "#000066")) +
  stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 1", subtitle = "Weeks 12-19") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t1_dt_pe@R2[1], digits = 1),"%)")) + 
  #ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t1_dt_pe@R2[2], digits = 2),"%)")) +
  ylab(paste0("PC2 (R2 = 5.90%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")

#Time interval 2
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
T2_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T3")]
T2_Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", T2_Merged_cohorts_wins_Detrended_exprset$Group)

pc_mc_lg_t2_dt_pe <- pca(T2_Merged_cohorts_wins_Detrended_exprset)
df_mc_lg_t2_dt_pe <- merge(scores(pc_mc_lg_t2_dt_pe), pData(T2_Merged_cohorts_wins_Detrended_exprset), by = 0)
pca_mc_lg_t2_dt_pe <- ggplot(df_mc_lg_t2_dt_pe, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#3399CC", "#000066")) +
  stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 2", subtitle = "Weeks 19-27") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t2_dt_pe@R2[1], digits = 1),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t2_dt_pe@R2[2], digits = 2),"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")

#Time interval 3
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
T3_Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T4")]
T3_Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", T3_Merged_cohorts_wins_Detrended_exprset$Group)

pc_mc_lg_t3_dt_pe <- pca(T3_Merged_cohorts_wins_Detrended_exprset)
df_mc_lg_t3_dt_pe <- merge(scores(pc_mc_lg_t3_dt_pe), pData(T3_Merged_cohorts_wins_Detrended_exprset), by = 0)
pca_mc_lg_t3_dt_pe <- ggplot(df_mc_lg_t3_dt_pe, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#3399CC", "#000066")) +
  stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 3", subtitle = "Weeks 27-34") +
  xlab(paste0("PC1 (R2 = ",round(100*pc_mc_lg_t3_dt_pe@R2[1], digits = 1),"%)")) + 
  ylab(paste0("PC2 (R2 = ",round(100*pc_mc_lg_t3_dt_pe@R2[2], digits = 2),"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "none")


# Legend PE ---------------------------------------------------------------

#Function for extracting legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

for_legend_pe <- ggplot(df_mc_lg_t3_dt_pe, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#3399CC", "#000066")) +
  #stat_ellipse(type = "t") +
  #labs(title = "Before detrending") + 
  labs(title = "Time interval 3", subtitle = "Weeks 27-34") +
  xlab(paste0("PC1 (R2 = ",100*pc_mc_lg_t3_dt_pe@R2[1],"%)")) + 
  ylab(paste0("PC2 (R2 = ",100*pc_mc_lg_t3_dt_pe@R2[2],"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  #labs(tag = "A)") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

legend_pe <- get_legend(for_legend_pe)


# Arrange figure ----------------------------------------------------------
Before_detrending_title <- textGrob("Before de-trending", gp = gpar(fontsize = 20))
After_detrending_title <- textGrob("After de-trending", gp = gpar(fontsize = 20))

pca_plots_b_dt <- ggarrange(pca_mc_lg_t1, pca_mc_lg_t2, pca_mc_lg_t3, legend_cohorts,
                            nrow = 1,
                            widths = c(1, 1, 1, 0.3),
                            labels = c("A)", " ", " ", " "),
                            font.label = list(size = 20))

pca_plots_a_dt <- ggarrange(pca_mc_lg_t1_dt, pca_mc_lg_t2_dt, pca_mc_lg_t3_dt, legend_cohorts,
                            nrow = 1,
                            widths = c(1, 1, 1, 0.3),
                            labels = c("B)", " ", " ", " "),
                            font.label = list(size = 20))

pca_plots_a_dt_pe <- ggarrange(pca_mc_lg_t1_dt_pe, pca_mc_lg_t2_dt_pe, pca_mc_lg_t3_dt_pe, legend_pe,
                               nrow = 1,
                               widths = c(1, 1, 1, 0.3),
                               labels = c("C)", " ", " ", " "),
                               font.label = list(size = 20))

pca_plots <- ggarrange(Before_detrending_title, pca_plots_b_dt, 
                       After_detrending_title, pca_plots_a_dt, pca_plots_a_dt_pe, 
                       ncol = 1,
                       heights = c(0.2, 1, 0.2, 1, 1))

ggsave(filename= "PCA_multicohort_review.png",
       plot = pca_plots,
       device = "png",
       path = "../Sample_clustering/",
       width = 35,
       height = 30,
       units = "cm",
       bg = "white",
       dpi = 800)
