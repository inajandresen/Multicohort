library(tidyverse)
library(gridExtra)
library(Biobase)
library(ggpubr)

# Histogram of sampling time ----------------------------------------------
load("../Data/Merged_cohorts_RFU_exprset.RData")

#Separapte the cohorts
Oslo_exprset <- Merged_cohorts_RFU_exprset[, pData(Merged_cohorts_RFU_exprset)$Cohort %in% "Oslo"]
Oslo_samp <- pData(Oslo_exprset)
Detroit_exprset <- Merged_cohorts_RFU_exprset[, pData(Merged_cohorts_RFU_exprset)$Cohort %in% "Detroit"]
Detroit_samp <- pData(Detroit_exprset)
Stanford_exprset <- Merged_cohorts_RFU_exprset[, pData(Merged_cohorts_RFU_exprset)$Cohort %in% "Stanford"]
Stanford_samp <- pData(Stanford_exprset)
Merged_samp <- pData(Merged_cohorts_RFU_exprset)


#Plot histogram to check for normal distrubution; separate the cohorts, vline at trimester
Oslo_hist_trim <- ggplot(data = Oslo_samp, aes(x = GAWeeks)) +
  geom_histogram(binwidth = 0.5, fill = "#00BA38", color = "#00BA38", alpha = 0.7) +
  ggtitle("Oslo") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10)) +
  #ylab("No. patients") +
  #xlab("Weeks of gestation") +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept = 12, size = 1, color = "black") +
  #geom_vline(xintercept = 19, size = 1, color = "black") +
  #geom_vline(xintercept = 27, size = 1, color = "black") +
  #geom_vline(xintercept = 34, size = 1, color = "black") +
  coord_cartesian(xlim = c(7, 42)) +
  ylim(c(0, 35))

ggsave(filename= "Oslo_GAweeks_histograms_hlines_binwidth0.5.png",
       plot = Oslo_hist_trim,
       device = "png",
       path = "../Results/Histograms_timepoints/",
       width = 20,
       height = 15,
       units = "cm")

Detroit_hist_trim <- ggplot(data = Detroit_samp, aes(x = GAWeeks)) +
  geom_histogram(binwidth = 0.5, fill = "#F8766D", color = "#F8766D", alpha = 0.7) +
  ggtitle("Detroit") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10)) +
  #ylab("No. patients") +
  #xlab("Weeks of gestation") +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept = 12, size = 1, color = "black") +
  #geom_vline(xintercept = 19, size = 1, color = "black") +
  #geom_vline(xintercept = 27, size = 1, color = "black") +
  #geom_vline(xintercept = 34, size = 1, color = "black") +
  coord_cartesian(xlim = c(7, 42)) +
  ylim(c(0, 35))

ggsave(filename= "Detroit_GAweeks_histograms_hlines_binwidth0.5.png",
       plot = Detroit_hist_trim,
       device = "png",
       path = "../Results/Histograms_timepoints/",
       width = 20,
       height = 15,
       units = "cm")


Stanford_hist_trim <- ggplot(data = Stanford_samp, aes(x = GAWeeks)) +
  geom_histogram(binwidth = 0.5, fill = "#619CFF", color = "#619CFF", alpha = 0.7) +
  ggtitle("Stanford") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10)) +
  #ylab("No. patients") +
  #xlab("Weeks of gestation") +
  ylab("") +
  xlab("") +
  #geom_vline(xintercept = 12, size = 1, color = "black") +
  #geom_vline(xintercept = 19, size = 1, color = "black") +
  #geom_vline(xintercept = 27, size = 1, color = "black") +
  #geom_vline(xintercept = 34, size = 1, color = "black") +
  coord_cartesian(xlim = c(7, 42)) +
  ylim(c(0, 35))

ggsave(filename= "Stanford_GAweeks_histograms_hlines_binwidth0.5.png",
       plot = Stanford_hist_trim,
       device = "png",
       path = "../Results/Histograms_timepoints/",
       width = 20,
       height = 15,
       units = "cm")

All_hist_trim <- ggplot(data = Merged_samp, aes(x = GAWeeks, fill = Cohort, color = Cohort)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7) +
  #ggtitle("All cohorts") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10)) +
  #ylab("No. of participants") +
  #xlab("Weeks of gestation") +
  ylab("") +
  xlab("") +
  geom_density() +
  #geom_vline(xintercept = 12, size = 1, color = "black") +
  #geom_vline(xintercept = 19, size = 1, color = "black") +
  #geom_vline(xintercept = 27, size = 1, color = "black") +
  #geom_vline(xintercept = 34, size = 1, color = "black") +
  coord_cartesian(xlim = c(7, 42))


ggsave(filename= "All_cohorts_GAweeks_histograms_binwidth0.5_no_lines.png",
       plot = All_hist_trim,
       device = "png",
       path = "../Results/Histograms_timepoints/",
       width = 20,
       height = 15,
       units = "cm")

#title <- text_grob("Timepoints", size = 15, face = "bold")

hist_plots <- grid.arrange(Oslo_hist_trim, Stanford_hist_trim, Detroit_hist_trim, All_hist_trim,
                           layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                           left = "No. Patients",
                           bottom = "Gestational Week")


ggsave(filename= "GAweeks_histograms_hlines_binwidth0.5_No_lines.png",
       plot = hist_plots,
       device = "png",
       path = "../Results/Histograms_timepoints/",
       width = 30,
       height = 20,
       units = "cm")

#Remove time point 1 and 5
Merged_samp2 <- Merged_samp %>%
  filter(GAWeeks >11.9) %>%
  filter(GAWeeks <33.9)

All_hist_trim2 <- ggplot(data = Merged_samp2, aes(x = GAWeeks, fill = Cohort, color = Cohort)) +
  geom_histogram(binwidth = 0.3, alpha = 0.7) +
  #ggtitle("All cohorts") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank()) +
  ylab("No. of participants") +
  xlab("Weeks of gestation") +
  #ylab("") +
  #xlab("") +
  geom_density() +
  #geom_vline(xintercept = 12, size = 1, color = "black") +
  geom_vline(xintercept = 19, size = 1, color = "black") +
  geom_vline(xintercept = 27, size = 1, color = "black") +
  #geom_vline(xintercept = 34, size = 1, color = "black") +
  coord_cartesian(xlim = c(10, 35))

ggsave(filename= "All_cohorts_GAweeks_histograms_binwidth0.5_T2T3T4_2.png",
       plot = All_hist_trim2,
       device = "png",
       path = "../Results/Histograms_timepoints/",
       width = 20,
       height = 10,
       units = "cm")

ggsave(filename= "Figure1_2.png",
       plot = All_hist_trim2,
       device = "png",
       path = "../Figures",
       width = 20,
       height = 10,
       units = "cm",
       dpi = 300)


# Add timepoints -----------------------------------------------
#Add timpoints based on GAweeks. T1: <12, T2: 12-<19, T3: 19-<27, T4: 27-<34, T5 >34 
load("../Data/Merged_cohorts_RFU_exprset.RData")
Merged_cohorts_RFU_exprset$Timepoint <- cut(Merged_cohorts_RFU_exprset$GAWeeks, c(-Inf,11.9,18.9,26.9,33.9,Inf), c("T1", "T2", "T3", "T4", "T5"))
Merged_cohorts_RFU_exprset$Timepoint <- as.character(Merged_cohorts_RFU_exprset$Timepoint)

samp_merged_cohorts <- pData(Merged_cohorts_RFU_exprset)
save(Merged_cohorts_RFU_exprset, file = "../Data/Merged_cohorts_RFU_exprset.RData")



# PCA plot ----------------------------------------------------------------
load("../Data/Merged_cohorts_RFU_exprset.RData")

pc1 <- pca(Merged_cohorts_RFU_exprset)
df1 <- merge(scores(pc1), pData(Merged_cohorts_RFU_exprset), by = 0)
pca1 <- ggplot(df1, aes(PC1, PC2, color = Cohort, shape = Group)) +
  geom_point(size = 2) +
  labs(title = "Merged cohorts",
       subtitle = "Log2RFU") + 
  xlab(paste0("PC1 (R2 = ",100*pc1@R2[1],"%)")) + 
  ylab(paste0("PC2 (R2 = ",100*pc1@R2[2],"%)")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  theme_classic() +
  theme(legend.position = "right",
        title = element_text(size = 9))

ggsave(filename= "PCA_merged_cohorts_log2RFU.png",
       plot = pca1,
       device = "png",
       path = "../Results/Merging_cohorts/",
       width = 10,
       height = 10,
       units = "cm")
