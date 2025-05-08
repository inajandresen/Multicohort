rm(list=ls())
library(splines)
library(mgcv)
library(tidyverse)
library(Biobase)
library(readxl)
library(reshape2)
library(ggforce)
library(gridExtra)

#SL005034 RAN 
#SL007373 PPID 
#SL000572 SAA1
#SL005217 Siglec-6
#SL000525 MMp7


sig_prot <- c("SL005034", "SL007373", "SL000572", "SL005217", "SL000525")

# Extract protein data ----------------------------------------------------
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")

#Plot merged cohorts 
load("../Line_graphs_DE_proteins/Pred_LOPE_DE_proteins.RData")
Pred_df_LOPE_sig_proteins <- filter(Pred_LOPE_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
load("../Line_graphs_DE_proteins/Pred_Control_DE_proteins.RData")
Pred_df_control_sig_proteins <- filter(Pred_control_df, grepl(paste(sig_prot, collapse = "|"), SomaId))


prot <- fData(Merged_cohorts_wins_Detrended_exprset)
prot2 <- prot %>%
  mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
prot2$TargetFullName_EntrezGeneSymbol <- paste(prot2$TargetFullName, prot2$EntrezGeneSymbol2)
prot2 <- prot2 %>% 
  dplyr::select(SomaId, TargetFullName_EntrezGeneSymbol)
prot2 <- prot2[sig_prot , ]



All_df_merged_cohorts <- rbind(LOPE_df, Control_df)
All_df_merged_cohorts_long <- All_df_merged_cohorts %>%
  dplyr::select(-ID, -Group, -GAWeeks, -Cohort, -Y)
All_df_merged_cohorts_long <- melt(All_df_merged_cohorts_long, id = "SampleID")
colnames(All_df_merged_cohorts_long) <- c("SampleID", "SomaId", "log2MoM")

All_df_merged_cohorts_long <- merge(All_df_merged_cohorts_long, prot2, by = "SomaId")
All_df_samp_merged_cohorts <- All_df_merged_cohorts %>%
  dplyr::select(SampleID, ID, Group, GAWeeks, Cohort)
All_df_merged_cohorts_long <- merge(All_df_samp_merged_cohorts, All_df_merged_cohorts_long, by = "SampleID")


plot_merged_cohorts <- ggplot(data = All_df_merged_cohorts_long, aes(x = GAWeeks, y = log2MoM)) +
  geom_point(aes(shape = Group, color = Group), alpha = 0.6) +
  scale_color_manual(values = c("grey74", "grey40")) +
  facet_wrap(~TargetFullName_EntrezGeneSymbol, ncol = 1) +
  #facet_wrap(~TargetFullName_EntrezGeneSymbol, ncol = 1, scales = "free") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_line(data = Pred_df_control_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "Black", size=1.5) +
  geom_line(data = Pred_df_LOPE_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "solid", color = "Black", size=1.5) + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "nope",
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  xlab("Gestational weeks")

#LOPE: "grey40"
#Controls: "grey74", 




# Plot the cohorts seperately ---------------------------------------------
swr = function(string, nwrap=35) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

load("../Line_graphs_DE_proteins/Pred_Stanford_LOPE_DE_proteins.RData")
Pred_Stanford_LOPE_sig_proteins <- filter(Pred_Stanford_LOPE_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
#Pred_Stanford_LOPE_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(Pred_Stanford_LOPE_sig_proteins$TargetFullName_EntrezGeneSymbol)

load("../Line_graphs_DE_proteins/Pred_Stanford_control_DE_proteins.RData")
Pred_Stanford_control_sig_proteins <- filter(Pred_Stanford_control_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
#Pred_Stanford_control_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(Pred_Stanford_control_sig_proteins$TargetFullName_EntrezGeneSymbol)

load("../Line_graphs_DE_proteins/Pred_Oslo_LOPE_DE_proteins.RData")
Pred_Oslo_LOPE_sig_proteins <- filter(Pred_Oslo_LOPE_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
#Pred_Oslo_LOPE_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(Pred_Oslo_LOPE_sig_proteins$TargetFullName_EntrezGeneSymbol)

load("../Line_graphs_DE_proteins/Pred_Oslo_control_DE_proteins.RData")
Pred_Oslo_control_sig_proteins <- filter(Pred_Oslo_control_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
#Pred_Oslo_control_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(Pred_Oslo_control_sig_proteins$TargetFullName_EntrezGeneSymbol)

load("../Line_graphs_DE_proteins/Pred_Detroit_LOPE_DE_proteins.RData")
Pred_Detroit_LOPE_sig_proteins <- filter(Pred_Detroit_LOPE_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
#Pred_Detroit_LOPE_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(Pred_Detroit_LOPE_sig_proteins$TargetFullName_EntrezGeneSymbol)

load("../Line_graphs_DE_proteins/Pred_Detroit_control_DE_proteins.RData")
Pred_Detroit_control_sig_proteins <- filter(Pred_Detroit_control_df, grepl(paste(sig_prot, collapse = "|"), SomaId))
#Pred_Detroit_control_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(Pred_Detroit_control_sig_proteins$TargetFullName_EntrezGeneSymbol)


All_df <- rbind(Detroit_control_df, Detroit_LOPE_df)
All_df <- rbind(All_df, Oslo_control_df)
All_df <- rbind(All_df, Oslo_LOPE_df)
All_df <- rbind(All_df, Stanford_control_df)
All_df <- rbind(All_df, Stanford_LOPE_df)

All_df_long <- All_df %>%
  #dplyr::select(-ID, -Group, -GAWeeks, -Cohort, -Y)
  dplyr::select(-ID, -Group, -GAWeeks, -Cohort)
All_df_long <- melt(All_df_long, id = "SampleID")
colnames(All_df_long) <- c("SampleID", "SomaId", "log2MoM")

All_df_long <- merge(All_df_long, prot2, by = "SomaId")
All_df_samp <- All_df %>%
  dplyr::select(SampleID, ID, Group, GAWeeks, Cohort)
All_df_long <- merge(All_df_samp, All_df_long, by = "SampleID")
All_df_long_sig_proteins <- filter(All_df_long, grepl(paste(sig_prot, collapse = "|"), SomaId))
#All_df_long_sig_proteins$TargetFullName_EntrezGeneSymbol <- swr(All_df_long_sig_proteins$TargetFullName_EntrezGeneSymbol)

plot_sep_cohorts <- ggplot(data = All_df_long_sig_proteins, aes(x = GAWeeks, y = log2MoM)) +
  geom_point(aes(shape = Group, color = Group), alpha = 0.6) +
  scale_color_manual(values = c("grey74", "grey40")) +
  #facet_wrap(~TargetFullName_EntrezGeneSymbol, ncol = 1, scales = "free") +
  facet_wrap(~TargetFullName_EntrezGeneSymbol, ncol = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_line(data = Pred_Detroit_control_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#F8766D", size=1.5) +
  geom_line(data = Pred_Detroit_LOPE_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#F8766D", size=1.5) + 
  geom_line(data = Pred_Oslo_control_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#00BA38", size=1.5) +
  geom_line(data = Pred_Oslo_LOPE_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#00BA38", size=1.5) + 
  geom_line(data = Pred_Stanford_control_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#619CFF", size=1.5) +
  geom_line(data = Pred_Stanford_LOPE_sig_proteins, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#619CFF", size=1.5) + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  xlab("Gestational weeks")






# Make a collection of legends, seperate cohorts --------------------------------------------
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]
Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", Merged_cohorts_wins_Detrended_exprset$Group)
exprs <- exprs(Merged_cohorts_wins_Detrended_exprset)
prot <- fData(Merged_cohorts_wins_Detrended_exprset)
samp <- pData(Merged_cohorts_wins_Detrended_exprset)

#Format the expression matrix
exprs <- melt(exprs)
colnames(exprs) <- c("SomaId", "SampleID", "lgRFU")

#Format the samp data
samp <- samp %>%
  dplyr::select(Trimester, Timepoint, GAWeeks, Group, Cohort, Cohort_ID) %>%
  rownames_to_column(var = "SampleID")

#Merge expression data and samp
exprs <- merge(samp, exprs, by = "SampleID")

#Select the relevant columns from prot
prot <- prot %>%
  dplyr::select(SomaId, TargetFullName, AptName, EntrezGeneSymbol)
prot <- prot %>%
  mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
prot$TargetFullName_EntrezGeneSymbol <- paste(prot$TargetFullName, prot$EntrezGeneSymbol2)

#Merge expression data and prot
exprs <- merge(prot, exprs, by = "SomaId")



#Function for extracting legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#Make a legend for PE and Control, points
legend_PE_ctrl_plot_points <- ggplot(data =  filter(exprs, AptName == "MMP7.2789.26"), aes(x = GAWeeks, y = lgRFU)) +
  geom_point(aes(shape = Group, color = Group), alpha = 0.5) +
  scale_color_manual(values = c("grey70", "grey10")) +
  facet_wrap(~TargetFullName_EntrezGeneSymbol, ncol = 3, scales = "free") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15))

legend_PE_ctrl_points <- get_legend(legend_PE_ctrl_plot_points)

#Make a legend for PE and Control, lines
legend_PE_ctrl_plot_lines <- ggplot(data =  filter(exprs, AptName == "MMP7.2789.26"), aes(x = GAWeeks, y = lgRFU)) +
  geom_line(aes(linetype = Group), size = 1.4) +
  scale_linetype_manual(values=c("twodash", "solid"))+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15))

legend_PE_ctrl_lines <- get_legend(legend_PE_ctrl_plot_lines)

#Make a Detroit legend
exprs$Cohort_group <- paste(exprs$Cohort, exprs$Group, sep = " ") 
detroit_legend_plot <- ggplot(data =  filter(exprs, Cohort == "Detroit", AptName == "MMP7.2789.26"), aes(x = GAWeeks, y = lgRFU, group = Cohort_group)) +
  geom_line(aes(linetype = Cohort_group), color = "#F8766D", size = 1.4) +
  scale_linetype_manual(values=c("twodash", "solid"))+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15))

legend_detroit <- get_legend(detroit_legend_plot)

#Make a Oslo legend
Oslo_legend_plot <- ggplot(data =  filter(exprs, Cohort == "Oslo", AptName == "MMP7.2789.26"), aes(x = GAWeeks, y = lgRFU, group = Cohort_group)) +
  geom_line(aes(linetype = Cohort_group), color = "#00BA38", size = 1.4) +
  scale_linetype_manual(values=c("twodash", "solid"))+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15))

legend_oslo <- get_legend(Oslo_legend_plot)


#Make a stanford legend
stanford_legend_plot <- ggplot(data =  filter(exprs, Cohort == "Stanford", AptName == "MMP7.2789.26"), aes(x = GAWeeks, y = lgRFU, group = Cohort_group)) +
  geom_line(aes(linetype = Cohort_group), color = "#619CFF", size = 1.4) +
  scale_linetype_manual(values=c("twodash", "solid"))+
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=6))) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15))

legend_stanford <- get_legend(stanford_legend_plot)

#Arrange the legends
lay <- rbind(c(2, 3, 4, 5, 6))

legends_arrange <- arrangeGrob(legend_PE_ctrl_points, legend_PE_ctrl_lines, legend_detroit, legend_oslo, legend_stanford,
                               layout_matrix = lay)


# Arrange the plots -------------------------------------------------------
plots_arrange <- grid.arrange(plot_merged_cohorts, plot_sep_cohorts, ncol = 2)

#Combine the plots and the legend
Line_plots <- grid.arrange(plots_arrange, legends_arrange, ncol = 1,
                           heights = c(5, 0.2)) 

ggsave(filename= "Line_plots_RAN_SAA1_PPID_SIGLEG6_MMP7.png",
       plot = Line_plots,
       device = "png",
       path = "../Line_graphs_DE_proteins/",
       width = 27,
       height = 40,
       units = "cm")

ggsave(filename= "Figure5_2.png",
       plot = Line_plots,
       device = "png",
       path = "../Figures",
       width = 27,
       height = 40,
       units = "cm",
       dpi = 350)


