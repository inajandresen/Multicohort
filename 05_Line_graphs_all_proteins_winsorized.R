rm(list=ls())
library(tidyverse)
library(ggforce)
library(reshape2)
library(readxl)
library(Biobase)
library(scales)
library(gridExtra)
#library(grid)
library(ggplotify)

load("../Data/Merged_cohorts_wins_lgRFU_exprset.RData")
exprs <- exprs(Merged_cohorts_wins_lgRFU_exprset)
prot <- fData(Merged_cohorts_wins_lgRFU_exprset)
samp <- pData(Merged_cohorts_wins_lgRFU_exprset)

#Format the expression matrix
exprs <- melt(exprs)
colnames(exprs) <- c("SomaId", "SampleID", "lgRFU")

#Format the samp data
samp <- samp %>%
  select(Trimester, Timepoint, GAWeeks, Group, Cohort, Cohort_ID) %>%
  rownames_to_column(var = "SampleID")

#Merge expression data and samp
exprs <- merge(samp, exprs, by = "SampleID")

#Select the relevant columns from prot
prot <- prot %>%
  select(SomaId, TargetFullName, AptName, EntrezGeneSymbol)
prot <- prot %>%
  mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
prot$TargetFullName_EntrezGeneSymbol <- paste(prot$TargetFullName, prot$EntrezGeneSymbol2)

#Split over several lines if header is longer than 50 characters
swr = function(string, nwrap=45) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

prot$TargetFullName_EntrezGeneSymbol <- swr(prot$TargetFullName_EntrezGeneSymbol)


#Merge expression data and prot
exprs <- merge(prot, exprs, by = "SomaId")



#Plot Detroit
Detroit_exprs <- exprs %>%
  filter(Cohort == "Detroit")


pdf("../Line_plots_win/Lineplots_all_proteins_Detroit_wins.pdf",
    width = 10, height = 10)

for(i in 1:84) {    
  print(ggplot(data=Detroit_exprs, aes(x = GAWeeks, y = lgRFU, group = Cohort_ID, color = Group)) +
          geom_line() +
          geom_point() +
          scale_color_manual(values = c("grey74", "grey40")) +
          theme_classic() +
          xlab("Gestational weeks") +
          theme(legend.title = element_blank()) +
          facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 3, nrow = 4, scales = "free", page = i))
}

dev.off()

#Plot Oslo
Oslo_exprs <- exprs %>%
  filter(Cohort == "Oslo")

pdf("../Line_plots_win/Lineplots_all_proteins_Oslo_wins.pdf",
    width = 10, height = 10)

for(i in 1:84) {    
  print(ggplot(data=Oslo_exprs, aes(x = GAWeeks, y = lgRFU, group = Cohort_ID, color = Group)) +
          geom_line() +
          geom_point() +
          scale_color_manual(values = c("grey74", "grey40")) +
          theme_classic() +
          xlab("Gestational weeks") +
          theme(legend.title = element_blank()) +
          facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 3, nrow = 4, scales = "free", page = i))
}

dev.off()

#Plot Stanford
Stanford_exprs <- exprs %>%
  filter(Cohort == "Stanford")

pdf("../Line_plots_win/Lineplots_all_proteins_Stanford_wins.pdf",
    width = 10, height = 10)

for(i in 1:84) {    
  print(ggplot(data=Stanford_exprs, aes(x = GAWeeks, y = lgRFU, group = Cohort_ID, color = Group)) +
          geom_line() +
          geom_point() +
          scale_color_manual(values = c("grey74", "grey40")) +
          theme_classic() +
          xlab("Gestational weeks") +
          theme(legend.title = element_blank()) +
          facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 3, nrow = 4, scales = "free", page = i))
}

dev.off()


