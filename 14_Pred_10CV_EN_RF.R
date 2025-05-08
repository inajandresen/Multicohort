rm(list=ls())
#
library(pROC)
#require(ROCR)
#require(caret)
#require(MASS)
library(ROCR)
library(caret)
library(MASS)
library(tidymodels) 
library(recipes)
library(limma)
library(randomForest)

modtype="EN" #RF #EN
tpoints=c("T2","T3","T4")

for(trim in tpoints){

load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
esetbig=exprs(Merged_cohorts_wins_Detrended_exprset)
anobig=pData(Merged_cohorts_wins_Detrended_exprset)

anosafe=anobig[anobig$Timepoint==trim,]
esetsafe=esetbig[,rownames(anosafe)]

all(rownames(anosafe)==colnames(esetsafe))
anosafe$G=factor(ifelse(anosafe$Group=="PE","D","C"))
#anosafe have G as "C", "D"


#split data in 10 folds (balanced by outcome and cohort)
anosafe$Cohort_G <- paste(anosafe$Cohort, anosafe$G, sep = "_")
set.seed(1);
anou=anosafe[!duplicated(anosafe$Cohort_ID),]
folds=createFolds(anou$Cohort_G,10)


lapply(folds,function(x){anou$Cohort_G[x]})

pile=list()



#parallel loop for training and testing
for(ite in 1:length(folds)){
  cat(ite);cat("\n")
  ano=anosafe
  
  esetn=esetsafe
  
  ano$subset_EOPE=ifelse(ano$Cohort_ID%in%anou$Cohort_ID[folds[[ite]]],"validation","discovery")
  anothis=ano
  
  ano_tr=ano[ano$subset_EOPE=="discovery",]
  esetn_tr=esetn[,rownames(ano_tr)]

  
  #####prepare test data 
  ano_t=anothis[anothis$subset_EOPE=="validation",]
  esetn_t=esetsafe[,rownames(ano_t)]
  
  if(TRUE){
  #select top genes by wilcoxon test
  ps=apply(esetn_tr,1,function(x){wilcox.test(x~ano_tr$G)$p.v})
  ps=sort(ps)
  sel=names(ps)[1:1006]
  }else{
  
  design <- model.matrix(~0+G,ano_tr) 
  colnames(design)<-gsub("G","",colnames(design))
  pe_v_control_cont = makeContrasts(
    D - C,
    levels = design
  )
  fit_noBayes = lmFit(esetn_tr, design)
  fit_pe_v_control = contrasts.fit(fit_noBayes, pe_v_control_cont)
  fit_pe_v_control = eBayes(fit_pe_v_control)
  deT1=topTable(fit_pe_v_control,coef=1, number = Inf)
  sel=rownames(deT1)[1:1006]
  }
  

  train=data.frame(G=factor(ano_tr$G),t(esetn_tr[sel,]))
  test=data.frame(G=factor(ano_t$G),t(esetn_t[sel,]))
  
  if(modtype=="EN"){
  ####select genes by lasso
  #EN
  p_recipe <- recipe(G ~ ., data = train) %>% 
    step_zv(all_numeric(), -all_outcomes()) 
  
  lasso <- logistic_reg(mixture = 0.5,penalty=0.5) %>% 
    set_engine("glmnet") %>% 
    set_mode("classification") 
  
  lasso_wf <- workflow() %>%
    add_recipe(p_recipe)
  
  lasso_fit <- lasso_wf %>%
    add_model(lasso) %>%
    fit(data = train)
  
  lasso_tun<- logistic_reg( mixture = tune(),penalty=tune()) %>% 
    set_engine("glmnet") %>% 
    set_mode("classification") 
  
  # Define the search grid for hyperparameter tuning 
  lasso_grid <- grid_regular( penalty(), mixture(),levels = 10) 
  
  # Tune the hyperparameters using split sample
  lasso_res <- tune_grid( lasso_wf %>% 
                            add_model(lasso_tun), resamples=apparent(train),
                          grid = lasso_grid,control=control_grid(parallel_over="everything"))
  
  
  # Get the best performing model 
  lowest_auc <- lasso_res %>%
    select_best("roc_auc")
  lasso_model_tun <- finalize_workflow(
    lasso_wf %>% add_model(lasso_tun),
    lowest_auc
  )%>%  fit(data = train)
  
  final<-lasso_model_tun%>%extract_fit_parsnip() %>% 
    tidy()
  
  #
  ###Predictions
  lasso_pred_tun <- as.matrix(predict(lasso_model_tun, test,type = "prob"))[,".pred_D"] 
  #final$term[final$estimate!=0][-1]
  
  }else{
  modb= randomForest(train[,-1], y=train$G,  ntree=500,importance=FALSE)
  lasso_pred_tun=predict(modb,test,type="prob")[,"D"]
  }
  
  outcs=as.numeric(test$G=="D")
  out=lasso_pred_tun
  outmat=cbind(outcs,out) #true outcomes for current test set side by side with predictited risk scores
  rownames(outmat) = rownames(test)
  pile[[ite]]<-list(outmat=outmat,feat=sel)
}


#create the full matrix of predicted risk scores and true outcomes
a=NULL
for (i in 1:length(pile)){
  a=rbind(a,pile[[i]]$outmat)
}


pdf(paste("../Prediction/Prediction_all_proteins/ROC_",trim,"_",modtype,".pdf",sep=""))
RC=roc(response=a[,1],predictor=a[,2],direction="<")
plot(1-RC$specificities,RC$sensitivities,col="black",type="l",lwd=2,xlab="False positive rate",
     ylab="Sensitivity")
abline(0,1)
aucsrf=round(ci.auc(RC),3)
#legend("bottomright",c(paste("Elastic Net AUC=",aucsrf[2],"(",aucsrf[1],"-",aucsrf[3],")",sep="")),lwd=2)

RC_dt=roc(response=a[grep("Detroit", rownames(a)),1],predictor=a[grep("Detroit", rownames(a)),2],direction="<")
points(1-RC_dt$specificities,RC_dt$sensitivities,col="#F8766D",type="l",lwd=2)
aucsrf_dt=round(ci.auc(RC_dt),3)

RC_std=roc(response=a[grep("Stanford", rownames(a)),1],predictor=a[grep("Stanford", rownames(a)),2],direction="<")
points(1-RC_std$specificities,RC_std$sensitivities,col="#619CFF",type="l",lwd=2)
aucsrf_std=round(ci.auc(RC_std),3)

RC_osl=roc(response=a[grep("Oslo", rownames(a)),1],predictor=a[grep("Oslo", rownames(a)),2],direction="<")
points(1-RC_osl$specificities,RC_osl$sensitivities,col="#00BA38",type="l",lwd=2)
aucsrf_osl=round(ci.auc(RC_osl),3)

legend <- c(paste("Overall AUC=",aucsrf[2],"(",aucsrf[1],"-",aucsrf[3],")",sep=""),
            paste("Detroit AUC=",aucsrf_dt[2],"(",aucsrf_dt[1],"-",aucsrf_dt[3],")",sep=""),
            paste("Stanford AUC=",aucsrf_std[2],"(",aucsrf_std[1],"-",aucsrf_std[3],")",sep=""),
            paste("Oslo AUC=",aucsrf_osl[2],"(",aucsrf_osl[1],"-",aucsrf_osl[3],")",sep=""))

legend("bottomright", legend = legend, col = c("black", "#F8766D", "#619CFF" ,"#00BA38"), lwd=2)
title(paste(modtype,trim))            

dev.off()


#Boxplot
a_mod <- a
ano_mod <- ano %>% dplyr::select(Cohort_G, Cohort)
ano_mod <- merge(a_mod, ano_mod, by = "row.names")

b1 <- ggplot(ano_mod, aes(x=Cohort_G, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#00BA38", "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted outcome") +
  ggtitle(paste(modtype,trim, sep=" ")) +
  theme(legend.title = element_blank())

#ggsave(filename= paste("Boxplots_all_proteins_",trim,"_",modtype,".png",sep=""),
 #      plot = b1,
#       device = "png",
#       path = "../Prediction/Prediction_all_proteins",
#       width = 20,
#      height = 10,
#      units = "cm")

###
#Cross validation stats

stat=table(unlist(lapply(pile,function(x){(x$feat)})))
stat=sort(stat,decreasing=TRUE)

freq=cbind(data.frame(stat))
#save(pile,freq,file=paste("../Prediction/Prediction_all_proteins/Predict_",trim,"",modtype,".RData",sep=""))


}



