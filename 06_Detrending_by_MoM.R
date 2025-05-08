library(splines)
library(mgcv)
library(Biobase)

load("../Data/Merged_cohorts_wins_lgRFU_exprset.RData")
bigeset=exprs(Merged_cohorts_wins_lgRFU_exprset)
bigano=pData(Merged_cohorts_wins_lgRFU_exprset)

for (coh in unique(bigano$Cohort)) {
  
  tmp=bigano[bigano$Cohort==coh,]
  eset=bigeset[,rownames(tmp)]
  tmp$GA=tmp$GAWeeks
  
  
  esetd=eset
  ans=rownames(eset)
  for( i in 1:nrow(eset)){
    an=ans[i]
    tmp$Y=eset[an,]
    #plot(Y~GA,data=tmp[tmp$Group=="Control",],pch=19,cex=0.6)
    tmp$group=factor(tmp$Group)
    
    #fit<-lm(Y ~ bs(GA,degree=2),data = tmp[tmp$Group=="Control",])
    tmpc<- tmp[tmp$Group=="Control",]
    tmpc$ID<-factor(tmpc$ID)
    
    
    
    fit3 <- gam(Y~s(GA,k=3)+s(ID,bs = 're'),
                data=tmpc,
                method="ML")
    
    
    fit4 <- gam(Y~s(GA,k=4)+s(ID,bs = 're'),
                data=tmpc,
                method="ML")
    
    fit5 <- gam(Y~s(GA,k=5)+s(ID,bs = 're'),
                data=tmpc,
                method="ML")
    
    fit6 <- gam(Y~s(GA,k=6)+s(ID,bs = 're'),
                data=tmpc,
                method="ML")
    
    fit=list(fit3,fit4,fit5,fit6)[[which.max(c(summary(fit3)$r.sq,
                                               summary(fit4)$r.sq,
                                               summary(fit5)$r.sq, 
                                               summary(fit6)$r.sq))]]
    
    # gam.check(fit)
    # ds=expand.grid(GA=seq(4,43,by=0.1))# 
    # ds$pred=predict.gam(fit,ds,exclude="s(ID)",newdata.guaranteed=TRUE)
    
    #points(pred~GA,data=ds,type="l",lwd=3,col="blue",  ylab="",xlab="") 
    tmp1<-tmp
    tmp1$ID<-NULL
    tmp$Ypred=predict.gam(fit,tmp1,exclude="s(ID)",newdata.guaranteed=TRUE)
    tmp$Yd=tmp$Y-tmp$Ypred
    #plot(Yd~GA,data=tmp[tmp$Group=="Control",],pch=19,cex=0.6)
    bigeset[an,rownames(tmp)] = tmp$Yd
  }
  
}

exprs(Merged_cohorts_wins_lgRFU_exprset)<-bigeset
boxplot(bigeset[219,]~bigano$Group+bigano$Cohort)
Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_lgRFU_exprset
save(Merged_cohorts_wins_Detrended_exprset,file="Merged_cohorts_wins_Detrended_exprset.RData")



