setwd("D:/data/R/STA/STA-master/STACMR-R")
source("staCMRsetup.R")
#在线性模型中，交互作用的H1成立时，通过数据模拟的方法，比较STA和ANOVA的I型错误
#假设是一个2（A）×2(B)×5(C)的完全随机设计。假设方差齐性
#—————————————————step 1: Data generation——————————————————————————
anova_vs_sta=function(se,n){# function of generation data, given se and n
  mu=5#Overall mean
  a=c(-1,1)#effects of factor A
  b=c(-0.5,0.5)#effects of factor B
  c=c(-2,-1,0,1,2)#effects of factor C
  na=n*length(b)*length(c)#Sample size in each level of factor A
  nb=n*length(c)#sample size in each A×B conditions
  m=length(a)*length(b)*length(c)#total number of conditions
  idx=1#NO. of observed values
  error_R=rnorm((m*n),mean = 0,sd=se)#Random errors in each observed value
  Y=NULL
  for (i in 1:length(a)) {#generated the data
    for (j in 1:length(b)) {
      for (k in 1:length(c)) {
        for (l in 1:n) {
          LV=b[j]+c[k]#the single latent variable
          if(i==1)Y[idx]=LV^3+error_R[idx]#observed values in a1 level
          if(i==2)Y[idx]=2^LV+error_R[idx]#observed values in a2 level
          idx=idx+1
        }
        
      }
      
    }
    
  }
  #…………………………ANOVA…………………………
  anovaDF=data.frame(A=c(rep(1,na),rep(2,na)),
                     B=c(rep(1,nb),rep(2,nb),rep(1,nb),rep(2,nb)),
                     C=rep(c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n)),4),
                     anovaY=Y)
  anovaResults=aov(anovaY~factor(A)*factor(B)*factor(C),anovaDF)
  Anova=summary(anovaResults)
  if((Anova[[1]]$`Pr(>F)`[4]<=0.05)||(Anova[[1]]$`Pr(>F)`[5]<=0.05)||(Anova[[1]]$`Pr(>F)`[7]<=0.05)){nANOVA=1} else {nANOVA=0}
  #…………………………sta…………………………
  staDF=data.frame(SNo=c(1:(n*m)),#被试编号，不进入计算
                    Cno=rep(rep(c(1:10),each=n),2),#B因素和C因素的共10个处理，其中1-5表示B1水平下C1-C5；6-10表示B2水平下C1-C5
                    DV=rep(c(1,2),each=na),
                    staY=Y)
  partialList=list(c(1:5),c(6:10),c(1,6),c(2,7),c(3,8),c(4,9),c(5,10))
  MR_results=staMR(staDF,partial=partialList)#偏序模型拟合值GOF
  CMR_results=staCMR(staDF,partial=partialList)#单调模型拟合值GOF
  staResult=staCMRFIT(staDF,partial = partialList,nsample = 10000)#单调回归检验
  truefit=CMR_results$fval-MR_results$fval
  truep=length(which(staResult$fits>truefit))/10000
  if(staResult$p<=0.05){nSTA=1}else{nSTA=0}
  return(c(nANOVA,nSTA,Anova[[1]]$`Pr(>F)`[4],Anova[[1]]$`Pr(>F)`[5],
           Anova[[1]]$`Pr(>F)`[6],Anova[[1]]$`Pr(>F)`[7],staResult$p,
           MR_fit=MR_results$fval,CMR_fit=CMR_results$fval,delta_fit=staResult$datafit,
           truep))

}
for (Verror in 1:10) {# standard deviation of random errors
  nalpha_ANOVA=0
  nalpha_STA=0
  resultsDF=data.frame(niter=0,#模拟次数
                       nANOVA=0,#ANOVAI型错误数
                       nSTA=0,#STAI型错误数
                       p_AinterB=0,#AB交互作用p值
                       p_AinterC=0,#AC交互作用p值
                       p_BinterC=0,#BC交互作用p值,此时没有用
                       p_AinterBinterC=0,#ABC交互作用p值
                       p_STA=0,
                       MR_fit=0,
                       CMR_fit=0,
                       delta_fit=0,
                       trueP_STA=0
  )
  for (n in 1:1000) {# analysis 1000 times
    nalpha=anova_vs_sta(n=10,se=Verror)#n could be replaced by 20,30,40,50
    nalpha_ANOVA=nalpha_ANOVA+nalpha[1]
    nalpha_STA=nalpha_STA+nalpha[2]
    resultsDF[n,]=c(n,nalpha)
    cat('进行第',n,'次取样','ANOVA拒绝H0次数为：',nalpha_ANOVA,'STA拒绝H0次数为：',nalpha_STA,'\n')
  }
  filenames=paste0('D:/data/R/STA/simulationData/NOLinear_1LV_E',Verror,'_N10.csv')
  write.csv(resultsDF,file = filenames) 
}

