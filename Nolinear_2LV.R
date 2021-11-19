setwd("D:/data/R/STA/STA-master/STACMR-R")
source("staCMRsetup.R")
#在线性模型中，交互作用的H1成立时，通过数据模拟的方法，比较STA和ANOVA的I型错误
#假设是一个2（A）×2(B)×5(C)的完全随机设计。假设方差齐性
#——————————————————第一步，生成数据——————————————————————————
anova_vs_sta=function(se,n){
  mu=5#总平均
  #se=5#误差标准差,设置为5的话通过GPOWER计算出power为0.88左右，如果太大的话power太高，不便于比较
  #n=50#各处理样本量
  a=c(-1,1)#A因素各效应值
  b=c(-0.5,0.5)#B因素各效应值
  c=c(-2,-1,0,1,2)#C因素各效应值
  #ab=array(c(-0.5,0.5,0.5,-0.5),dim=c(2,2))
  na=n*length(b)*length(c)#A各水平的样本量
  nb=n*length(c)#A各水平下，B各水平的样本量
  m=length(a)*length(b)*length(c)#处理数
  idx=1
  error_R=rnorm((m*n),mean = 0,sd=se)#随机抽取m×n的随机误差
  Y=NULL
  for (i in 1:length(a)) {#生成观测值
    for (j in 1:length(b)) {
      for (k in 1:length(c)) {
        for (l in 1:n) {
          LV1=b[j]+c[k]
          LV2=b[j]-c[k]
          if(i==1)Y[idx]=a[i]+LV1^3+error_R[idx]
          if(i==2)Y[idx]=a[i]+2^LV2+error_R[idx]
          idx=idx+1
        }
        
      }
      
    }
    
  }
  #…………………………ANOVA分析…………………………
  anovaDF=data.frame(A=c(rep(1,na),rep(2,na)),
                     B=c(rep(1,nb),rep(2,nb),rep(1,nb),rep(2,nb)),
                     C=rep(c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n)),4),
                     anovaY=Y)
  anovaResults=aov(anovaY~factor(A)*factor(B)*factor(C),anovaDF)
  Anova=summary(anovaResults)
  if((Anova[[1]]$`Pr(>F)`[4]<=0.05)||(Anova[[1]]$`Pr(>F)`[5]<=0.05)||(Anova[[1]]$`Pr(>F)`[7]<=0.05)){nANOVA=1} else {nANOVA=0}
  #………………sta分析…………………………
  staDF=data.frame(SNo=c(1:(n*m)),#被试编号，不进入计算
                    Cno=rep(rep(c(1:10),each=n),2),#B因素和C因素的共10个处理，其中1-5表示B1水平下C1-C5；6-10表示B2水平下C1-C5
                    DV=rep(c(1,2),each=na),
                    staY=Y)
  #staPLOT(staDF,groups=list(c(1:5),c(6:10)),grouplabels=list('B1','B2'),
   #       axislabels=list('A1','A2'))
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
for (Verror in 1:10) {
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
  for (n in 1:1000) {
    nalpha=anova_vs_sta(n=50,se=Verror)
    nalpha_ANOVA=nalpha_ANOVA+nalpha[1]
    nalpha_STA=nalpha_STA+nalpha[2]
    resultsDF[n,]=c(n,nalpha)
    cat('进行第',n,'次取样','ANOVA拒绝H0次数为：',nalpha_ANOVA,'STA拒绝H0次数为：',nalpha_STA,'\n')
  }
  filenames=paste0('D:/data/R/STA/simulationData/NOLinear_2LV_E',Verror,'_N50.csv')
  write.csv(resultsDF,file = filenames) 
}

#cat('A与B交互作用在1000次模拟中，犯α错误的概率为:',np/1000)
