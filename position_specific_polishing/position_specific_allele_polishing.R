library("fitdistrplus")

Position_based_allele_polishing <- function(i,indexTrain,x,g){
  indexTest=i
  vTest=x$Reads2[indexTest]/(x$Reads1[indexTest]+x$Reads2[indexTest])
  vTrain=unlist(V[indexTrain[i],]) 
  vTrain=as.numeric(vTrain[!is.na(vTrain)])
  p=NA
  if(length(vTrain)==0){
    p=NA
  }
  if(length(vTrain)>5){
    nonZero=length(vTrain)/g[i] 
    sh=as.numeric(fitdist(vTrain, distr = "weibull",method = "mle", lower = c(0, 0))[[1]][1])
    sc=as.numeric(fitdist(vTrain, distr = "weibull",method = "mle", lower = c(0, 0))[[1]][2])
    P = pweibull(vTest,shape = sh,scale = sc)
    p = 1-((1-nonZero)+(nonZero*P))
  }
  if(length(vTrain)<=5 & length(vTrain)!=0 ){
    mu=fitdistr(vTrain,"normal")$estimate[[1]]
    sd=fitdistr(vTrain,"normal")$estimate[[2]]
    a <- vTest
    s <- sd
    n <- length(vTrain)
    xbar <- mu
    p=1-pnorm(a,mean=xbar,sd=s)
  }
  return(p)
}

#Load position_specific_allele_polishing.Rdata
#In this example we will assign pvalues to determine whether variants allele frequencies (in dataframe x) 
#reach statistical significance when compared to their corresponding error distribution obtained from unrelated samples (in dataframe VAF)
#dataframe x was derived by using Varscan2
#Pvalue of NA means that a model cannot be determined due to the limited number of errors 

x$Pvalue_polishing=NA
indexTrain=match(paste(x$Chrom,x$Position,x$Ref,x$VarAllele),paste(VAF$Chrom,VAF$Position,VAF$Ref,VAF$VarAllele))
V=VAF[indexTrain,]
g=as.numeric(apply(V[,5:dim(V)[2]],
                   1, function(x) length(which(is.na(x))))) 
index=which(V[,5:dim(V)[2]]>=0.05,arr.ind = T) 
index[,2]=index[,2]+4
V[index]=NA
V=V[,5:dim(V)[2]]
g=g+as.numeric(apply(V, 1, function(x) length(which(!is.na(x)))))
indexTrain=1:dim(V)[1]
for(i in 1:dim(x)[1]){ 
  print(i)
  x$Pvalue_polishing[i]=Position_based_allele_polishing(i,indexTrain,x,g)
}

#PLEASE NOTICED THAT SMALL NUMBER OF CONTROL SAMPLES MIGHT NOT REPRESENT THE TRUE ERROR RATE. 
#IMPROVED RESULTS CAN BE ACHIEVED WHEN THE WEIBULL MODEL IS INITIATED, UTILIZING LARGER CONTROL DATASETS

