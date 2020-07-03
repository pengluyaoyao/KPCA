##############################################################
###############      KPCA CODE      ##########################
##############################################################
kernel<-function(x,y,sigma){
   diff<-x-y
   k<- exp(-t(diff)%*%diff/(2*sigma^2))
}

recerr<-function(test,data,sigma,sumalpha,alpha.norm,alphaKrow,Ksum){
   n<-nrow(data)
   k<-matrix(0,n,1)  
   for (i in 1:n){
         k[i]<-kernel(test, data[i,],sigma)
      } 
   f<-t(k)%*%alpha.norm-sumalpha*(sum(k)/n)+sumalpha*Ksum-alphaKrow
   err<-kernel(test,test,sigma)-2*sum(k)/n+Ksum-f%*%t(f)
   return(err)
}

kpcabound<-function(data,scale.test,sigma,numev,outlier){
   n<-nrow(data)
   K<-matrix(0,n,n)

   for (i in 1:n){
      for (j in 1:n){
         K[i,j]<-kernel(data[i,],data[j,],sigma)
      }
   }

   Krow<-colSums(K)/n
   Ksum<-sum(Krow)/n

   for (i in 1:n){
      for (j in 1:n){
          K[i,j]<-K[i,j]-Krow[i]-Krow[j]+Ksum
      }
   }

   eigen<-eigen(K)
   lambda<-eigen$values
   alpha<-eigen$vectors
   alpha.norm<-matrix(0,n,numev)

   for (i in 1:numev){
      alpha.norm[,i]<-alpha[,i]/sqrt(lambda[i])
   }

   sumalpha = colSums(alpha.norm)
   alphaKrow<-Krow%*%alpha.norm

   err<-matrix(0,n,1)
      for (i in 1:n){
          test<-data[i,]
          err[i]<-recerr(test,data,sigma,sumalpha,alpha.norm,alphaKrow,Ksum)
    }

   serr<-err[order(err),]
   maxerr<-serr[n-outlier]
   
   m<-nrow(scale.test)/2
   err<-matrix(0,nrow(scale.test),1)
   ind<-matrix(c(rep(1,m),rep(0,m)),nrow(scale.test),1)
   tp<-0 
   fp<-0 
   for (i in 1:nrow(scale.test)){
         test0<-scale.test[i,]
         err[i]<-recerr(test0,data,sigma,sumalpha,alpha.norm,alphaKrow,Ksum)
         if(err[i]> maxerr && ind[i]==0){
           tp=tp+1
         } else if (err[i]> maxerr && ind[i]==1) {
           fp=fp+1
         } else {
           tp=tp
           fp=fp
         }
    }
   return(list(serr=serr,maxerr=maxerr,tp=tp,fp=fp,err=err))
}



sigma<-c(0.09,0.1,0.12,0.14,0.16)
evec<-c(20,25)
ptm<-proc.time()
for (s in 1:length(sigma)){
   for (e in 1:length(evec)){
      result<-kpcabound(scale.train,sigma[s],evec[e],0)
      print(c(result$fp,result$tp,result$maxerr))
  }
}
proc.time()-ptm



##########################################################################
##########################################################################
setwd('C:/Users/luyao_peng/Desktop/kpca_file/data')
 raw.narr<-read.csv('narrative.csv')
 narr<-raw.narr[,c(-1,-2)]
 regular<-narr[which(narr[,1]==1),]
 copy<-narr[which(narr[,1]==0),]
 set.seed(12345)
 test.reg.ind<-sample(1:nrow(regular),nrow(copy))
 test.reg<-regular[test.reg.ind,-1]
 test<-rbind(test.reg,copy[,-1])

 left.reg<-regular[-test.reg.ind,-1]
 set.seed(333)
 train.ind<-sample(1:nrow(left.reg),3000)
 train<-left.reg[train.ind,]
 
 
 train<-as.matrix(train)
 test<-as.matrix(test)



library(matrixStats)
library(ROCR)
###################################Informative 
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 info1<-read.csv("informative1.csv")
 info0<-read.csv("informative0.csv")
 n<-nrow(info1)
 m<-nrow(info0)
 train.info1<-info1[1:3000,]
 left.info1.ind<-seq(1,n)[-c(1:3000)]
 left.info1<-info1[left.info1.ind,]
 set.seed(123)
 test.info1<-left.info1[sample(1:nrow(left.info1),m),]
 test<-rbind(test.info1,info0)
 train<-as.matrix(train.info1)
 test<-as.matrix(test)

 sds<-colSds(train)
 means<-colMeans(train)

 scale.test<-matrix(0,nrow(test),ncol(test))
 for(i in 1:ncol(test)){
    scale.test[,i]<-(test[,i]-means[i])/sds[i]
 }


 scale.train<-scale(train)

###################Argumentative Z dendro sampling
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 argu1<-read.csv('argumentative1.csv')
 argu0<-read.csv('argumentative0.csv')
 test.novel<-argu0
 index<-seq(1,nrow(argu1),1)
 argu1.w.ind<-cbind(index,argu1)
 set.seed(1234555)
 test.reg.ind<-sample(1: nrow(argu1.w.ind),nrow(argu0))
 test.reg<-argu1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-argu1.w.ind[-test.reg.ind,]
 scale.left.reg<-scale(left.reg[,-1])
 clusters<-hclust(dist(scale.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=2000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,2000,1)
 train<-matrix(0,2000,ncol(argu1))
 for (i in 1:2000){
    cand<-cut.ind[which(cut.ind[,2]==i),]    
    cand<-as.matrix(cand)
    set.seed(1234444)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(argu1[train.ind[i],])    
 }
 
 
 sds<-colSds(train)
 means<-colMeans(train)
 
 scale.test<-matrix(0,nrow(test),ncol(test))
 for(i in 1:ncol(test)){
    scale.test[,i]<-(test[,i]-means[i])/sds[i]
 }

 scale.train<-scale(train)
 scale.test<-as.matrix(scale.test)
 scale.train<-as.matrix(scale.train)
####################################Scaled Narrative in Matlab##################
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/data')
 test.normal<-read.csv("test_normal.csv")
 test.novel<-read.csv("test_novel.csv")
 train.normal<-read.csv("train_normal.csv")
 test<-rbind(test.normal,test.novel)
 scale.test<-as.matrix(test)
 scale.train<-as.matrix(train.normal)
  m<-nrow(test.novel)



#################dendrogram sampling for narrative Z transform before dendrogram sampling############################
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 narr1<-read.csv('narrative1.csv')
 narr0<-read.csv('narrative0.csv')
 test.novel<-narr0
 index<-seq(1,nrow(narr1),1)
 narr1.w.ind<-cbind(index,narr1)
 set.seed(12345)
 test.reg.ind<-sample(1: nrow(narr1.w.ind),nrow(narr0))
 test.reg<-narr1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-narr1.w.ind[-test.reg.ind,]
 scale.left.reg<-scale(left.reg[,-1])
 clusters<-hclust(dist(scale.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=3000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,3000,1)
 train<-matrix(0,3000,ncol(narr1))
 for (i in 1:3000){
    cand<-cut.ind[which(cut.ind[,2]==i),]    
    cand<-as.matrix(cand)
    set.seed(12341)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(narr1[train.ind[i],])    
 }
 
 
 sds<-colSds(train)
 means<-colMeans(train)
 
 scale.test<-matrix(0,nrow(test),ncol(test))
 for(i in 1:ncol(test)){
    scale.test[,i]<-(test[,i]-means[i])/sds[i]
 }

 scale.train<-scale(train)
###################01 scale
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 narr1<-read.csv('narrative1.csv')
 narr0<-read.csv('narrative0.csv')
 test.novel<-narr0
 index<-seq(1,nrow(narr1),1)
 narr1.w.ind<-cbind(index,narr1)
 set.seed(12345)
 test.reg.ind<-sample(1: nrow(narr1.w.ind),nrow(narr0))
 test.reg<-narr1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-as.matrix(narr1.w.ind[-test.reg.ind,])
 
 unit.left.reg<-matrix(0,nrow(left.reg),14)
 mins<-colMins(left.reg[,-1])
 for (i in 1:14){
    unit.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
 }

 clusters<-hclust(dist(unit.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=3000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,3000,1)
 train<-matrix(0,3000,ncol(narr1))
 for (i in 1:3000){
    cand<-cut.ind[which(cut.ind[,2]==i),]    
    cand<-as.matrix(cand)
    set.seed(12341)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(narr1[train.ind[i],])    
 }

 unit.train<-matrix(0,3000,14)
 mins2<-colMins(train)
 maxs2<-matrix(0,1,14)
 for (i in 1:14){
    maxs2[i]<-max(train[,i]-mins2[i])
    unit.train[,i]<-(train[,i]-mins2[i])/maxs2[i]
 }
 
 unit.test<-matrix(0,368,14)
 for (i in 1:14){
     unit.test[,i]<-(test[,i]-mins2[i])/maxs2[i]
 }

 scale.test<-unit.test
 scale.train<-unit.train


#################dendrogram sampling for informative Z Transform before dendrogram sampling##############################
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 info1<-read.csv('informative1.csv')
 info0<-read.csv('informative0.csv')
 test.novel<-info0
 index<-seq(1,nrow(info1),1)
 info1.w.ind<-cbind(index,info1)
 set.seed(123498)
 test.reg.ind<-sample(1: nrow(info1.w.ind),nrow(info0))
 test.reg<-info1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-info1.w.ind[-test.reg.ind,]
 scale.left.reg<-scale(left.reg[,-1])
 clusters<-hclust(dist(scale.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=3000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,3000,1)
 train<-matrix(0,3000,ncol(info1))
 for (i in 1:3000){
    cand<-cut.ind[which(cut.ind[,2]==i),]
    set.seed(1234198)
    cand<-as.matrix(cand)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(info1[train.ind[i],])    
 }

 sds<-colSds(train)
 means<-colMeans(train)
 
 scale.test<-matrix(0,nrow(test),ncol(test))
 for(i in 1:ncol(test)){
    scale.test[,i]<-(test[,i]-means[i])/sds[i]
 }


 scale.train<-scale(train) 

#######01 scale for info dendro sampling 3000
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 info1<-read.csv('informative1.csv')
 info0<-read.csv('informative0.csv')
 test.novel<-info0
 index<-seq(1,nrow(info1),1)
 info1.w.ind<-cbind(index,info1)
 set.seed(123498)
 test.reg.ind<-sample(1: nrow(info1.w.ind),nrow(info0))
 test.reg<-info1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-as.matrix(info1.w.ind[-test.reg.ind,])
 
 unit.left.reg<-matrix(0,nrow(left.reg),14)
 mins<-colMins(left.reg[,-1])
 for (i in 1:14){
    unit.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
 }

 clusters<-hclust(dist(unit.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=3000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,3000,1)
 train<-matrix(0,3000,ncol(info1))
 for (i in 1:3000){
    cand<-cut.ind[which(cut.ind[,2]==i),]    
    cand<-as.matrix(cand)
    set.seed(1234198)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(info1[train.ind[i],])    
 }

 unit.train<-matrix(0,3000,14)
 mins2<-colMins(train)
 maxs2<-matrix(0,1,14)
 for (i in 1:14){
    maxs2[i]<-max(train[,i]-mins2[i])
    unit.train[,i]<-(train[,i]-mins2[i])/maxs2[i]
 }
 
 unit.test<-matrix(0,936,14)
 for (i in 1:14){
     unit.test[,i]<-(test[,i]-mins2[i])/maxs2[i]
 }

 scale.test<-unit.test
 scale.train<-unit.train

#################01 scale informative dendro sampling 4000
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 info1<-read.csv('informative1.csv')
 info0<-read.csv('informative0.csv')
 test.novel<-info0
 index<-seq(1,nrow(info1),1)
 info1.w.ind<-cbind(index,info1)
 set.seed(123498)
 test.reg.ind<-sample(1: nrow(info1.w.ind),nrow(info0))
 test.reg<-info1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-as.matrix(info1.w.ind[-test.reg.ind,])
 
scale.left.reg<-matrix(0,nrow(left.reg),14)
 mins<-colMins(left.reg[,-1])
 for (i in 1:14){
    scale.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
 }

 clusters<-hclust(dist(scale.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=4000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,4000,1)
 train<-matrix(0,4000,ncol(info1))
 for (i in 1:4000){
    cand<-cut.ind[which(cut.ind[,2]==i),]    
    cand<-as.matrix(cand)
    set.seed(1234198)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(info1[train.ind[i],])    
 }

 unit.train<-matrix(0,4000,14)
 mins2<-colMins(train)
 maxs2<-matrix(0,1,14)
 for (i in 1:14){
    maxs2[i]<-max(train[,i]-mins2[i])
    unit.train[,i]<-(train[,i]-mins2[i])/maxs2[i]
 }
 
 unit.test<-matrix(0,936,14)
 for (i in 1:14){
     unit.test[,i]<-(test[,i]-mins2[i])/maxs2[i]
 }

 scale.test<-unit.test
 scale.train<-unit.train


#########01 scale argumentative
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 argu1<-read.csv('argumentative1.csv')
 argu0<-read.csv('argumentative0.csv')
 test.novel<-argu0
 index<-seq(1,nrow(argu1),1)
 argu1.w.ind<-cbind(index,argu1)
 set.seed(1234555)
 test.reg.ind<-sample(1: nrow(argu1.w.ind),nrow(argu0))
 test.reg<-argu1.w.ind[test.reg.ind,]
 test<-rbind(test.reg[,-1],test.novel)
 left.reg<-as.matrix(argu1.w.ind[-test.reg.ind,])

 unit.left.reg<-matrix(0,nrow(left.reg),14)
 mins<-colMins(left.reg[,-1])
 for (i in 1:14){
    unit.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
 }

 clusters<-hclust(dist(unit.left.reg),method="complete")
 ##plot(clusters,hang=-1,xlab="Columns",ylab="Distance")
 cut<-cutree(clusters,k=2000)
 cut<-as.matrix(cut)
 cut.ind<-cbind(left.reg[,1],cut)
 train.ind<-matrix(0,2000,1)
 train<-matrix(0,2000,ncol(info1))
 for (i in 1:2000){
    cand<-cut.ind[which(cut.ind[,2]==i),]    
    cand<-as.matrix(cand)
    set.seed(1234444)
    train.ind[i]<-sample(cand[,1],1)
    train[i,]<-t(argu1[train.ind[i],])    
 }

 unit.train<-matrix(0,2000,14)
 mins2<-colMins(train)
 maxs2<-matrix(0,1,14)
 for (i in 1:14){
    maxs2[i]<-max(train[,i]-mins2[i])
    unit.train[,i]<-(train[,i]-mins2[i])/maxs2[i]
 }
 
 unit.test<-matrix(0,nrow(test),ncol(test))
 for (i in 1:ncol(test)){
     unit.test[,i]<-(test[,i]-mins2[i])/maxs2[i]
 }

 scale.test<-unit.test
 scale.train<-unit.train

##########
dendro_group<-1
train2<-cbind(dendro_group,narr1[1,])
  for (i in 2:100){
    cand<-cut.ind[which(cut.ind[,2]==i),1]    
    cand<-matrix(cand)
    dendro_group<-rep(i,nrow(cand))
    train<-cbind(dendro_group,narr1[cand[,1],]) 
    train2<-rbind(train2,train)   
 }


