##########data preperation
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 narr1<-read.csv('narrative1.csv')
 narr0<-read.csv('narrative0.csv')
 test.novel<-narr0

 data.narr<-function(data1,data0,sampling,scale){
     narr1<-data1
     narr0<-data0
     index<-seq(1,nrow(narr1),1)
     narr1.w.ind<-cbind(index,narr1)
     set.seed(12345)
     test.reg.ind<-sample(1: nrow(narr1.w.ind),nrow(narr0))
     test.reg<-narr1.w.ind[test.reg.ind,]
     test<-rbind(test.reg[,-1],test.novel)
     left.reg<-narr1.w.ind[-test.reg.ind,]
     if (scale=="Z"){
        scale.left.reg<-scale(left.reg[,-1])
     } else {
        left.reg<-as.matrix(left.reg)
        scale.left.reg<-matrix(0,nrow(left.reg),14)
        mins<-colMins(left.reg[,-1])
        for (i in 1:14){
            scale.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
        }
     }
  
     if (sampling=="dendrogram"){
        clusters<-hclust(dist(scale.left.reg),method="complete")
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
    } else {
        set.seed(333)
        train.ind<-sample(1:nrow(left.reg),3000)
        train<-left.reg[train.ind,-1]
    }
  
    train<-as.matrix(train)
    test<-as.matrix(test)
     
    if (scale=="Z"){
        sds<-colSds(train)
        means<-colMeans(train)
        scale.test<-matrix(0,nrow(test),ncol(test))
        for(i in 1:ncol(test)){
            scale.test[,i]<-(test[,i]-means[i])/sds[i]
        }
        scale.train<-scale(train)

    } else {
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
    }
  return(list(scale.train,scale.test))
}

narr<-data.narr(narr1,narr0,sampling="random",scale="01")
scale.test<-narr[[2]]
scale.train<-narr[[1]]

############################informative
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 info1<-read.csv('informative1.csv')
 info0<-read.csv('informative0.csv')
 test.novel<-info0
 data.info<-function(data1,data0,sampling,scale){
     info1<-data1
     info0<-data0     
     index<-seq(1,nrow(info1),1)
     info1.w.ind<-cbind(index,info1)
     set.seed(123498)
     test.reg.ind<-sample(1: nrow(info1.w.ind),nrow(info0))
     test.reg<-info1.w.ind[test.reg.ind,]
     test<-rbind(test.reg[,-1],test.novel)
     left.reg<-as.matrix(info1.w.ind[-test.reg.ind,]) 
     
     if (scale=="Z"){
        scale.left.reg<-scale(left.reg[,-1])
     } else {
        scale.left.reg<-matrix(0,nrow(left.reg),14)
        mins<-colMins(left.reg[,-1])
        for (i in 1:14){
            scale.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
        }
     }
  
     if (sampling=="dendrogram"){
        clusters<-hclust(dist(scale.left.reg),method="complete")
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
    } else {
        set.seed(333)
        train.ind<-sample(1:nrow(left.reg),4000)
        train<-left.reg[train.ind,-1]
    }
  
    train<-as.matrix(train)
    test<-as.matrix(test)
     
    if (scale=="Z"){
        sds<-colSds(train)
        means<-colMeans(train)
        scale.test<-matrix(0,nrow(test),ncol(test))
        for(i in 1:ncol(test)){
            scale.test[,i]<-(test[,i]-means[i])/sds[i]
        }
        scale.train<-scale(train)

    } else {
        unit.train<-matrix(0,4000,ncol(test))
        mins2<-colMins(train)
        maxs2<-matrix(0,1,ncol(test))
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
    }
  return(list(scale.train,scale.test))
}

info<-data.info(info1,info0,sampling="dendrogram",scale="01")
scale.test<-info[[2]]
scale.train<-info[[1]]

#########################Argumentative
 setwd('C:/Users/luyao_peng/Desktop/kpca_file/OutlierDetection')
 argu1<-read.csv('argumentative1.csv')
 argu0<-read.csv('argumentative0.csv')
 test.novel<-argu0
 data.info<-function(data1,data0,sampling,scale){
     argu1<-data1
     argu0<-data0     
     index<-seq(1,nrow(argu1),1)
     argu1.w.ind<-cbind(index,argu1)
     set.seed(1234555)
     test.reg.ind<-sample(1: nrow(argu1.w.ind),nrow(argu0))
     test.reg<-argu1.w.ind[test.reg.ind,]
     test<-rbind(test.reg[,-1],test.novel)
     left.reg<-argu1.w.ind[-test.reg.ind,] 
     
     if (scale=="Z"){
        scale.left.reg<-scale(left.reg[,-1])
     } else {
        scale.left.reg<-matrix(0,nrow(left.reg),14)
        left.reg<-as.matrix(left.reg)
        mins<-colMins(left.reg[,-1])
        for (i in 1:14){
            scale.left.reg[,i]<-(left.reg[,i+1]-mins[i])/max(left.reg[,i+1]-mins[i])
        }
     }
  
     if (sampling=="dendrogram"){
        clusters<-hclust(dist(scale.left.reg),method="complete")
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
    } else {
        set.seed(333)
        train.ind<-sample(1:nrow(left.reg),2000)
        train<-left.reg[train.ind,-1]
    }
  
    train<-as.matrix(train)
    test<-as.matrix(test)
     
    if (scale=="Z"){
        sds<-colSds(train)
        means<-colMeans(train)
        scale.test<-matrix(0,nrow(test),ncol(test))
        for(i in 1:ncol(test)){
            scale.test[,i]<-(test[,i]-means[i])/sds[i]
        }
        scale.train<-scale(train)

    } else {
        unit.train<-matrix(0,2000,ncol(test))
        mins2<-colMins(train)
        maxs2<-matrix(0,1,ncol(test))
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
    }
  return(list(scale.train,scale.test))
}

argu<-data.info(argu1,argu0,sampling="dendrogram",scale="01")
scale.test<-argu[[2]]
scale.train<-argu[[1]]
