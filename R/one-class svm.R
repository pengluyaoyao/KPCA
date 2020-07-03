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

#class<-rep(1,3000)
#model <- svm(scale.train, class,type='one-classification',gamma=0.5,cross=10)
#pred <- predict(model, scale.train)

set.seed(3333)
cv<-sample(1:3000,3000)
cv.scale.train<-cbind(cv,scale.train)

gamma<-seq(0.0001,0.2,0.001)

cv.gamma<-function(cv.scale.train,gamma,nu,k){
for (g in 1:length(gamma)){
correct<-0
false<-0
n<-nrow(cv.scale.train)
m<-n/k
n1<-seq(1,n,m)
n2<-seq(m,n,m)
for(i in 1:length(n1)){    
    train<-cv.scale.train[which(cv>n2[i] | cv<n1[i] ),-1]
    test<-cv.scale.train[which(cv>=n1[i] & cv <= n2[i]),-1]
    class<-rep(1,nrow(train))
    model<-svm(train,class,scale=FALSE,type='one-classification',
       kernel='radial',gamma=gamma[g],nu=nu)
    predict<-predict(model,test)
    predict<-ifelse(predict=='TRUE',1,0)
    correct.new<-sum(predict==1)
    false.new<-nrow(test)-correct.new
    correct<-correct+correct.new
    false<-false+false.new
   }
  correct<-correct/(k*m)
  false<-false/(k*m) 
  print(c(gamma[g],false,correct)) 
  }
}


result<-cv.gamma(cv.scale.train,gamma,0.1,10)




###model on test dataset 
class<-rep(1,3000)
model <- svm(scale.train, class,type='one-classification',gamma=0.0001,nu=0.1)
pred <- predict(model, scale.test,probability=TRUE)
class1<-matrix(1,103,1)
class0<-matrix(0,103,1)
class<-rbind(class1,class0)
predict<-cbind(pred,class)
colnames(predict)<-c("pred","class")
acc<-ifelse(pred==class,1,0)
total.acc<-sum(acc)/206
acc.c1<-sum(acc[1:103])/103
acc.c2<-1-sum(acc[104:206])/103




########informative 3000
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

set.seed(112)
cv<-sample(1:3000,3000)
cv.scale.train<-cbind(cv,scale.train)
gamma<-seq(0.0001,1.0,0.01)
result<-cv.gamma(cv.scale.train,gamma,10)

class<-rep(1,3000)
model <- svm(scale.train, class,type='one-classification',gamma=0.0001,nu=0.1)
pred <- predict(model, scale.test)
class1<-matrix(1,468,1)
class0<-matrix(0,468,1)
class<-rbind(class1,class0)
predict<-cbind(pred,class)
colnames(predict)<-c("pred","class")
acc<-ifelse(pred==class,1,0)
acc<-sum(acc)/936
acc.c1<-sum(predict[c(1:468),1])/468
acc.c2<-1-sum(predict[c(469:936),1])/468

#####informative4000
set.seed(112)
cv<-sample(1:4000,4000)
cv.scale.train<-cbind(cv,scale.train)
gamma<-seq(0.2,0.5,0.01)
result<-cv.gamma(cv.scale.train,gamma,0.1,10)

class<-rep(1,4000)
model <- svm(scale.train, class,type='one-classification',gamma=0.5,nu=0.1)
pred <- predict(model, scale.test)
class1<-matrix(1,468,1)
class0<-matrix(0,468,1)
class<-rbind(class1,class0)
predict<-cbind(pred,class)
colnames(predict)<-c("pred","class")
acc<-ifelse(pred==class,1,0)
acc<-sum(acc)/936
acc.c1<-sum(predict[c(1:468),1])/468
acc.c2<-1-sum(predict[c(469:936),1])/468


####argumentative
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

set.seed(111)
cv<-sample(1:2000,2000)
cv.scale.train<-cbind(cv,scale.train)

gamma<-seq(0.0001,1.0,0.01)
result<-cv.gamma(cv.scale.train,gamma,10)

class<-rep(1,2000)
model <- svm(scale.train, class,type='one-classification',gamma=0.00005,nu=0.1)
pred <- predict(model, scale.test)
class1<-matrix(1,103,1)
class0<-matrix(0,103,1)
class<-rbind(class1,class0)
predict<-cbind(pred,class)
colnames(predict)<-c("pred","class")
acc<-ifelse(pred==class,1,0)
acc<-sum(acc)/206
acc.c1<-sum(predict[c(1:103),1])/103
acc.c2<-1-sum(predict[c(104:206),1])/103
