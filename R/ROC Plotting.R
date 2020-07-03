##############using training data as thr to plot ROC
#thr<-result$serr[order(-result$serr)]
#tpr<-numeric(3000)
#fpr<-numeric(3000)

##############KPCA using testing data as thr to plot ROC

result<-kpcabound(scale.train,scale.test,1.1,10,0) #narrative

n<-184
err<-result$err
ind<-matrix(c(rep(1,184),rep(0,184)),nrow(scale.test),1)
thr<-result$err[order(-result$err)]
tpr<-numeric(length(thr))
fpr<-numeric(length(thr))
acc<-numeric(length(thr))
for (i in 1:length(thr)){
 tp=0
 fp=0
 for (j in 1:nrow(scale.test)){
      if(err[j]> thr[i] && ind[j]==0){
           tp=tp+1
         } else if (err[j]> thr[i] && ind[j]==1) {
           fp=fp+1
         } else {
           tp=tp
           fp=fp
         }
    }
 tpr[i]<-tp/184
 fpr[i]<-fp/184
 acc[i]<-(tp+n-fp)/368
}
par(cex.axis=2, cex.lab=2, mar=(c(5, 8, 4, 2)+0.1), lab=c(7,3,7))
plot(fpr,tpr,type='l',col="blue",xlab="False Positive Rate",ylab="True Positive Rate")
 
# label<-ind
# pred <- prediction(result$err, label)
# perf<- performance(pred, 'tpr', 'fpr')
# plot(perf)

result<-kpcabound(scale.train,scale.test,0.1,15,0) #informative
#4000 information# result<-kpcabound(scale.train,scale.test,0.09,6,0)

n<-nrow(test.novel)
err<-result$err
ind<-matrix(c(rep(1,n),rep(0,n)),nrow(scale.test),1)
thr<-result$err[order(-result$err)]
tpr<-numeric(length(thr))
fpr<-numeric(length(thr))
acc<-numeric(length(thr))
for (i in 1:length(thr)){
 tp=0
 fp=0
 for (j in 1:nrow(scale.test)){
      if(err[j]> thr[i] && ind[j]==0){
           tp=tp+1
         } else if (err[j]> thr[i] && ind[j]==1) {
           fp=fp+1
         } else {
           tp=tp
           fp=fp
         }
    }
 tpr[i]<-tp/n
 fpr[i]<-fp/n
 acc[i]<-(tp+n-fp)/(2*n)
}
windows()
par(cex.axis=2, cex.lab=2, mar=(c(5, 8, 4, 2)+0.1), lab=c(7,3,7))
plot(fpr,tpr,type='l',col="blue",xlab="False Positive Rate",ylab="True Positive Rate")


result<-kpcabound(scale.train,scale.test,0.12,10,0) #argumentative
n<-nrow(test.novel) 
err<-result$err
ind<-matrix(c(rep(1,n),rep(0,n)),nrow(scale.test),1)
thr<-result$err[order(-result$err)]
tpr<-numeric(length(thr))
fpr<-numeric(length(thr))
acc<-numeric(length(thr))
for (i in 1:length(thr)){
 tp=0
 fp=0
 for (j in 1:nrow(scale.test)){
      if(err[j]> thr[i] && ind[j]==0){
           tp=tp+1
         } else if (err[j]> thr[i] && ind[j]==1) {
           fp=fp+1
         } else {
           tp=tp
           fp=fp
         }
    }
 tpr[i]<-tp/n
 fpr[i]<-fp/n
 acc[i]<-(tp+n-fp)/(2*n)
}
windows()
par(cex.axis=2, cex.lab=1.8, mar=(c(5, 8, 4, 2)+0.1), lab=c(7,3,7))
plot(fpr,tpr,type='l',col="blue",xlab="False Positive Rate",ylab="True Positive Rate")


#####################OneClass SVM ROC tuning NU
nu<-seq(0.001,1,0.005)
tpr1<-numeric(length(nu))
fpr1<-numeric(length(nu))
acc<-numeric(length(nu))
class<-rep(1,3000)
for (i in 1:length(nu)){
  model <- svm(scale.train, class,type='one-classification',gamma=0.0001,nu=nu[i])
  pred <- predict(model, scale.test)
  class1<-matrix(1,184,1)
  class0<-matrix(0,184,1)
  true.class<-rbind(class1,class0)
  predict<-cbind(pred,true.class)
  colnames(predict)<-c("pred","true.class")
  correct.predict<-ifelse(pred==true.class,1,0)
  acc[i]<-sum(correct.predict)/(2*184)
  fpr1[i]<-1-sum(correct.predict[1:184])/184
  tpr1[i]<-sum(correct.predict[185:368])/184
}

par(new = TRUE)
plot(fpr1,tpr1,axes=FALSE,type='l',lty=1,col='red',xlab='False Positive Rate',ylab='True Positive Rate',xlim=c(0,1),ylim=c(0,1))
legend(0.4,0.2,c("KPCA","1-class SVM"),lty=c(1,1),col=c("blue","red"))
title(main="ROC of KPCA and 1-class SVM for Narrative",cex=2)


