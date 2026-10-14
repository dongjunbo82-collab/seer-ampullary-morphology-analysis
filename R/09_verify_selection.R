# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
out<-'coding_selection'
x<-readRDS(file.path(out,'selection_fit_private.rds'));q<-x$q
fit<-glm(x$formula,data=q,family=binomial(),control=glm.control(maxit=100))
stopifnot(max(abs(coef(fit)-x$fit$coefficients))<1e-10)
z<-x$results[x$results$estimand=='probability',]
check<-numeric(nrow(z))
for(i in seq_len(nrow(z))){
 nd<-q;v<-z$variable[i];nd[[v]]<-factor(rep(z$level[i],nrow(q)),levels=levels(q[[v]]))
 check[i]<-mean(predict(fit,newdata=nd,type='response'))
}
stopifnot(max(abs(check-z$estimate))<1e-12)
stopifnot(nrow(q)==1812,sum(q$typed)==289,all(q$stage!='Unknown'),all(q$surgery!='Unknown'))
write.csv(data.frame(check=c('formula versus matrix coefficients','predict versus matrix standardized probabilities','model records','detailed records'),
 value=c(max(abs(coef(fit)-x$fit$coefficients)),max(abs(check-z$estimate)),nrow(q),sum(q$typed))),file.path(out,'independent_validation.csv'),row.names=FALSE)
cat('PASS independent formula/predict verification\n')
