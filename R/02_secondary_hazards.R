# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
library(survival)
out<-'descriptive'
d<-readRDS(file.path(out,'broad_derived_with_ids.rds'))
rows<-list(); balance<-list(); risks<-list()
for(start in c(2010,2018)) {
 q<-d[d$eligible & d$year>=start & d$year<=2022 & (d$typed3|d$nos),]
 q$typed<-as.integer(q$typed3); q<-droplevels(q)
 ps<-glm(typed~splines::ns(year,3)+splines::ns(age,3)+sex+race+stage,data=q,family=binomial())
 p<-predict(ps,type='response'); stopifnot(length(p)==nrow(q),all(is.finite(p)))
 q$weight<-ifelse(q$typed==1,1-p,p)
 x<-model.matrix(~age+year+sex+race+stage,q)[,-1,drop=FALSE]
 for(j in seq_len(ncol(x))) {
  a<-x[q$typed==1,j]; b<-x[q$typed==0,j]; sd0<-sqrt((var(a)+var(b))/2)
  balance[[length(balance)+1]]<-data.frame(start=start,variable=colnames(x)[j],SMD_before=(mean(a)-mean(b))/sd0,SMD_after=(weighted.mean(a,q$weight[q$typed==1])-weighted.mean(b,q$weight[q$typed==0]))/sd0)
 }
 for(interval in list(c(0,Inf),c(0,12),c(12,36),c(36,60))) {
  a<-interval[1]; b<-interval[2]; z<-q[q$time>a,]; z$tt<-pmin(z$time,b)-a; z$ev<-as.integer(z$status==1&z$time<=b)
  fit<-coxph(Surv(tt,ev)~typed+cluster(id),data=z,weights=weight)
  sm<-summary(fit); r<-sm$conf.int['typed',]
  rows[[length(rows)+1]]<-data.frame(start=start,interval_start=a,interval_end=b,n=nrow(z),cancer_events=sum(z$ev),HR=r[1],lower95=r[3],upper95=r[4],effective_n_typed=sum(z$weight[z$typed==1])^2/sum(z$weight[z$typed==1]^2),effective_n_NOS=sum(z$weight[z$typed==0])^2/sum(z$weight[z$typed==0]^2))
 }
 # Weighted Aalen-Johansen probabilities; descriptive point estimates only.
 for(g in 0:1) {
  z<-q[q$typed==g,]; z$event<-factor(z$status,levels=0:2,labels=c('censor','cancer','other'))
  fit<-survfit(Surv(time,event)~1,data=z,weights=weight,id=id)
  for(t in c(36,60)) {
   j<-max(which(fit$time<=t)); idx<-match('cancer',fit$states)
   risks[[length(risks)+1]]<-data.frame(start=start,typed=g,months=t,weighted_cancer_CIF=fit$pstate[j,idx],unweighted_at_risk=sum(z$time>=t),support_warning=sum(z$time>=t)<20)
  }
 }
}
write.csv(do.call(rbind,balance),file.path(out,'overlap_balance.csv'),row.names=FALSE)
write.csv(do.call(rbind,rows),file.path(out,'overlap_weighted_piecewise_cox.csv'),row.names=FALSE)
write.csv(do.call(rbind,risks),file.path(out,'overlap_weighted_cif.csv'),row.names=FALSE)
print(do.call(rbind,rows)); print(do.call(rbind,risks))
