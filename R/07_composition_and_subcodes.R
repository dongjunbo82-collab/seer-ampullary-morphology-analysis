# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
library(survival)
set.seed(2026101403)
out <- 'composition'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- readRDS('descriptive/broad_derived_with_ids.rds')
d$hist_num <- as.integer(substr(d[['Histologic Type ICD-O-3']],1,4))
aj <- function(t,st,w,tau=36){
  oo<-order(t);t<-t[oo];st<-st[oo];w<-w[oo];first<-!duplicated(t)
  den<-rev(cumsum(rev(w)))[first];times<-t[first]
  da<-as.vector(rowsum(w*(st!=0),t,reorder=FALSE));dc<-as.vector(rowsum(w*(st==1),t,reorder=FALSE))
  use<-times<=tau;ha<-da[use]/den[use];hc<-dc[use]/den[use]
  sum(c(1,head(cumprod(1-ha),-1))*hc)
}
q0 <- d[d$eligible & d$year>=2010 & d$year<=2019 & (d$typed3|d$nos),]
print(with(q0,table(typed3,surgery,status)))
write.csv(as.data.frame(with(q0,table(code=hist_num,surgery,status))),file.path(out,'support_cells.csv'),row.names=FALSE)
spec<-data.frame(name=c('primary','exclude8255'),start=2010,end=2019,horizon=36,
  target='Pooled empirical covariates, with common pooled marginal surgery mixture',
  comparator='NOS 8140',B=1000,interpretation='Hypothetical noncausal composition standardization')
write.csv(spec,file.path(out,'analysis_specification.csv'),row.names=FALSE)
prepare <- function(narrow=FALSE){
 q<-q0[q0$surgery%in%c('No','Yes') & (!narrow | q0$hist_num!=8255),]
 q$G<-as.integer(q$typed3);q$S<-as.integer(q$surgery=='Yes');q<-droplevels(q)
 fm<-~G*S + splines::ns(age,3)+splines::ns(year,3)+sex+race+stage
 X<-model.matrix(fm,q)[,-1,drop=FALSE]
 tt<-terms(fm); bs<-list()
 for(g in 0:1)for(s in 0:1){z<-q;z$G<-g;z$S<-s;bs[[paste0(g,s)]]<-model.matrix(tt,z)[,-1,drop=FALSE]}
 # Prediction bases must use the same knots as the fitted matrices.
 stopifnot(max(abs(bs[['00']][,grep('ns',colnames(X)),drop=FALSE]-X[,grep('ns',colnames(X)),drop=FALSE]))<1e-12)
 list(q=q,X=X,bs=bs)
}
calculate <- function(pr,ix=seq_len(nrow(pr$q))){
 q<-pr$q[ix,];X<-pr$X[ix,,drop=FALSE];warn<-character()
 ans<-withCallingHandlers({
   fits<-lapply(1:2,function(cause)coxph(Surv(pmin(q$time,36),as.integer(q$status==cause&q$time<=36))~ridge(X,theta=0.1),
      ties='breslow',singular.ok=TRUE,control=coxph.control(iter.max=100)))
   co<-lapply(fits,function(f){b<-coef(f);b[is.na(b)]<-0;if(any(!is.finite(b)))stop('nonfinite coefficient');b})
   bh<-lapply(fits,function(f)basehaz(f,centered=FALSE))
   times<-sort(unique(c(bh[[1]]$time,bh[[2]]$time)));times<-times[times<=36]
   inc<-lapply(bh,function(h){j<-findInterval(times,h$time);v<-c(0,h$hazard)[j+1];diff(c(0,v))})
   risks<-sapply(pr$bs,function(mat){
      xx<-mat[ix,,drop=FALSE];r1<-exp(drop(xx%*%co[[1]]));r2<-exp(drop(xx%*%co[[2]]))
      if(any(!is.finite(c(r1,r2))))stop('nonfinite predicted hazard')
      ss<-rep(1,nrow(q));ff<-rep(0,nrow(q))
      for(j in seq_along(times)){
        h1<-r1*inc[[1]][j];h2<-r2*inc[[2]][j];ht<-h1+h2
        share<-ifelse(ht>0,h1/ht,0);ff<-ff+ss*(-expm1(-ht))*share;ss<-ss*exp(-ht)
      }
      mean(ff)
   })
   p0<-mean(q$S[q$G==0]);p1<-mean(q$S[q$G==1]);p<-mean(q$S)
   nos_mix<-(1-p0)*risks['00']+p0*risks['01'];det_mix<-(1-p1)*risks['10']+p1*risks['11']
   nos_common<-(1-p)*risks['00']+p*risks['01'];det_common<-(1-p)*risks['10']+p*risks['11']
   rd0<-unname(det_mix-nos_mix);rd1<-unname(det_common-nos_common)
   # Relative attenuation is undefined when the baseline contrast is negligible or reverses.
   ratio<-if(abs(rd0)>=.01 && rd0*rd1>=0)1-rd1/rd0 else NA_real_
   list(vals=c(NOS_group_mixture=unname(nos_mix),detailed_group_mixture=unname(det_mix),RD_group_mixture=rd0,
     NOS_common_mixture=unname(nos_common),detailed_common_mixture=unname(det_common),RD_common_mixture=rd1,
     contrast_change=rd1-rd0,relative_attenuation=ratio),stratum_risks=risks,p=c(NOS=p0,detailed=p1,pooled=p))
 },warning=function(e){warn<<-c(warn,conditionMessage(e));invokeRestart('muffleWarning')})
 ans$warnings<-paste(unique(warn),collapse='; ');ans
}
results<-draws<-diag<-list()
cl<-NULL
if(Sys.getenv('SEER_MODE')=='full'){
 cl<-parallel::makePSOCKcluster(4)
 parallel::clusterEvalQ(cl,library(survival))
 parallel::clusterExport(cl,'calculate')
}
for(narrow in c(FALSE,TRUE)){
 nm<-if(narrow)'exclude8255'else'primary';pr<-prepare(narrow);q<-pr$q;pt<-calculate(pr)
 chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);B<-bootstrap_n(1000L)
 bs<-matrix(NA_real_,B,length(pt$vals));colnames(bs)<-names(pt$vals);warnings<-errors<-character(B)
 samples<-lapply(seq_len(B),function(b)unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE))
 worker<-function(ix){tryCatch({v<-calculate(pr,ix);list(vals=v$vals,warnings=v$warnings,error='')},error=function(e)list(vals=NULL,warnings='',error=conditionMessage(e)))}
 if(!is.null(cl))parallel::clusterExport(cl,'pr')
 for(batch in split(seq_len(B),ceiling(seq_len(B)/100))){
  vs<-if(is.null(cl))lapply(samples[batch],worker)else parallel::parLapplyLB(cl,samples[batch],worker)
  for(k in seq_along(batch)){b<-batch[k];v<-vs[[k]];if(!is.null(v$vals))bs[b,]<-v$vals;warnings[b]<-v$warnings;errors[b]<-v$error}
  cat(nm,max(batch),'of',B,'\n');flush.console()
 }
 ci<-apply(bs,2,function(z)if(sum(is.finite(z))>=950)quantile(z,c(.025,.975),na.rm=TRUE)else c(NA,NA))
 results[[nm]]<-data.frame(scenario=nm,estimand=names(pt$vals),estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],
  finite_replicates=colSums(is.finite(bs)),B=B,warning_replicates=sum(nzchar(warnings)),failed=sum(nzchar(errors)))
 draws[[nm]]<-data.frame(scenario=nm,replicate=1:B,bs,warnings=warnings,error=errors)
 diag[[nm]]<-data.frame(scenario=nm,group=c('NOS','Detailed'),n=c(sum(q$G==0),sum(q$G==1)),
  surgery_pct=100*pt$p[1:2],pooled_surgery_pct=100*pt$p[3],at36=c(sum(q$G==0&q$time>=36),sum(q$G==1&q$time>=36)),
  no_surgery_n=c(sum(q$G==0&q$S==0),sum(q$G==1&q$S==0)),no_surgery_at36=c(sum(q$G==0&q$S==0&q$time>=36),sum(q$G==1&q$S==0&q$time>=36)),
  modeled_no_surgery_risk=pt$stratum_risks[c('00','10')],modeled_surgery_risk=pt$stratum_risks[c('01','11')],point_warnings=pt$warnings)
 write.csv(do.call(rbind,results),file.path(out,'gcomp_results.csv'),row.names=FALSE)
 write.csv(do.call(rbind,draws),file.path(out,'gcomp_bootstrap.csv'),row.names=FALSE)
 write.csv(do.call(rbind,diag),file.path(out,'gcomp_diagnostics.csv'),row.names=FALSE)
 print(results[[nm]])
}
if(!is.null(cl))parallel::stopCluster(cl)
# Descriptive subtype comparison uses the same supported diagnosis window.
# No covariate-adjusted subtype model or inferential claim is introduced.
subres<-list();subdraw<-list()
for(surg in c(FALSE,TRUE)){
 q<-q0[q0$hist_num%in%c(8144,8163) & (!surg|q0$surgery=='Yes'),]
 chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks)
 calc<-function(ix){z<-q[ix,];r<-sapply(c(8144,8163),function(code){a<-z$hist_num==code;if(!any(a))return(NA_real_);aj(z$time[a],z$status[a],rep(1,sum(a)))});c(intestinal=r[1],pancreatobiliary=r[2],RD=r[2]-r[1])}
 pt<-calc(seq_len(nrow(q)));bs<-t(replicate(bootstrap_n(1000L),calc(unlist(chunks[sample(ids,length(ids),TRUE)],use.names=FALSE))))
 ci<-apply(bs,2,quantile,c(.025,.975),na.rm=TRUE)
 for(k in 1:2){code<-c(8144,8163)[k];z<-q[q$hist_num==code,]
  subres[[length(subres)+1]]<-data.frame(stratum=if(surg)'Recorded surgery'else'All',code=code,n=nrow(z),cancer_events_by36=sum(z$status==1&z$time<=36),other_events_by36=sum(z$status==2&z$time<=36),at36=sum(z$time>=36),surgery_n=sum(z$surgery=='Yes'),risk=pt[k],lower95=ci[1,k],upper95=ci[2,k])
 }
 subdraw[[length(subdraw)+1]]<-data.frame(stratum=if(surg)'Recorded surgery'else'All',estimand=names(pt),estimate=pt,lower95=ci[1,],upper95=ci[2,])
}
write.csv(do.call(rbind,subres),file.path(out,'subcode_support_risks.csv'),row.names=FALSE)
write.csv(do.call(rbind,subdraw),file.path(out,'subcode_contrasts.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('COMPLETE\n')
