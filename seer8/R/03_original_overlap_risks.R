# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
library(survival)
set.seed(20261012)
out<-'original_overlap'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-readRDS('descriptive/broad_derived_with_ids.rds')
B<-bootstrap_n(1000L)
aj<-function(time,status,weight,horizon) {
 keep<-time<=horizon; tt<-sort(unique(time[keep])); ss<-1; ff<-0
 # weighted risk sets account for censoring at each observed time
 ix<-match(time,tt); at<-sum(weight)
 sums<-function(mask) {v<-numeric(length(tt)); z<-tapply(weight[mask],ix[mask],sum); if(length(z))v[as.integer(names(z))]<-z; v}
 exits<-sums(keep); e1<-sums(keep&status==1); ea<-sums(keep&status!=0)
 for(i in seq_along(tt)) {if(at>1e-12){ff<-ff+ss*e1[i]/at; ss<-ss*(1-ea[i]/at)};at<-at-exits[i]}
 ff
}
scenarios<-data.frame(name=c('modern_baseline','recent_baseline','modern_surgery','modern_no_surgery','recent_surgery','recent_no_surgery','modern_surgery_balanced'),
 start=c(2010,2018,2010,2010,2018,2018,2010),surg=c('All','All','Yes','No','Yes','No','All'),horizon=c(60,36,60,60,36,36,60),extra=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))
results<-list(); balances<-list(); boots<-list(); composition<-list(); pointQA<-list()
for(k in seq_len(nrow(scenarios))) {
 sc<-scenarios[k,]; q<-d[d$eligible & d$year>=sc$start & d$year<=2022 & (d$typed3|d$nos),]
 if(sc$surg!='All') q<-q[q$surgery==sc$surg,]
 q$typed<-as.integer(q$typed3); q<-droplevels(q)
 terms<-c('splines::ns(year,3)','splines::ns(age,3)',c('sex','race','stage')[sapply(q[c('sex','race','stage')],function(z)length(unique(z))>1)])
 if(sc$extra) terms<-c(terms,'surgery')
 X<-model.matrix(as.formula(paste('~',paste(terms,collapse='+'))),q)
 calc<-function(idx) {
  z<-q[idx,]; xx<-X[idx,,drop=FALSE]; warns<-character()
  fit<-withCallingHandlers(glm.fit(xx,z$typed,family=binomial(),control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
  if(!fit$converged||length(unique(z$typed))!=2)stop('nonconvergence or missing exposure')
  p<-fit$fitted.values; wt<-ifelse(z$typed==1,1-p,p)
  vals<-sapply(0:1,function(g){a<-z$typed==g;aj(z$time[a],z$status[a],wt[a],sc$horizon)})
  if(any(!is.finite(vals))||any(tapply(wt,z$typed,sum)<1e-8))stop('invalid weights')
  list(vals=c(NOS=vals[1],typed=vals[2],RD=vals[2]-vals[1]),weights=wt,warnings=length(warns))
 }
 pt<-calc(seq_len(nrow(q))); wt<-pt$weights
 # independent weighted multi-state calculation for every point estimate
 for(g in 0:1) {
  z<-q[q$typed==g,]; z$wt<-wt[q$typed==g]; z$ev<-factor(z$status,levels=0:2,labels=c('censor','cancer','other'));z$record<-seq_len(nrow(z))
  sf<-survfit(Surv(time,ev)~1,data=z,weights=wt,id=record)
  j<-max(which(sf$time<=sc$horizon)); chk<-sf$pstate[j,match('cancer',sf$states)]
  stopifnot(abs(chk-pt$vals[g+1])<1e-10)
  pointQA[[length(pointQA)+1]]<-data.frame(scenario=sc$name,typed=g,AJ_difference=abs(chk-pt$vals[g+1]))
  composition[[length(composition)+1]]<-data.frame(scenario=sc$name,typed=g,n=nrow(z),cancer_events=sum(z$status==1),other_events=sum(z$status==2),at_risk=sum(z$time>=sc$horizon),ESS=sum(z$wt)^2/sum(z$wt^2),unweighted_risk=aj(z$time,z$status,rep(1,nrow(z)),sc$horizon),age_mean=mean(z$age),surgery_pct=100*mean(z$surgery=='Yes'))
 }
 for(j in seq_len(ncol(X))[-1]) {
  a<-q$typed==1; sd0<-sqrt((var(X[a,j])+var(X[!a,j]))/2)
  balances[[length(balances)+1]]<-data.frame(scenario=sc$name,term=colnames(X)[j],SMD_before=(mean(X[a,j])-mean(X[!a,j]))/sd0,SMD_after=(weighted.mean(X[a,j],wt[a])-weighted.mean(X[!a,j],wt[!a]))/sd0)
 }
 chunks<-split(seq_len(nrow(q)),q$id); ids<-names(chunks)
 bs<-matrix(NA_real_,B,3); warns<-integer(B); errors<-character(B)
 for(b in seq_len(B)) {
  ix<-unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE)
  val<-tryCatch(calc(ix),error=function(e){errors[b]<<-conditionMessage(e);NULL})
  if(!is.null(val)){bs[b,]<-val$vals;warns[b]<-val$warnings}
 }
 ok<-complete.cases(bs);ci<-apply(bs[ok,,drop=FALSE],2,quantile,c(.025,.975))
 results[[k]]<-data.frame(scenario=sc$name,start=sc$start,horizon=sc$horizon,surgery_stratum=sc$surg,estimand=c('NOS','typed','RD'),estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],B=B,valid=sum(ok),warning_replicates=sum(warns>0),failed=sum(!ok))
 boots[[k]]<-data.frame(scenario=sc$name,replicate=seq_len(B),NOS=bs[,1],typed=bs[,2],RD=bs[,3],warnings=warns,error=errors)
 cat('FINISHED',sc$name,'n=',nrow(q),'typed=',sum(q$typed),'valid=',sum(ok),'RD=',pt$vals[3],'\n');flush.console()
}
wr<-function(x,n)write.csv(do.call(rbind,x),file.path(out,paste0(n,'.csv')),row.names=FALSE)
wr(results,'standardized_risks_bootstrap');wr(composition,'stratum_counts');wr(balances,'balance');wr(boots,'bootstrap_replicates');wr(pointQA,'independent_weighted_AJ_check')
saveRDS(scenarios,file.path(out,'analysis_scenarios.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
