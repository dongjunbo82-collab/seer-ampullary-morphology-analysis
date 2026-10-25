# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
library(survival)
set.seed(20261013)
out<-'sensitivity';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-readRDS('descriptive/broad_derived_with_ids.rds')
# Separate era-specific grade sources, no cross-era grade concatenation.
old<-d[['Grade Recode (thru 2017)']]
d$grade_old<-factor(ifelse(grepl('^Well',old),'I',ifelse(grepl('^Moderately',old),'II',ifelse(grepl('^Poorly|^Undifferentiated',old),'III/IV','Unknown'))))
gr<-function(x)factor(ifelse(x%in%c('1','2','3'),x,'Unknown'),levels=c('1','2','3','Unknown'))
d$grade_clin<-gr(d[['Grade Clinical (2018+)']]);d$grade_path<-gr(d[['Grade Pathological (2018+)']])
q0<-d[d$eligible&d$year>=2010&d$year<=2022&(d$typed3|d$nos),]
ga<-list()
for(nm in c('grade_old','grade_clin','grade_path')) {
 q<-q0[if(nm=='grade_old')q0$year<=2017 else q0$year>=2018,]
 x<-as.data.frame(table(group=ifelse(q$typed3,'Detailed','NOS'),surgery=q$surgery,grade=q[[nm]]));x$source<-nm;ga[[nm]]<-x
}
write.csv(do.call(rbind,ga),file.path(out,'grade_by_group_surgery.csv'),row.names=FALSE)
yr<-do.call(rbind,lapply(split(q0,list(q0$year,q0$typed3),drop=TRUE),function(z)data.frame(year=z$year[1],typed=as.integer(z$typed3[1]),n=nrow(z),max_followup=max(z$time),at36=sum(z$time>=36),at60=sum(z$time>=60),censored_before36=sum(z$status==0&z$time<36),censored_before60=sum(z$status==0&z$time<60))))
write.csv(yr,file.path(out,'horizon_support_by_year.csv'),row.names=FALSE)
sc<-data.frame(name=c('all_reference','all_reference_surgery','old_surgery_base','old_surgery_grade','recent_clinical_grade','recent_surgery_pathgrade','recent_surgery_path_complete','mature_base','mature_ipcw','mature_surgery_base','mature_surgery_ipcw','supported3_base','supported3_ipcw'),
 start=c(2010,2010,2010,2010,2018,2018,2018,2010,2010,2010,2010,2010,2010),
 end=c(2022,2022,2017,2017,2022,2022,2022,2017,2017,2017,2017,2019,2019),
 surg=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE),
 comparator=c('other_reference','other_reference',rep('NOS',11)),
 grade=c('none','none','none','grade_old','grade_clin','grade_path','grade_path',rep('none',6)),
 complete=c(rep(FALSE,6),TRUE,rep(FALSE,6)),horizon=c(60,60,60,60,36,36,36,60,60,60,60,36,36),
 ipcw=c(rep(FALSE,8),TRUE,FALSE,TRUE,FALSE,TRUE))
write.csv(sc,file.path(out,'analysis_specification.csv'),row.names=FALSE)
aj<-function(time,status,w,tau){s<-1;f<-0;for(t in sort(unique(time[status!=0&time<=tau]))){at<-time>=t;den<-sum(w[at]);da<-sum(w[time==t&status!=0]);dc<-sum(w[time==t&status==1]);f<-f+s*dc/den;s<-s*(1-da/den)};f}
prepare<-function(a){
 q<-d[d$eligible&d$year>=a$start&d$year<=a$end&if(a$comparator=='NOS')(d$typed3|d$nos)else d$current,]
 if(a$surg)q<-q[q$surgery=='Yes',]
 if(a$complete)q<-q[as.character(q[[a$grade]])!='Unknown',]
 q$typed<-as.integer(q$typed3);q<-droplevels(q)
 terms<-c('splines::ns(year,3)','splines::ns(age,3)',c('sex','race','stage')[sapply(q[c('sex','race','stage')],function(v)length(unique(v))>1)])
 if(a$grade!='none')terms<-c(terms,a$grade)
 X<-model.matrix(as.formula(paste('~',paste(terms,collapse='+'))),q)
 # Censoring model includes coding group, baseline variables, and recorded surgery.
 ct<-c('typed',terms,if(!a$surg)'surgery')
 C<-model.matrix(as.formula(paste('~',paste(ct,collapse='+'))),q)[,-1,drop=FALSE]
 list(q=q,X=X,C=C)
}
calculate<-function(pr,a,ix=seq_len(nrow(pr$q)),details=FALSE){
 q<-pr$q[ix,];X<-pr$X[ix,,drop=FALSE];C<-pr$C[ix,,drop=FALSE];warn<-character()
 ans<-withCallingHandlers({
 fit<-glm.fit(X,q$typed,family=binomial(),control=glm.control(maxit=100))
 if(!fit$converged||length(unique(q$typed))<2)stop('invalid coding fit')
 w<-ifelse(q$typed==1,1-fit$fitted.values,fit$fitted.values)
 cf<-NULL;gmin<-1;maxcw<-1
 if(a$ipcw){
   tt<-pmin(q$time,a$horizon);cc<-as.integer(q$status==0&q$time<=a$horizon)
   cf<-coxph(Surv(tt,cc)~C,ties='breslow',singular.ok=TRUE,control=coxph.control(iter.max=50))
   if(any(is.infinite(coef(cf))))stop('invalid censor fit')
   bh<-basehaz(cf,centered=TRUE);lp<-cf$linear.predictors
 }
 vals<-sapply(0:1,function(g){
   z<-q[q$typed==g,];wg<-w[q$typed==g]
   if(sum(wg)<1e-8)stop('zero overlap weight')
   if(!a$ipcw)return(aj(z$time,z$status,wg,a$horizon))
   ss<-1;ff<-0
   for(t in sort(unique(z$time[z$status!=0&z$time<=a$horizon]))){
     jj<-which(bh$time<t);hh<-if(length(jj))bh$hazard[max(jj)]else 0
     G<-exp(-hh*exp(lp[q$typed==g]));at<-z$time>=t
     gmin<<-min(gmin,G[at]);maxcw<<-max(maxcw,1/G[at])
     if(any(G[at]<1e-8))stop('censoring support collapse')
     ww<-wg/G;den<-sum(ww[at]);dc<-sum(ww[z$time==t&z$status==1]);da<-sum(ww[z$time==t&z$status!=0]);ff<-ff+ss*dc/den;ss<-ss*(1-da/den)
   };ff
 })
 list(vals=c(control=vals[1],detailed=vals[2],RD=vals[2]-vals[1]),w=w,minG=gmin,maxCW=maxcw)
 },warning=function(e){warn<<-c(warn,conditionMessage(e));invokeRestart('muffleWarning')})
 ans$warnings<-length(warn);ans$warning_text<-paste(unique(warn),collapse='; ');ans
}
res<-diag<-draws<-list();B<-bootstrap_n(1000L)
for(k in seq_len(nrow(sc))){
 a<-sc[k,];pr<-prepare(a);pt<-calculate(pr,a);q<-pr$q
 for(g in 0:1){ii<-q$typed==g;ww<-pt$w[ii];at<-ii&q$time>=a$horizon;wa<-pt$w[at]
 diag[[length(diag)+1]]<-data.frame(scenario=a$name,typed=g,n=sum(ii),cancer=sum(q$status[ii]==1),other=sum(q$status[ii]==2),at_risk=sum(at),ESS=sum(ww)^2/sum(ww^2),at_risk_ESS=if(length(wa))sum(wa)^2/sum(wa^2)else 0,min_weight=min(ww),max_weight=max(ww),max_normalized_share=max(ww)/sum(ww),minG=pt$minG,maxCW=pt$maxCW,point_warnings=pt$warning_text)
 }
 # Constant time-varying censor weights reduce exactly to static AJ when G is one.
 if(!a$ipcw)for(g in 0:1){z<-q[q$typed==g,];z$w<-pt$w[q$typed==g];z$ev<-factor(z$status,0:2,c('censor','cancer','other'));z$record<-seq_len(nrow(z));sf<-survfit(Surv(time,ev)~1,data=z,weights=w,id=record);j<-max(which(sf$time<=a$horizon));stopifnot(abs(sf$pstate[j,match('cancer',sf$states)]-pt$vals[g+1])<1e-10)}
 chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);bs<-matrix(NA_real_,B,3);wn<-integer(B);err<-character(B)
 for(b in seq_len(B)){ix<-unlist(chunks[sample(ids,length(ids),TRUE)],use.names=FALSE);v<-tryCatch(calculate(pr,a,ix),error=function(e){err[b]<<-conditionMessage(e);NULL});if(!is.null(v)){bs[b,]<-v$vals;wn[b]<-v$warnings}}
 ok<-complete.cases(bs);ci<-apply(bs[ok,,drop=FALSE],2,quantile,c(.025,.975))
 res[[k]]<-data.frame(scenario=a$name,estimand=names(pt$vals),estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],valid=sum(ok),warnings=sum(wn>0),B=B)
 draws[[k]]<-data.frame(scenario=a$name,replicate=seq_len(B),control=bs[,1],detailed=bs[,2],RD=bs[,3],warnings=wn,error=err)
 write.csv(do.call(rbind,res),file.path(out,'risks.csv'),row.names=FALSE);write.csv(do.call(rbind,diag),file.path(out,'diagnostics.csv'),row.names=FALSE);write.csv(do.call(rbind,draws),file.path(out,'bootstrap.csv'),row.names=FALSE)
 cat('DONE',a$name,'n',nrow(q),'RD',pt$vals[3],'valid',sum(ok),'\n');flush.console()
}
# Influence check for the recent surgery estimate, refitting weights for each omitted patient.
a<-sc[sc$name=='recent_surgery_pathgrade',];a$grade<-'none';pr<-prepare(a);pt<-calculate(pr,a);q<-pr$q
lv<-sapply(unique(q$id),function(id)calculate(pr,a,which(q$id!=id))$vals[3])
write.csv(data.frame(metric=c('full_RD','min_leave_one_patient_RD','max_leave_one_patient_RD','max_absolute_change'),value=c(pt$vals[3],min(lv),max(lv),max(abs(lv-pt$vals[3])))),file.path(out,'recent_surgery_influence.csv'),row.names=FALSE)
write.csv(data.frame(group=0:1,n=sapply(0:1,function(g)sum(q$typed==g)),ESS=sapply(0:1,function(g){w<-pt$w[q$typed==g];sum(w)^2/sum(w^2)}),max_share=sapply(0:1,function(g){w<-pt$w[q$typed==g];max(w)/sum(w)}),at36=sapply(0:1,function(g)sum(q$typed==g&q$time>=36))),file.path(out,'recent_surgery_weight_diagnostics.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('COMPLETE\n')
