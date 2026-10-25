# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
library(survival)
set.seed(2026101302)
out<-'primary';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-readRDS('descriptive/broad_derived_with_ids.rds')
d$hist_num<-as.integer(substr(d[['Histologic Type ICD-O-3']],1,4))
d$surg_num<-as.integer(d[['RX Summ--Surg Prim Site (1998-2022)']])
stopifnot(all(d$typed3== (d$hist_num%in%c(8144,8163,8255)&d$current)))
sc<-data.frame(name=c('supported3_all','supported3_surgery','modern1_all','modern1_surgery','narrow3_all','narrow3_surgery','strict3_surgery','old3_surgery','recent_supported3_surgery'),
 start=c(rep(2010,8),2018),end=c(2019,2019,2022,2022,2019,2019,2019,2017,2019),
 surg=c(FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE),
 narrow=c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE),
 strict=c(rep(FALSE,6),TRUE,FALSE,FALSE),
 horizon=c(36,36,12,12,36,36,36,36,36),comparator='NOS',grade='none',complete=FALSE,ipcw=FALSE)
write.csv(sc,file.path(out,'analysis_specification.csv'),row.names=FALSE)
ss<-readLines(file.path(Sys.getenv('SEER_CODE_ROOT'),'R','04_grade_comparator_ipcw.R'))
i<-grep('^aj<-',ss);j<-grep('^res<-diag<-draws',ss)
eval(parse(text=ss[i:(j-1)]))
# Vectorized, algebraically identical weighted AJ; checked against survival below.
aj<-function(time,status,w,tau){
 oo<-order(time);tt<-time[oo];st<-status[oo];ww<-w[oo]
 first<-!duplicated(tt);times<-tt[first]
 den<-rev(cumsum(rev(ww)))[first]
 da<-as.vector(rowsum(ww*(st!=0),tt,reorder=FALSE));dc<-as.vector(rowsum(ww*(st==1),tt,reorder=FALSE))
 use<-times<=tau & den>0;ha<-da[use]/den[use];hc<-dc[use]/den[use]
 if(!length(ha))return(0)
 sum(c(1,head(cumprod(1-ha),-1))*hc)
}
original_prepare<-prepare
prepare<-function(a){
 dd<-d
 if(a$narrow)dd<-dd[dd$hist_num!=8255,]
 if(a$strict)dd<-dd[dd$surg_num%in%c(20:49,60),]
 fun<-original_prepare;environment(fun)<-list2env(list(d=dd),parent=environment(original_prepare))
 fun(a)
}
res<-diag<-draws<-verification<-list();B<-bootstrap_n(1000L)
for(k in seq_len(nrow(sc))){
 a<-sc[k,];pr<-prepare(a);q<-pr$q;pt<-calculate(pr,a)
 for(g in 0:1){
  ii<-q$typed==g;ww<-pt$w[ii];wa<-pt$w[ii&q$time>=a$horizon]
  diag[[length(diag)+1]]<-data.frame(scenario=a$name,typed=g,n=sum(ii),cancer=sum(q$status[ii]==1),other=sum(q$status[ii]==2),at_risk=sum(ii&q$time>=a$horizon),ESS=sum(ww)^2/sum(ww^2),at_risk_ESS=if(length(wa))sum(wa)^2/sum(wa^2)else 0,max_share=max(ww)/sum(ww),point_warnings=pt$warning_text)
  z<-q[ii,];z$w<-ww;z$ev<-factor(z$status,0:2,c('censor','cancer','other'));z$record<-seq_len(nrow(z))
  sf<-survfit(Surv(time,ev)~1,data=z,weights=w,id=record);jj<-max(which(sf$time<=a$horizon))
  diff<-abs(sf$pstate[jj,match('cancer',sf$states)]-pt$vals[g+1]);stopifnot(diff<1e-10)
  verification[[length(verification)+1]]<-data.frame(scenario=a$name,typed=g,independent_AJ_error=diff)
 }
 chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);bs<-matrix(NA_real_,B,3);wn<-integer(B);err<-character(B)
 for(b in seq_len(B)){
  ix<-unlist(chunks[sample(ids,length(ids),TRUE)],use.names=FALSE)
  v<-tryCatch(calculate(pr,a,ix),error=function(e){err[b]<<-conditionMessage(e);NULL})
  if(!is.null(v)){bs[b,]<-v$vals;wn[b]<-v$warnings}
 }
 ok<-complete.cases(bs);ci<-apply(bs[ok,,drop=FALSE],2,quantile,c(.025,.975))
 res[[k]]<-data.frame(scenario=a$name,estimand=names(pt$vals),estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],valid=sum(ok),warnings=sum(wn[ok]>0),B=B)
 draws[[k]]<-data.frame(scenario=a$name,replicate=seq_len(B),control=bs[,1],detailed=bs[,2],RD=bs[,3],warnings=wn,error=err)
 write.csv(do.call(rbind,res),file.path(out,'risks.csv'),row.names=FALSE)
 write.csv(do.call(rbind,diag),file.path(out,'diagnostics.csv'),row.names=FALSE)
 write.csv(do.call(rbind,draws),file.path(out,'bootstrap.csv'),row.names=FALSE)
 cat('DONE',a$name,'N',nrow(q),'RD',pt$vals[3],'CI',ci[,3],'valid',sum(ok),'warn',sum(wn[ok]>0),'\n');flush.console()
}
write.csv(do.call(rbind,verification),file.path(out,'point_validation.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
cat('COMPLETE bounded analysis set\n')
