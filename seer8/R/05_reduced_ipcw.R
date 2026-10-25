# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
library(survival)
lines<-readLines(file.path(Sys.getenv('SEER_CODE_ROOT'),'R','04_grade_comparator_ipcw.R'));end<-grep('^res<-diag',lines)[1]-1
eval(parse(text=lines[seq_len(end)]))
set.seed(20261014)
# Diagnostic-driven reduction, not selected on risk contrasts: sparse categorical
# censoring strata caused separation in the full model. Retain coding group,
# linear calendar year and age, and recorded-surgery indicator when applicable.
# Only 15, 12, and 19 censoring events occur before the respective horizons;
# unknown-surgery strata contain no censoring events. Preserve full-model fits.
sel<-sc[sc$ipcw,];res2<-diag2<-draws2<-list();B<-bootstrap_n(1000L)
for(k in seq_len(nrow(sel))){
 a<-sel[k,];pr<-prepare(a);q<-pr$q
 q$recorded_surgery<-as.integer(q$surgery=='Yes')
 terms<-c('typed','year','age',if(!a$surg)'recorded_surgery')
 pr$C<-model.matrix(as.formula(paste('~',paste(terms,collapse='+'))),q)[,-1,drop=FALSE]
 pt<-calculate(pr,a)
 for(g in 0:1){ii<-q$typed==g;ww<-pt$w[ii];at<-ii&q$time>=a$horizon;wa<-pt$w[at]
 diag2[[length(diag2)+1]]<-data.frame(scenario=a$name,typed=g,n=sum(ii),cancer=sum(q$status[ii]==1),other=sum(q$status[ii]==2),at_risk=sum(at),ESS=sum(ww)^2/sum(ww^2),at_risk_ESS=if(length(wa))sum(wa)^2/sum(wa^2)else 0,min_weight=min(ww),max_weight=max(ww),max_normalized_share=max(ww)/sum(ww),minG=pt$minG,maxCW=pt$maxCW,point_warnings=pt$warning_text)
 }
 chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);bs<-matrix(NA_real_,B,3);wn<-integer(B);err<-character(B)
 for(b in seq_len(B)){ix<-unlist(chunks[sample(ids,length(ids),TRUE)],use.names=FALSE);v<-tryCatch(calculate(pr,a,ix),error=function(e){err[b]<<-conditionMessage(e);NULL});if(!is.null(v)){bs[b,]<-v$vals;wn[b]<-v$warnings}}
 ok<-complete.cases(bs);ci<-apply(bs[ok,,drop=FALSE],2,quantile,c(.025,.975))
 res2[[k]]<-data.frame(scenario=a$name,estimand=names(pt$vals),estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],valid=sum(ok),warnings=sum(wn>0),B=B)
 draws2[[k]]<-data.frame(scenario=a$name,replicate=seq_len(B),control=bs[,1],detailed=bs[,2],RD=bs[,3],warnings=wn,error=err)
 write.csv(do.call(rbind,res2),file.path(out,'ipcw_reduced_risks.csv'),row.names=FALSE)
 write.csv(do.call(rbind,diag2),file.path(out,'ipcw_reduced_diagnostics.csv'),row.names=FALSE)
 write.csv(do.call(rbind,draws2),file.path(out,'ipcw_reduced_bootstrap.csv'),row.names=FALSE)
 cat('REDUCED',a$name,'RD',pt$vals[3],'valid',sum(ok),'warnings',sum(wn>0),'\n');flush.console()
}
