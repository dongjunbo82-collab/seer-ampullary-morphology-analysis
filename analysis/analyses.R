# Supplementary code for the SEER8 and expanded SEER17 analyses.
# SEER8 coding analyses: Rscript analyses.R seer8_export.txt output_directory
# SEER17 analyses: Rscript analyses.R seer17_export.txt output_directory seer17
# Both modes require authorized SEER data. No patient records are redistributed.
# SEER17 overlaps SEER8 and is not independent geographic validation.
# Outputs are aggregate results and bootstrap estimates only.
run_seer17 <- function() {
options(stringsAsFactors=FALSE,scipen=999)
library(survival)
args<-commandArgs(TRUE)
input<-if(length(args))args[1] else stop('Supply the authorized SEER17 case-listing TXT.')
out<-if(length(args)>1)args[2] else 'seer17_results'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
wr<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='NA')
raw<-read.delim(input,colClasses='character',check.names=FALSE,quote='"',na.strings=c('','NA'),fileEncoding='UTF-8')
stopifnot(ncol(raw)==60,!anyDuplicated(paste(raw[['Patient ID']],raw[['Record number recode']])),all(raw[['Primary Site - labeled']]=='C24.1-Ampulla of Vater'),all(raw[['Behavior recode for analysis']]=='Malignant'))
d<-data.frame(id=raw[['Patient ID']],year=as.integer(raw[['Year of diagnosis']]),hist=as.integer(raw[['Histologic Type ICD-O-3']]))
d$age<-as.numeric(sub(' years','',sub('85\\+','85',raw[['Age recode with single ages and 85+']])));stopifnot(all(d$age>=20),all(d$year %in% 2018:2022))
micro<-raw[['Diagnostic Confirmation']] %in% c('Positive histology','Positive exfoliative cytology, no positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies','Positive microscopic confirm, method not specified')
codes<-c(8140,8144,8163,8210,8211,8255,8260:8263,8480,8481,8490)
d$current<-micro & d$hist %in% codes;d$typed3<-d$hist %in% c(8144,8163,8255);d$typed2<-d$hist %in% c(8144,8163);d$nos<-d$hist==8140
ss<-raw[['Combined Summary Stage with Expanded Regional Codes (2004+)']];d$stage<-'Unknown'
d$stage[ss=='Localized only']<-'Localized';d$stage[ss %in% c('Regional by direct extension only','Regional lymph nodes involved only','Regional by both direct extension and lymph node involvement')]<-'Regional';d$stage[ss=='Distant site(s)/node(s) involved']<-'Distant'
d$stage<-factor(d$stage,levels=c('Localized','Regional','Distant','Unknown'))
d$surg_code<-suppressWarnings(as.integer(raw[['RX Summ--Surg Prim Site (1998-2022)']]))
d$surgery<-factor(ifelse(d$surg_code==0,'No',ifelse(d$surg_code>=10 & d$surg_code<=90,'Yes','Unknown')),levels=c('No','Yes','Unknown'))
d$surgery[is.na(d$surgery)]<-'Unknown';d$sex<-factor(raw[['Sex']])
race<-raw[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']];race[is.na(race)|grepl('Unknown|American Indian',race)]<-'Other/unknown';d$race<-factor(race)
d$time<-suppressWarnings(as.numeric(raw[['Survival months']]));d$time[d$time<0|d$time>=9999]<-NA;d$time[d$time==0]<-0.5
vital<-raw[['Vital status recode (study cutoff used)']];css<-raw[['SEER cause-specific death classification']]
d$status<-NA_integer_;d$status[vital=='Alive']<-0;d$status[vital=='Dead' & css=='Dead (attributable to this cancer dx)']<-1;d$status[vital=='Dead' & css=='Alive or dead of other cause']<-2
d$eligible<-!is.na(d$time)&!is.na(d$status)
wr(data.frame(step=c('All malignant age20+ C241','Microscopic confirmation','Reference adenocarcinoma','Detailed broad','Detailed narrow','NOS','Reference survival eligible','Reference missing stage','Reference unknown surgery'),N=c(nrow(d),sum(micro),sum(d$current),sum(d$current&d$typed3),sum(d$current&d$typed2),sum(d$current&d$nos),sum(d$current&d$eligible),sum(d$current&d$stage=='Unknown'),sum(d$current&d$surgery=='Unknown'))),'01_cohort_flow')
wr(data.frame(input=basename(input),md5=unname(tools::md5sum(input)),database='SEER Research 17 Nov2025',diagnosis='2018-2022',seed=2026102307,B=1000),'00_provenance')
d<-droplevels(subset(d,current))

# Descriptive sources for Figure 2b and Table S13B; no random resampling.
annual <- do.call(rbind,lapply(2018:2022,function(y) do.call(rbind,lapply(c('Yes','No','Unknown','All'),function(g){
 z<-d[d$year==y & (g=='All' | d$surgery==g),];n<-nrow(z);k<-sum(z$typed3)
 ci<-if(n>0)binom.test(k,n)$conf.int else c(NA,NA)
 data.frame(year=y,surgery=g,N=n,detailed=k,pct=100*k/n,lower=100*ci[1],upper=100*ci[2])
}))))
wr(annual,'18_annual_surgery_coverage')
missing <- do.call(rbind,lapply(c('Detailed','NOS','Other reference','Full reference'),function(g){
 z<-switch(g,Detailed=d[d$typed3,],NOS=d[d$nos,],'Other reference'=d[!d$typed3 & !d$nos,],'Full reference'=d)
 us<-z$stage=='Unknown';ut<-z$surgery=='Unknown';it<-is.na(z$time);ie<-is.na(z$status)
 data.frame(group=g,N=nrow(z),unknown_stage=sum(us),unknown_surgery=sum(ut),both_unknown=sum(us&ut),excluded_coding=sum(us|ut),invalid_survival_time=sum(it),unclassified_event=sum(ie),excluded_survival=sum(!z$eligible))
}))
wr(missing,'19_missingness_by_coding_group')

support<-do.call(rbind,lapply(c(8144,8163),function(h)do.call(rbind,lapply(2018:2022,function(y){z<-subset(d,hist==h&year==y&eligible);data.frame(hist=h,year=y,N=nrow(z),cancer_events=sum(z$status==1),other_events=sum(z$status==2),at12=sum(z$time>=12),at36=sum(z$time>=36),at60=sum(z$time>=60),max_months=max(z$time))}))))
wr(support,'02_subtype_followup_support');print(support)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
wr(data.frame(analysis=c('Coding broad primary','Coding narrow','Coding omit stage','Coding unknown stage category','Subtype 3year primary','Subtype 1year secondary','Subtype surgery 3year sensitivity','Subtype strict resection 3year sensitivity','Detailed versus NOS 3year secondary'),population=c(rep('Reference 2018-2022',4),'8144 and 8163, 2018-2020','8144 and 8163, 2018-2022','8144 and 8163 surgery, 2018-2020','8144 and 8163 codes30/40/60, 2018-2020','Detailed broad and NOS, 2018-2020'),method=c(rep('Logistic marginal standardization',4),rep('Overlap weighted Aalen-Johansen',5)),adjustment=c('Surgery + year linear + stage + age ns3 + sex + race','Same, narrow outcome retains8255 in reference','Same except no stage','Same, Unknown as stage category','Year linear + age ns3 + sex + race + stage + surgery','Same','Same except surgery constant','Same except surgery constant','Year linear + age ns3 + sex + race, stage; with/without surgery'),B=1000),'03_analysis_specification')
coverage<-function(q,v)do.call(rbind,lapply(unique(as.character(q[[v]])),function(l){z<-q[as.character(q[[v]])==l,];data.frame(variable=v,level=l,N=nrow(z),broad=sum(z$typed3),narrow=sum(z$typed2),broad_pct=100*mean(z$typed3),narrow_pct=100*mean(z$typed2))}))
wr(do.call(rbind,lapply(c('surgery','stage','year'),function(v)coverage(d,v))),'04_observed_coverage')
wr(as.data.frame(with(d,table(year,hist))),'05_annual_histology')
comp<-function(z){c(N=nrow(z),mean_age=mean(z$age),sd_age=sd(z$age),mean_year=mean(z$year),female_pct=100*mean(z$sex=='Female'),localized_pct=100*mean(z$stage=='Localized'),regional_pct=100*mean(z$stage=='Regional'),distant_pct=100*mean(z$stage=='Distant'),unknown_stage_pct=100*mean(z$stage=='Unknown'),surgery_pct=100*mean(z$surgery=='Yes'),unknown_surgery_pct=100*mean(z$surgery=='Unknown'))}
groups<-list(Reference=d,Detailed=subset(d,typed3),NOS=subset(d,nos),Intestinal8144=subset(d,hist==8144),Pancreatobiliary8163=subset(d,hist==8163))
wr(data.frame(characteristic=names(comp(d)),sapply(groups,comp),row.names=NULL),'06_composition_all_reference')
wr(as.data.frame(with(subset(d,typed3|nos),table(group=ifelse(typed3,'Detailed','NOS'),surg_code))),'07_procedure_composition')
set.seed(2026102307);B<-1000L
bootix<-function(q){chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);function()unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE)}
coding<-function(q,scenario,narrow=FALSE,omitstage=FALSE){
 q<-droplevels(q);q$g<-as.integer(if(narrow)q$typed2 else q$typed3);q$surgery<-factor(q$surgery,levels=c('No','Yes'))
 ab<-splines::ns(q$age,3);q$a1<-ab[,1];q$a2<-ab[,2];q$a3<-ab[,3]
 fm<-as.formula(paste('g~surgery+year+',if(!omitstage)'stage+' else '', 'a1+a2+a3+sex+race'))
 X<-model.matrix(fm,q);stopifnot(qr(X)$rank==ncol(X))
 mats<-lapply(c('No','Yes'),function(s){z<-q;z$surgery<-factor(s,levels=c('No','Yes'));model.matrix(fm,z)})
 calc<-function(ix){warn<-character();fit<-withCallingHandlers(glm.fit(X[ix,,drop=FALSE],q$g[ix],family=binomial(),control=glm.control(maxit=100)),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')});if(!fit$converged||any(!is.finite(fit$coefficients)))stop('failed coding model');v<-sapply(mats,function(m)mean(plogis(m[ix,,drop=FALSE]%*%fit$coefficients)));list(v=c(v,v[2]-v[1]),warn=paste(unique(warn),collapse=';'))}
 pt<-calc(seq_len(nrow(q)));draw<-bootix(q);bs<-matrix(NA_real_,B,3);warn<-err<-character(B)
 for(b in seq_len(B)){z<-tryCatch(calc(draw()),error=function(e){err[b]<<-conditionMessage(e);NULL});if(!is.null(z)){bs[b,]<-z$v;warn[b]<-z$warn}}
 ok<-complete.cases(bs);ci<-apply(bs[ok,,drop=FALSE],2,quantile,c(.025,.975));wr(data.frame(replicate=1:B,bs,warnings=warn,error=err),paste0('bootstrap_coding_',scenario))
 # Independent formula glm/predict point calculation.
 ff<-glm(fm,data=q,family=binomial());pr<-sapply(c('No','Yes'),function(s){z<-q;z$surgery<-factor(s,levels=c('No','Yes'));mean(predict(ff,z,type='response'))});stopifnot(max(abs(pr-pt$v[1:2]))<1e-9)
 cat('CODING',scenario,'N',nrow(q),'RD',pt$v[3],'valid',sum(ok),'\n');flush.console()
 data.frame(scenario,estimand=c('No surgery','Surgery','Yes minus No'),estimate=pt$v,lower95=ci[1,],upper95=ci[2,],N=nrow(q),detailed=sum(q$g),B,valid=sum(ok),warned=sum(nzchar(warn)),failed=sum(!ok),point_warning=pt$warn,independent_max_error=max(abs(pr-pt$v[1:2])))
}
known<-subset(d,stage!='Unknown' & surgery!='Unknown');ks<-subset(d,surgery!='Unknown')
cr<-rbind(coding(known,'broad'),coding(known,'narrow',TRUE),coding(ks,'omit_stage',omitstage=TRUE),coding(ks,'unknown_stage_category'))
wr(cr,'08_adjusted_coding_probabilities');print(cr)
iq<-droplevels(known);iq$g<-as.integer(iq$typed3)
f0<-glm(g~surgery+year+stage+splines::ns(age,3)+sex+race,data=iq,family=binomial());f1<-update(f0,.~.+surgery:year);ia<-anova(f0,f1,test='LRT')
wr(data.frame(test='surgery by linear diagnosis year coding interaction',N=nrow(iq),df=ia$Df[2],statistic=ia$Deviance[2],p=ia$`Pr(>Chi)`[2]),'09_coding_year_interaction')
aj<-function(time,status,weight,horizon){
 keep<-time<=horizon;tt<-sort(unique(time[keep]));ss<-1;ff<-0;ix<-match(time,tt);at<-sum(weight)
 sums<-function(mask){v<-numeric(length(tt));z<-tapply(weight[mask],ix[mask],sum);if(length(z))v[as.integer(names(z))]<-z;v}
 exits<-sums(keep);e1<-sums(keep&status==1);ea<-sums(keep&status!=0)
 for(i in seq_along(tt)){if(at>1e-12){ff<-ff+ss*e1[i]/at;ss<-ss*(1-ea[i]/at)};at<-at-exits[i]};c(cancer=ff,allcause=1-ss)
}
sc<-data.frame(name=c('subtype_36_primary','subtype_12_full','subtype_36_surgery','subtype_36_strict','detailed_NOS_36_baseline','detailed_NOS_36_surgery_balanced','detailed_NOS_36_surgery'),subtype=c(rep(TRUE,4),rep(FALSE,3)),end=c(2020,2022,2020,2020,2020,2020,2020),horizon=c(36,12,36,36,36,36,36),surg=c('All','All','Yes','Strict','All','All','Yes'),balance_surgery=c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE))
wr(sc,'10_survival_scenarios')
risks<-unweighted<-counts<-balances<-checks<-censoring<-list()
for(k in seq_len(nrow(sc))){
 s<-sc[k,];full<-subset(d,eligible & year<=s$end & if(s$subtype)hist %in% c(8144,8163) else typed3|nos)
 if(s$surg=='Yes')full<-subset(full,surgery=='Yes');if(s$surg=='Strict')full<-subset(full,surg_code %in% c(30,40,60))
 full$g<-if(s$subtype)as.integer(full$hist==8144) else as.integer(full$typed3)
 # Unweighted comparison retains unknown stage/surgery; adjusted comparison explicitly complete cases.
 for(g in 0:1){z<-full[full$g==g,];v<-aj(z$time,z$status,rep(1,nrow(z)),s$horizon);unweighted[[length(unweighted)+1]]<-data.frame(scenario=s$name,g,N=nrow(z),cancer=v[1],allcause=v[2])}
 q<-droplevels(subset(full,stage!='Unknown' & surgery!='Unknown'))
 ab<-splines::ns(q$age,3);q$a1<-ab[,1];q$a2<-ab[,2];q$a3<-ab[,3]
 terms<-c('year','a1','a2','a3','sex','race','stage');if(s$balance_surgery)terms<-c(terms,'surgery')
 terms<-terms[sapply(terms,function(t)length(unique(q[[t]]))>1)]
 X<-model.matrix(as.formula(paste('~',paste(terms,collapse='+'))),q);stopifnot(qr(X)$rank==ncol(X))
 calc<-function(ix,weighted=TRUE){
 z<-q[ix,];warn<-character();fit<-withCallingHandlers(glm.fit(X[ix,,drop=FALSE],z$g,family=binomial(),control=glm.control(maxit=100)),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
 # Aliased coefficients caused by an absent bootstrap category do not invalidate
 # fitted propensities for the same resample. No out-of-sample prediction is used.
 if(!fit$converged||any(!is.finite(fit$fitted.values))||length(unique(z$g))<2)stop('failed propensity fit')
 p<-fit$fitted.values;wt<-ifelse(z$g==1,1-p,p)
 vals<-sapply(0:1,function(g){a<-z$g==g;if(sum(wt[a])<1e-8||!any(z$time[a]>=s$horizon))stop('no horizon support');aj(z$time[a],z$status[a],wt[a],s$horizon)})
 list(v=c(vals[1,],vals[1,2]-vals[1,1],vals[2,],vals[2,2]-vals[2,1]),wt=wt,p=p,warn=paste(unique(warn),collapse=';'))
 }
 pt<-calc(seq_len(nrow(q)));draw<-bootix(q);bs<-matrix(NA_real_,B,6);warn<-err<-character(B)
 for(b in 1:B){z<-tryCatch(calc(draw()),error=function(e){err[b]<<-conditionMessage(e);NULL});if(!is.null(z)){bs[b,]<-z$v;warn[b]<-z$warn}}
 ok<-complete.cases(bs);ci<-apply(bs[ok,,drop=FALSE],2,quantile,c(.025,.975))
 risks[[k]]<-data.frame(scenario=s$name,horizon=s$horizon,estimand=c('cancer_g0','cancer_g1','cancer_RD_g1_minus_g0','allcause_g0','allcause_g1','allcause_RD_g1_minus_g0'),estimate=pt$v,lower95=ci[1,],upper95=ci[2,],N=nrow(q),excluded_unknown=nrow(full)-nrow(q),valid=sum(ok),warned=sum(nzchar(warn)),failed=sum(!ok),point_warning=pt$warn)
 wr(data.frame(replicate=1:B,bs,warnings=warn,error=err),paste0('bootstrap_survival_',s$name))
 for(g in 0:1){
  a<-q$g==g;z<-q[a,];z$wt<-pt$wt[a];z$ev<-factor(z$status,levels=0:2,labels=c('censor','cancer','other'));z$record<-seq_len(nrow(z));sf<-survfit(Surv(time,ev)~1,data=z,weights=wt,id=record);j<-max(which(sf$time<=s$horizon));chk<-sf$pstate[j,match('cancer',sf$states)];stopifnot(abs(chk-pt$v[g+1])<1e-9)
  km<-survfit(Surv(time,status!=0)~1,data=z,weights=wt);j2<-max(which(km$time<=s$horizon));stopifnot(abs(1-km$surv[j2]-pt$v[g+4])<1e-9)
  checks[[length(checks)+1]]<-data.frame(scenario=s$name,g,AJ_error=abs(chk-pt$v[g+1]),KM_error=abs(1-km$surv[j2]-pt$v[g+4]))
  counts[[length(counts)+1]]<-data.frame(scenario=s$name,g,N=nrow(z),cancer_before_horizon=sum(z$status==1&z$time<=s$horizon),other_before_horizon=sum(z$status==2&z$time<=s$horizon),censored_before_horizon=sum(z$status==0&z$time<s$horizon),at_risk=sum(z$time>=s$horizon),weighted_at_risk=sum(z$wt[z$time>=s$horizon]),ESS=sum(z$wt)^2/sum(z$wt^2),min_propensity=min(pt$p[a]),max_propensity=max(pt$p[a]),surgery_pct=100*mean(z$surgery=='Yes'))
  rev<-survfit(Surv(time,status==0)~1,data=z);j3<-which(rev$time<=s$horizon);cv<-if(length(j3))rev$surv[max(j3)] else 1
  censoring[[length(censoring)+1]]<-data.frame(scenario=s$name,g,reverse_KM_probability_followup_at_horizon=cv)
 }
 for(j in seq_len(ncol(X))[-1]){a<-q$g==1;sd0<-sqrt((var(X[a,j])+var(X[!a,j]))/2);balances[[length(balances)+1]]<-data.frame(scenario=s$name,term=colnames(X)[j],SMD_before=if(sd0>0)(mean(X[a,j])-mean(X[!a,j]))/sd0 else NA,SMD_after=if(sd0>0)(weighted.mean(X[a,j],pt$wt[a])-weighted.mean(X[!a,j],pt$wt[!a]))/sd0 else NA)}
 cat('SURVIVAL',s$name,'N',nrow(q),'RD',pt$v[3],'valid',sum(ok),'\n');flush.console()
}
wr(do.call(rbind,risks),'11_overlap_weighted_risks');wr(do.call(rbind,unweighted),'12_unweighted_risks');wr(do.call(rbind,counts),'13_survival_support_ESS');wr(do.call(rbind,balances),'14_covariate_balance');wr(do.call(rbind,checks),'15_independent_AJ_KM_check');wr(do.call(rbind,censoring),'16_censoring_support')
# Subtype Cox HR is not reported: the exploratory full-covariate fit had
# an infinite race coefficient (Other/unknown has zero cancer events).
# The prespecified overlap-weighted absolute risks above do not require PH.
q<-subset(d,eligible & hist %in% c(8144,8163)&year<=2020&stage!='Unknown'&surgery!='Unknown')
wr(as.data.frame(with(q,table(race,status))),'17_subtype_race_event_support')
cat('ANALYSIS COMPLETE\n')
}
mode_args<-commandArgs(TRUE)
if(length(mode_args)>=3 && mode_args[3]=='seer17'){run_seer17();quit(save='no',status=0)}

# SEER8 coding and sensitivity analyses follow.
# Earlier statistical code archive (v1.0.0): https://doi.org/10.5281/zenodo.23374749
# Supplementary Code: coding audit and modern coding-probability analysis.
# R 4.6.1; base R and recommended splines only. No installation is performed.
# Usage: Rscript analyses.R authorized_case_listing.txt output_directory
# Input: UTF-8, tab-separated SEER*Stat export with the original variable names.
# Extraction: SEER Research Data, 8 Registries, Nov 2025 Sub (1975-2023);
# malignant C24.1, known age >=20, all histologies and diagnostic confirmations.
# Do not restrict to first primary or to complete survival. Do not redistribute records.
# Scope: reproduces the new coding audit/probability results, not historical survival models.
options(stringsAsFactors=FALSE,scipen=999)
args<-commandArgs(TRUE)
if(length(args)<1)stop('Supply an authorized case-listing export path.')
out<-if(length(args)>=2)args[2] else 'coding_selection_results'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
wr<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='NA')
raw<-read.delim(args[1],colClasses='character',check.names=FALSE,quote='"',na.strings=c('','NA'),fileEncoding='UTF-8')
key<-paste(raw[['Patient ID']],raw[['Record number recode']])
stopifnot(!anyDuplicated(key),all(raw[['Primary Site - labeled']]=='C24.1-Ampulla of Vater'),
 all(raw[['Behavior recode for analysis']]=='Malignant'))
d<-data.frame(id=raw[['Patient ID']],year=as.integer(raw[['Year of diagnosis']]),hist=as.integer(raw[['Histologic Type ICD-O-3']]))
d$age<-as.numeric(sub(' years','',sub('85\\+','85',raw[['Age recode with single ages and 85+']])))
stopifnot(all(d$age>=20))
micro<-raw[['Diagnostic Confirmation']] %in% c('Positive histology','Positive exfoliative cytology, no positive histology',
 'Pos hist AND immunophenotyping AND/OR pos genetic studies','Positive microscopic confirm, method not specified')
codes<-c(8140,8144,8163,8210,8211,8255,8260:8263,8480,8481,8490)
d$current<-micro & d$hist %in% codes
d$typed3<-d$current & d$hist %in% c(8144,8163,8255)
ss<-raw[['Combined Summary Stage with Expanded Regional Codes (2004+)']]
d$stage<-'Unknown'
d$stage[ss=='Localized only']<-'Localized'
d$stage[ss %in% c('Regional by direct extension only','Regional lymph nodes involved only',
 'Regional by both direct extension and lymph node involvement')]<-'Regional'
d$stage[ss=='Distant site(s)/node(s) involved']<-'Distant'
sx<-suppressWarnings(as.integer(raw[['RX Summ--Surg Prim Site (1998-2022)']]))
d$surgery<-ifelse(d$year>=1998 & d$year<=2022 & !is.na(sx) & sx==0,'No',
 ifelse(d$year>=1998 & d$year<=2022 & !is.na(sx) & sx>=10 & sx<=90,'Yes','Unknown'))
d$sex<-factor(raw[['Sex']]);race<-raw[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']]
race[is.na(race)|grepl('Unknown|American Indian',race)]<-'Other/unknown';d$race<-factor(race)
d$era<-ifelse(d$year<=2017,'2010-2017','2018-2022')
tab<-function(q,vars){x<-as.data.frame(do.call(table,c(lapply(vars,function(v)q[[v]]),list(useNA='ifany'))));names(x)<-c(vars,'n');x[x$n>0,]}
era<-cut(d$year,c(1974,1987,1997,2009,2017,2022,2023),labels=c('1975-1987','1988-1997','1998-2009','2010-2017','2018-2022','2023'))
wr(do.call(rbind,lapply(levels(era),function(e){z<-d[era==e,];data.frame(era=e,reference=sum(z$current),
 intestinal8144=sum(z$current & z$hist==8144),pancreatobiliary8163=sum(z$current & z$hist==8163),mixed8255=sum(z$current & z$hist==8255))})),'histology_by_era')
mod<-subset(d,current & year>=2010 & year<=2022)
coverage<-function(q,v)do.call(rbind,lapply(unique(as.character(q[[v]])),function(lev){z<-q[as.character(q[[v]])==lev,];data.frame(variable=v,level=lev,n=nrow(z),detailed=sum(z$typed3),coverage=mean(z$typed3))}))
wr(do.call(rbind,lapply(c('surgery','era','stage'),function(v)coverage(mod,v))),'modern_observed_coverage')
wr(do.call(rbind,lapply(c('2010-2017','2018-2022'),function(e){z<-coverage(mod[mod$era==e,],'surgery');z$era<-e;z})),'modern_surgery_era_coverage')
wr(tab(mod,c('era','stage','surgery','typed3')),'modern_joint_support')
# Reproduce the stage/surgery audit in the survival-eligible detailed subset.
tm<-suppressWarnings(as.numeric(raw[['Survival months']]))
vital<-raw[['Vital status recode (study cutoff used)']];css<-raw[['SEER cause-specific death classification']]
eventok<-(!is.na(vital)&vital=='Alive')|(!is.na(vital)&vital=='Dead'&!is.na(css)&css %in% c('Dead (attributable to this cancer dx)','Alive or dead of other cause'))
eligible<-!is.na(tm)&tm>=0&tm<9999&eventok
td<-raw[d$typed3 & d$year>=2010 & d$year<=2022 & d$stage=='Distant' & eligible,]
wr(tab(td,c('Year of diagnosis','RX Summ--Surg Prim Site (1998-2022)')),'distant_surgery_by_year')
wr(tab(td,c('Year of diagnosis','Derived EOD 2018 M Recode (2018+)','Derived EOD 2018 T Recode (2018+)','RX Summ--Surg Prim Site (1998-2022)')),'distant_m_t_surgery')
q<-mod[mod$stage!='Unknown' & mod$surgery!='Unknown',]
q$typed<-as.integer(q$typed3);q$surgery<-factor(q$surgery,levels=c('No','Yes'))
q$era<-factor(q$era,levels=c('2010-2017','2018-2022'));q$stage<-factor(q$stage,levels=c('Localized','Regional','Distant'))
ab<-splines::ns(q$age,df=3);q$age1<-ab[,1];q$age2<-ab[,2];q$age3<-ab[,3]
form<-typed ~ surgery + era + stage + age1 + age2 + age3 + sex + race
X<-model.matrix(form,q);stopifnot(qr(X)$rank==ncol(X))
settings<-list();labs<-data.frame(variable=character(),level=character())
for(v in c('surgery','era','stage'))for(lev in levels(q[[v]])){
 z<-q;z[[v]]<-factor(rep(lev,nrow(z)),levels=levels(q[[v]]))
 settings[[length(settings)+1]]<-model.matrix(form,z)
 labs<-rbind(labs,data.frame(variable=v,level=lev))
}
fitfun<-function(ix){
 warns<-character();fit<-withCallingHandlers(glm.fit(X[ix,,drop=FALSE],q$typed[ix],family=binomial(),control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 if(!fit$converged||any(!is.finite(fit$coefficients)))stop('Nonconverged or aliased fit')
 ps<-sapply(settings,function(xx)mean(plogis(xx[ix,,drop=FALSE]%*%fit$coefficients)))
 list(vals=c(ps,ps[2]-ps[1],ps[4]-ps[3],ps[6]-ps[5],ps[7]-ps[5]),fit=fit,warnings=paste(unique(warns),collapse='; '))
}
pt<-fitfun(seq_len(nrow(q)))
set.seed(2026101403)
chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);B<-1000L
draws<-matrix(NA_real_,B,length(pt$vals));warn<-err<-character(B)
for(b in seq_len(B)){
 ix<-unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE)
 z<-tryCatch(fitfun(ix),error=function(e){err[b]<<-conditionMessage(e);NULL})
 if(!is.null(z)){draws[b,]<-z$vals;warn[b]<-z$warnings}
}
ok<-complete.cases(draws);ci<-apply(draws[ok,,drop=FALSE],2,quantile,c(.025,.975))
labels<-rbind(labs,data.frame(variable=c('surgery','era','stage','stage'),level=c('Yes minus No','2018-2022 minus 2010-2017','Regional minus Localized','Distant minus Localized')))
res<-data.frame(labels,estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],estimand=c(rep('probability',nrow(labs)),rep('probability_difference',4)))
wr(res,'adjusted_coding_probabilities')
wr(data.frame(replicate=1:B,draws,warnings=warn,error=err),'selection_bootstrap')
wr(data.frame(metric=c('modern_reference','model_records','model_detailed','excluded_unknown_stage','excluded_unknown_surgery','excluded_union','patients','coefficients','bootstrap_valid','bootstrap_warned','bootstrap_failed'),
 value=c(nrow(mod),nrow(q),sum(q$typed),sum(mod$stage=='Unknown'),sum(mod$surgery=='Unknown'),nrow(mod)-nrow(q),length(ids),ncol(X),sum(ok),sum(nzchar(warn)),sum(!ok))),'selection_diagnostics')
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
print(res)


# Reuse verified reference records constructed from the authorized export above.
d<-mod
d$typed2<-d$hist %in% c(8144L,8163L)
d$stage<-factor(d$stage,levels=c("Distant","Localized","Regional","Unknown"))
d$surgery<-factor(d$surgery,levels=c("No","Unknown","Yes"))
stopifnot(nrow(d)==1904,sum(d$typed3)==289,sum(d$year>=2018)==803)
set.seed(2026102301)
writeout<-function(x,n)wr(x,paste0('seer8_',n))
writeout(data.frame(scenario=c('modern_broad','modern_narrow','all_no_stage'),window=c('2018-2022','2018-2022','2010-2022'),outcome=c('8144/8163/8255','8144/8163','8144/8163/8255'),adjustment=c('surgery + linear year + stage + age ns3 + sex + race','surgery + linear year + stage + age ns3 + sex + race','surgery + era + age ns3 + sex + race'),B=1000),'specification')
run<-function(q,narrow=FALSE,omitstage=FALSE){
 q<-droplevels(q);q$g<-as.integer(if(narrow)q$typed2 else q$typed3)
 q$surgery<-factor(q$surgery,levels=c('No','Yes'))
 ab<-splines::ns(q$age,df=3);q$a1<-ab[,1];q$a2<-ab[,2];q$a3<-ab[,3]
 era<-length(unique(q$year))>5
 fm<-as.formula(paste('g ~ surgery +',if(era)'factor(era)' else 'year','+',if(!omitstage)'stage +' else '', 'a1+a2+a3+sex+race'))
 X<-model.matrix(fm,q);stopifnot(qr(X)$rank==ncol(X))
 mats<-lapply(c('No','Yes'),function(s){z<-q;z$surgery<-factor(s,levels=c('No','Yes'));model.matrix(fm,z)})
 calc<-function(ix){warn<-character();f<-withCallingHandlers(glm.fit(X[ix,,drop=FALSE],q$g[ix],family=binomial(),control=glm.control(maxit=100)),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')});if(!f$converged||any(!is.finite(f$coefficients)))stop('nonfinite fit');pr<-sapply(mats,function(m)mean(plogis(m[ix,,drop=FALSE]%*%f$coefficients)));list(v=c(pr,pr[2]-pr[1]),warn=length(warn)>0)}
 pt<-calc(seq_len(nrow(q)));chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks)
 out<-matrix(NA_real_,1000,3);warn<-logical(1000)
 for(b in 1:1000){ix<-unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE);z<-tryCatch(calc(ix),error=function(e)NULL);if(!is.null(z)){out[b,]<-z$v;warn[b]<-z$warn}}
 ci<-apply(out,2,quantile,c(.025,.975),na.rm=TRUE)
 data.frame(estimand=c('No surgery','Surgery','Difference'),estimate=pt$v,lower=ci[1,],upper=ci[2,],N=nrow(q),detailed=sum(q$g),valid=sum(complete.cases(out)),warned=sum(warn))
}
modern<-subset(d,year>=2018 & stage!='Unknown' & surgery!='Unknown')
allknown<-subset(d,surgery!='Unknown')
r<-rbind(transform(run(modern),scenario='modern_broad'),transform(run(modern,TRUE),scenario='modern_narrow'),transform(run(allknown,omitstage=TRUE),scenario='all_no_stage'))
writeout(r,'standardized');print(r)
q<-droplevels(subset(d,stage!='Unknown' & surgery!='Unknown'));q$g<-as.integer(q$typed3)
f0<-glm(g~surgery+factor(era)+stage+splines::ns(age,3)+sex+race,data=q,family=binomial())
f1<-update(f0,.~.+surgery:factor(era));a<-anova(f0,f1,test='LRT')
writeout(data.frame(test='surgery by era interaction',N=nrow(q),df=a$Df[2],statistic=a$Deviance[2],p=a$`Pr(>Chi)`[2]),'interaction');print(a)

ann<-do.call(rbind,lapply(2010:2022,function(y)do.call(rbind,lapply(c("No","Yes"),function(s){z<-subset(d,year==y&surgery==s);data.frame(year=y,surgery=s,N=nrow(z),broad=sum(z$typed3),narrow=sum(z$typed2))}))))
writeout(ann,"annual_coverage")

# Reconstruct survival eligibility only for the added descriptive composition tables.
rawyear<-as.integer(raw[['Year of diagnosis']])
rawhist<-as.integer(raw[['Histologic Type ICD-O-3']])
ix<-which(rawyear>=2010 & rawyear<=2022 & micro & rawhist %in% codes)
stopifnot(length(ix)==nrow(d))
tm<-suppressWarnings(as.numeric(raw[['Survival months']]))
tm[tm>=9999 | tm<0]<-NA_real_
css<-raw[['SEER cause-specific death classification']]
vital<-raw[['Vital status recode (study cutoff used)']]
st<-rep(NA_integer_,nrow(raw))
st[which(vital=='Alive')]<-0L
st[which(vital=='Dead' & css=='Dead (attributable to this cancer dx)')]<-1L
st[which(vital=='Dead' & css=='Alive or dead of other cause')]<-2L
d$eligible<-(!is.na(tm)&!is.na(st))[ix]
d$nos<-d$hist==8140
d$surg_code<-raw[['RX Summ--Surg Prim Site (1998-2022)']][ix]
ref<-subset(d,eligible);nos<-subset(ref,nos);det<-subset(ref,typed3)
summary<-function(z){n<-nrow(z);c(N=n,age_mean=mean(z$age),age_sd=sd(z$age),year_mean=mean(z$year),female=sum(z$sex=='Female'),white=sum(z$race=='Non-Hispanic White'),black=sum(z$race=='Non-Hispanic Black'),asian=sum(z$race=='Non-Hispanic Asian or Pacific Islander'),hispanic=sum(z$race=='Hispanic (All Races)'),other=sum(z$race=='Other/unknown'),localized=sum(z$stage=='Localized'),regional=sum(z$stage=='Regional'),distant=sum(z$stage=='Distant'),unknownstage=sum(z$stage=='Unknown'),surgery=sum(z$surgery=='Yes'),nosurgery=sum(z$surgery=='No'),unknownsurgery=sum(z$surgery=='Unknown'))}
comp<-data.frame(characteristic=names(summary(ref)),reference=summary(ref),detailed=summary(det),NOS=summary(nos),row.names=NULL)
writeout(comp,'composition');print(comp)
proc<-as.data.frame(with(subset(ref,typed3|nos),table(group=ifelse(typed3,'Detailed','NOS'),code=surg_code)));writeout(proc,'procedure_composition')
