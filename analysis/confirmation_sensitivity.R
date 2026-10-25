# Same reference histology set; no microscopic-confirmation restriction.

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
d<-droplevels(subset(d,hist %in% codes))
set.seed(2026102401);B<-1000L
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

known<-subset(d,stage!='Unknown' & surgery!='Unknown')
res<-rbind(coding(known,'unrestricted_broad'),coding(known,'unrestricted_narrow',TRUE))
wr(res,'20_confirmation_sensitivity')
desc<-data.frame(group=c('Reference','Detailed','NOS'),N=c(nrow(d),sum(d$typed3),sum(d$nos)),surgery=c(sum(d$surgery=='Yes'),sum(d$surgery[d$typed3]=='Yes'),sum(d$surgery[d$nos]=='Yes')))
desc$surgery_pct<-100*desc$surgery/desc$N
wr(desc,'21_confirmation_composition')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
print(res);print(desc)
