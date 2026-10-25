# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE, scipen=999)
out <- 'coding_selection'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
writeout <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='NA')
sourcefile <- Sys.getenv('SEER_INPUT')
raw <- read.delim(sourcefile,colClasses='character',check.names=FALSE,quote='"',na.strings=c('','NA'),fileEncoding='UTF-8')
d <- readRDS('descriptive/broad_derived_with_ids.rds')
stopifnot(nrow(raw)==5695,identical(raw,d[,names(raw)]))
key <- function(x) paste(x[['Patient ID']],x[['Record number recode']])
stopifnot(!anyDuplicated(key(raw)),all(raw[['Primary Site - labeled']]=='C24.1-Ampulla of Vater'),
 all(raw[['Behavior recode for analysis']]=='Malignant'))
# Independently reconstruct every inclusion and modern mapping from exported labels.
hist <- as.integer(raw[['Histologic Type ICD-O-3']]); yr <- as.integer(raw[['Year of diagnosis']])
codes <- c(8140,8144,8163,8210,8211,8255,8260:8263,8480,8481,8490)
micro <- raw[['Diagnostic Confirmation']] %in% c('Positive histology','Positive exfoliative cytology, no positive histology',
 'Pos hist AND immunophenotyping AND/OR pos genetic studies','Positive microscopic confirm, method not specified')
ref <- hist %in% codes & micro; typed <- ref & hist %in% c(8144,8163,8255)
stopifnot(identical(ref,d$current),identical(typed,d$typed3),sum(ref)==4811,sum(typed)==383)
ss <- raw[['Combined Summary Stage with Expanded Regional Codes (2004+)']]
stage <- rep('Unknown',nrow(raw))
stage[ss=='Localized only'] <- 'Localized'
stage[ss %in% c('Regional by direct extension only','Regional lymph nodes involved only',
 'Regional by both direct extension and lymph node involvement')] <- 'Regional'
stage[ss=='Distant site(s)/node(s) involved'] <- 'Distant'
stopifnot(identical(stage[yr>=2004],as.character(d$stage[yr>=2004])))
sx <- suppressWarnings(as.integer(raw[['RX Summ--Surg Prim Site (1998-2022)']]))
validsx <- c(0,10:14,20:27,30,40,41,50,60,90,99)
stopifnot(all(sx[yr>=1998 & yr<=2022 & !is.na(sx)] %in% validsx))
surgery <- ifelse(yr>=1998 & yr<=2022 & !is.na(sx) & sx==0,'No',
 ifelse(yr>=1998 & yr<=2022 & !is.na(sx) & sx>=10 & sx<=90,'Yes','Unknown'))
stopifnot(identical(surgery,as.character(d$surgery)))
d$surg_code <- raw[['RX Summ--Surg Prim Site (1998-2022)']]
d$era <- ifelse(yr<=2017,'2010-2017','2018-2022')
mod <- subset(d,current & year>=2010 & year<=2022)
survmod <- subset(mod,eligible & (typed3|nos))
td <- subset(survmod,typed3 & stage=='Distant')
stopifnot(nrow(td)==91,sum(td$surgery=='Yes')==83)
writeout(data.frame(check=c('raw_records','reference','detailed','raw_mapping_mismatches','modern_reference',
 'modern_detailed_all','modern_detailed_survival','modern_detailed_distant','distant_recorded_surgery'),
 value=c(nrow(d),sum(ref),sum(typed),0,nrow(mod),sum(mod$typed3),sum(survmod$typed3),nrow(td),sum(td$surgery=='Yes'))),'audit_summary')
tab <- function(q,cols) {x<-as.data.frame(do.call(table,c(lapply(cols,function(v)q[[v]]),list(useNA='ifany'))));names(x)<-c(cols,'n');x[x$n>0,]}
writeout(tab(td,c('year','surg_code')),'distant_surgery_by_year')
writeout(tab(subset(mod,typed3),c('stage','surgery','era')),'detailed_stage_surgery_era')
writeout(tab(td,c('year','hist','surg_code','Combined Summary Stage with Expanded Regional Codes (2004+)',
 'Derived AJCC M, 7th ed (2010-2015)','Derived EOD 2018 M Recode (2018+)',
 'Derived EOD 2018 T Recode (2018+)')),'distant_crosswalk_aggregate')
mm <- subset(td,year>=2018)
writeout(tab(mm,c('Derived EOD 2018 M Recode (2018+)','Derived EOD 2018 T Recode (2018+)','surg_code')),'distant_2018plus_m_t_surgery')
stagemap <- tab(subset(d,year>=2004),c('Combined Summary Stage with Expanded Regional Codes (2004+)','stage'))
writeout(stagemap,'stage_mapping')
writeout(tab(subset(d,year>=1998 & year<=2022),c('surg_code','surgery')),'surgery_mapping')
eras <- cut(d$year,c(1974,1987,1997,2009,2017,2022,2023),labels=c('1975-1987','1988-1997','1998-2009','2010-2017','2018-2022','2023'))
hc <- do.call(rbind,lapply(levels(eras),function(e){q<-d[eras==e,]; data.frame(era=e,reference=sum(q$current),
 intestinal8144=sum(q$current&q$hist==8144),pancreatobiliary8163=sum(q$current&q$hist==8163),mixed8255=sum(q$current&q$hist==8255))}))
hc$detailed <- rowSums(hc[,3:5]);writeout(hc,'histology_by_era')
coverage <- function(q,v) do.call(rbind,lapply(unique(as.character(q[[v]])),function(lev){z<-q[as.character(q[[v]])==lev,];data.frame(variable=v,level=lev,n=nrow(z),detailed=sum(z$typed3),coverage=mean(z$typed3))}))
writeout(do.call(rbind,lapply(c('surgery','era','stage'),function(v)coverage(mod,v))),'modern_observed_coverage')
writeout(do.call(rbind,lapply(c('2010-2017','2018-2022'),function(e){z<-coverage(mod[mod$era==e,],'surgery');z$era<-e;z})),'modern_surgery_era_coverage')
writeout(tab(mod,c('era','stage','surgery','typed3')),'modern_joint_support')
writeout(data.frame(file=sourcefile,md5=unname(tools::md5sum(sourcefile))),'source_manifest')
spec <- data.frame(item=c('denominator','window','outcome','exclusions','covariates','standardization','uncertainty','interpretation'),
 value=c('All microscopically confirmed reference adenocarcinomas, including NOS and other variants','2010-2022',
 'Detailed codes 8144/8163/8255 versus all other reference codes','Unknown summary stage or unknown recorded surgery; retained in descriptive counts',
 'Surgery + era (2010-2017 / 2018-2022) + summary stage + age natural spline df3 + sex + race/origin; no interactions',
 'Set each factor level in turn and average fitted probabilities over pooled eligible covariates; paired marginal probability differences',
 '1000 patient-level percentile bootstrap; refit one logistic model; fixed spline basis and factor levels',
 'Exploratory conditional associations in the known-stage known-surgery reference population; not causal or a prediction tool'))
writeout(spec,'selection_analysis_specification')
cat('AUDIT_PASS\n');print(read.csv(file.path(out,'audit_summary.csv')));print(hc)
saveRDS(mod,file.path(out,'modern_reference_private.rds'))
# Audit assertions above must pass before fitting the published model.
set.seed(2026101403)
q <- mod[mod$stage!='Unknown' & mod$surgery!='Unknown',]
q$typed <- as.integer(q$typed3)
q$surgery <- factor(q$surgery,levels=c('No','Yes'))
q$era <- factor(q$era,levels=c('2010-2017','2018-2022'))
q$stage <- factor(q$stage,levels=c('Localized','Regional','Distant'))
agebasis<-splines::ns(q$age,df=3);q$age1<-agebasis[,1];q$age2<-agebasis[,2];q$age3<-agebasis[,3]
form <- typed ~ surgery + era + stage + age1 + age2 + age3 + sex + race
X <- model.matrix(form,q);stopifnot(qr(X)$rank==ncol(X))
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
 vals<-c(ps,ps[2]-ps[1],ps[4]-ps[3],ps[6]-ps[5],ps[7]-ps[5])
 list(vals=vals,fit=fit,warnings=paste(unique(warns),collapse='; '))
}
pt<-fitfun(seq_len(nrow(q)));stopifnot(pt$warnings=='')
chunks<-split(seq_len(nrow(q)),q$id);ids<-names(chunks);B<-bootstrap_n(1000L)
draws<-matrix(NA_real_,B,length(pt$vals));warn<-err<-character(B)
for(b in seq_len(B)){
 ix<-unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE)
 z<-tryCatch(fitfun(ix),error=function(e){err[b]<<-conditionMessage(e);NULL})
 if(!is.null(z)){draws[b,]<-z$vals;warn[b]<-z$warnings}
}
ok<-complete.cases(draws);ci<-apply(draws[ok,,drop=FALSE],2,quantile,c(.025,.975))
labels<-rbind(labs,data.frame(variable=c('surgery','era','stage','stage'),level=c('Yes minus No','2018-2022 minus 2010-2017','Regional minus Localized','Distant minus Localized')))
res<-data.frame(labels,estimate=pt$vals,lower95=ci[1,],upper95=ci[2,],estimand=c(rep('probability',nrow(labs)),rep('probability_difference',4)))
writeout(res,'adjusted_coding_probabilities')
writeout(data.frame(replicate=1:B,draws,warnings=warn,error=err),'selection_bootstrap')
writeout(data.frame(metric=c('modern_reference','model_records','model_detailed','excluded_unknown_stage','excluded_unknown_surgery','excluded_union','patients','coefficients','bootstrap_valid','bootstrap_warned','bootstrap_failed','min_fitted','max_fitted'),
 value=c(nrow(mod),nrow(q),sum(q$typed),sum(mod$stage=='Unknown'),sum(mod$surgery=='Unknown'),nrow(mod)-nrow(q),length(ids),ncol(X),sum(ok),sum(nzchar(warn)),sum(!ok),min(pt$fit$fitted.values),max(pt$fit$fitted.values))),'selection_diagnostics')
writeout(data.frame(term=colnames(X),coefficient=pt$fit$coefficients),'selection_coefficients')
saveRDS(list(q=q,formula=form,fit=pt$fit,agebasis=agebasis,results=res),file.path(out,'selection_fit_private.rds'))
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
print(res);cat('SELECTION_COMPLETE valid=',sum(ok),' warnings=',sum(nzchar(warn)),' failed=',sum(!ok),'\n')
