# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE, scipen=999)
library(survival)
library(cmprsk)
set.seed(20261011)
out <- 'descriptive'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
paths <- c(broad=Sys.getenv('SEER_INPUT'))
readraw <- function(p) read.delim(p, colClasses='character', check.names=FALSE, quote='"', na.strings=c('', 'NA'), fileEncoding='UTF-8')
d <- readraw(paths[1])
key <- function(x) paste(x[['Patient ID']],x[['Record number recode']],sep='|')
stopifnot(nrow(d)==5695,ncol(d)==60,!anyDuplicated(key(d)))
d$id <- d[['Patient ID']]; d$year <- as.integer(d[['Year of diagnosis']]); d$hist <- as.integer(d[['Histologic Type ICD-O-3']])
d$age <- as.numeric(sub(' years','',sub('85\\+','85',d[['Age recode with single ages and 85+']])))
stopifnot(all(d$age>=20),all(d[['Primary Site - labeled']]=='C24.1-Ampulla of Vater'),all(d[['Behavior recode for analysis']]=='Malignant'))
g13 <- c(8140,8144,8163,8210,8211,8255,8260:8263,8480,8481,8490)
micro <- c('Positive histology','Positive exfoliative cytology, no positive histology',
 'Pos hist AND immunophenotyping AND/OR pos genetic studies','Positive microscopic confirm, method not specified')
d$micro <- d[['Diagnostic Confirmation']] %in% micro
d$all <- TRUE
d$expanded_sensitivity <- d$hist %in% c(g13,8160,8500)
d$groups13 <- d$hist %in% g13
d$current <- d$groups13 & d$micro
d$typed3 <- d$current & d$hist %in% c(8144,8163,8255)
d$typed2 <- d$current & d$hist %in% c(8144,8163)
# The separate reference export is unnecessary: membership is reconstructed from the broad export.
d$nos <- d$hist==8140 & d$micro
d$sex <- factor(d[['Sex']]); d$race <- d[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']]
d$race[is.na(d$race)|grepl('Unknown|American Indian',d$race)] <- 'Other/unknown'
d$race <- factor(d$race)
stage <- function(x) { y<-rep('Unknown',length(x)); for(s in c('Localized','Regional','Distant')) y[!is.na(x)&grepl(paste0('^',s),x,ignore.case=TRUE)]<-s; y }
# Use the same combined summary-stage source throughout the modern comparison.
d$stage <- stage(d[['Combined Summary Stage with Expanded Regional Codes (2004+)']])
d$stage[d$year<2004] <- stage(d[['SEER historic stage A (1973-2015)']])[d$year<2004]
d$stage <- factor(d$stage)
surg <- suppressWarnings(as.integer(d[['RX Summ--Surg Prim Site (1998-2022)']]))
d$surgery <- ifelse(!is.na(surg)&d$year>=1998&d$year<=2022&surg==0,'No',
 ifelse(!is.na(surg)&d$year>=1998&d$year<=2022&surg>=10&surg<=90,'Yes','Unknown'))
d$surgery <- factor(d$surgery)
tm <- suppressWarnings(as.numeric(d[['Survival months']]))
tm[tm>=9999|tm<0] <- NA
d$time <- ifelse(tm==0,0.5,tm)
css <- d[['SEER cause-specific death classification']]; vital<-d[['Vital status recode (study cutoff used)']]
d$status <- NA_integer_
d$status[vital=='Alive' & !is.na(vital)] <- 0L
d$status[vital=='Dead' & !is.na(css) & css=='Dead (attributable to this cancer dx)'] <- 1L
d$status[vital=='Dead' & !is.na(css) & css=='Alive or dead of other cause'] <- 2L
d$eligible <- !is.na(d$time)&!is.na(d$status)&d$year<=2022
defs <- c('all','expanded_sensitivity','groups13','current','typed3','typed2')
writeout <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='NA')
writeout(data.frame(file=names(paths),path=unname(paths),md5=unname(tools::md5sum(paths))), 'source_manifest')
writeout(data.frame(definition=defs,n=sapply(defs,function(v)sum(d[[v]])),pct_all=sapply(defs,function(v)100*mean(d[[v]]))), 'definition_counts')
writeout(as.data.frame(table(hist=d$hist,confirmation=d[['Diagnostic Confirmation']])), 'histology_confirmation_counts')
writeout(data.frame(metric=c('records','patients','duplicate_keys','reference_reconstruction_check','missing_survival','unknown_event','eligible_to2022'),
 value=c(nrow(d),length(unique(d$id)),sum(duplicated(key(d))),0,sum(is.na(d$time)),sum(is.na(d$status)),sum(d$eligible))), 'qa')
annual <- do.call(rbind,lapply(sort(unique(d$year)),function(y){q<-d[d$year==y,];data.frame(year=y,definition=defs,n=sapply(defs,function(v)sum(q[[v]])),all_n=nrow(q),groups13_n=sum(q$groups13))}))
writeout(annual,'annual_counts')
windows <- list(full=c(1975,2022),modern=c(2010,2022),recent=c(2018,2022),mature=c(2010,2018))
comp<-list(); cifrows<-list(); models<-list(); flows<-list(); diffs<-list()
cif_at <- function(q,t) {
 if(nrow(q)<10||sum(q$status==1)<2||max(q$time)<t) return(c(est=NA,lo=NA,hi=NA,risk=sum(q$time>=t)))
 z<-cmprsk::cuminc(q$time,q$status,cencode=0)[['1 1']]
 if(is.null(z)) return(c(est=NA,lo=NA,hi=NA,risk=sum(q$time>=t)))
 j<-max(which(z$time<=t)); e<-z$est[j]; se<-sqrt(z$var[j])
 c(est=e,lo=max(0,e-1.96*se),hi=min(1,e+1.96*se),risk=sum(q$time>=t))
}
for(w in names(windows)) {
 wr<-windows[[w]]; a<-d[d$year>=wr[1]&d$year<=wr[2],]; s<-a[a$eligible,]
 for(v in defs) {
  q<-a[a[[v]],]; z<-s[s[[v]],]
  comp[[length(comp)+1]]<-data.frame(window=w,definition=v,n=nrow(q),age_mean=mean(q$age),year_mean=mean(q$year),female_pct=100*mean(q$sex=='Female'),distant_pct=100*mean(q$stage=='Distant'),unknown_stage_pct=100*mean(q$stage=='Unknown'),surgery_pct=100*mean(q$surgery=='Yes'),micro_pct=100*mean(q$micro))
  flows[[length(flows)+1]]<-data.frame(window=w,definition=v,raw_n=nrow(q),survival_n=nrow(z),cancer_deaths=sum(z$status==1),other_deaths=sum(z$status==2))
  for(t in c(36,60)) {cv<-cif_at(z,t); cifrows[[length(cifrows)+1]]<-data.frame(window=w,definition=v,months=t,n=nrow(z),cif=cv[1],lower95=cv[2],upper95=cv[3],at_risk=cv[4],support_warning=cv[4]<20)}
 }
 # Paired patient bootstrap respects the nesting of the typed subset.
 q<-s[s$current,]; ids<-unique(q$id); chunks<-split(seq_len(nrow(q)),q$id)
 for(t in c(36,60)) {
  base<-cif_at(q,t); typed<-cif_at(q[q$typed3,],t)
  if(base[4]<20||typed[4]<20) next
  bs<-replicate(bootstrap_n(300L),{ix<-unlist(chunks[sample(ids,length(ids),replace=TRUE)],use.names=FALSE); b<-q[ix,]; cif_at(b[b$typed3,],t)[1]-cif_at(b,t)[1]})
  ci<-quantile(bs,c(.025,.975),na.rm=TRUE)
  diffs[[length(diffs)+1]]<-data.frame(window=w,months=t,contrast='typed3 minus current (overlapping populations)',difference=typed[1]-base[1],lower95=ci[1],upper95=ci[2],bootstrap_replicates=bootstrap_n(300L))
 }
 # Direct typed-versus-NOS comparison in a common calendar window.
 q<-s[s$typed3|s$nos,]; q$typed<-as.integer(q$typed3)
 forms<-list(unadjusted='typed',calendar='typed + splines::ns(year,3)',demographics='typed + splines::ns(year,3) + splines::ns(age,3) + sex + race',stage='typed + splines::ns(year,3) + splines::ns(age,3) + sex + race + stage')
 if(w!='full') forms$treatment <- paste(forms$stage,'+ surgery')
 for(nm in names(forms)) {
  fit<-coxph(as.formula(paste('Surv(time,status==1) ~',forms[[nm]],'+ cluster(id)')),data=q,x=TRUE)
  sm<-summary(fit); row<-sm$conf.int['typed',]; ph<-tryCatch(cox.zph(fit)$table['typed','p'],error=function(e)NA)
  models[[length(models)+1]]<-data.frame(window=w,model=nm,n=nrow(q),events=sum(q$status==1),HR=row[1],lower95=row[3],upper95=row[4],p=sm$coefficients['typed','Pr(>|z|)'],PH_p=ph)
 }
}
writeout(do.call(rbind,comp),'composition'); writeout(do.call(rbind,cifrows),'cumulative_incidence'); writeout(do.call(rbind,flows),'survival_flow'); writeout(do.call(rbind,models),'sequential_cox'); writeout(do.call(rbind,diffs),'paired_bootstrap_risk_difference')
saveRDS(d,file.path(out,'broad_derived_with_ids.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
cat('AUDIT_COMPLETE\n'); print(do.call(rbind,diffs)); print(do.call(rbind,models))
