# Descriptive calculations only; no manuscript or figures are written.
options(stringsAsFactors=FALSE)
library(survival)
d<-readRDS('descriptive/broad_derived_with_ids.rds')
out<-'original_overlap'
q<-d[d$eligible&d$year>=2010&(d$typed3|d$nos),]
q$group<-ifelse(q$typed3,'Detailed codes','NOS')
tab<-do.call(rbind,lapply(split(q,q$group),function(z)data.frame(group=unique(z$group),n=nrow(z),patients=length(unique(z$id)),age_mean=mean(z$age),age_sd=sd(z$age),year_mean=mean(z$year),female=sum(z$sex=='Female'),localized=sum(z$stage=='Localized'),regional=sum(z$stage=='Regional'),distant=sum(z$stage=='Distant'),unknown_stage=sum(z$stage=='Unknown'),surgery=sum(z$surgery=='Yes'),no_surgery=sum(z$surgery=='No'),unknown_surgery=sum(z$surgery=='Unknown'))))
write.csv(tab,file.path(out,'modern_characteristics.csv'),row.names=FALSE)
era<-cut(d$year,c(1974,1987,1997,2009,2017,2022,2023),labels=c('1975-1987','1988-1997','1998-2009','2010-2017','2018-2022','2023'))
er<-do.call(rbind,lapply(levels(era),function(e){z<-d[era==e,];data.frame(era=e,all=nrow(z),current=sum(z$current),NOS=sum(z$nos),typed=sum(z$typed3))}))
write.csv(er,file.path(out,'era_coverage.csv'),row.names=FALSE)
flow<-data.frame(step=c('All malignant','Age/year/site eligible before survival exclusions','Valid survival and event','Detailed or NOS in modern era'),n=c(nrow(d),sum(d$year>=2010&d$year<=2022),sum(d$eligible&d$year>=2010),nrow(q)))
write.csv(flow,file.path(out,'modern_flow.csv'),row.names=FALSE)
