# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
options(stringsAsFactors=FALSE)
library(survival)
d<-readRDS('descriptive/broad_derived_with_ids.rds')
z<-d[d$year>=2010 & d$year<=2022,]
q<-z[z$eligible & (z$typed3|z$nos),]
out<-'original_overlap'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
a<-do.call(rbind,lapply(c(TRUE,FALSE),function(g){x<-q[q$typed3==g,];data.frame(group=ifelse(g,'Detailed','NOS'),n=nrow(x),person_months=sum(x$time),observed_median=median(x$time),reverse_KM_median=unname(summary(survfit(Surv(time,status==0)~1,data=x))$table['median']),unknown_race=sum(grepl('Unknown',x[[8]])),unknown_stage=sum(x$stage=='Unknown'),unknown_surgery=sum(x$surgery=='Unknown'),missing_age=sum(is.na(x$age)),missing_year=sum(is.na(x$year)),missing_sex=sum(is.na(x$sex)))}))
write.csv(a,file.path(out,'followup_missing_audit.csv'),row.names=FALSE)
write.csv(data.frame(total=nrow(z),invalid_time=sum(is.na(z$time)),unknown_status=sum(is.na(z$status)),both=sum(is.na(z$time)&is.na(z$status)),excluded=sum(!z$eligible)),file.path(out,'flow_exclusion_audit.csv'),row.names=FALSE)
print(a)
