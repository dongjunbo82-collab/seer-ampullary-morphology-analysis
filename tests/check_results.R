# Inspect locally generated full-mode outputs, never distribute those outputs.
args<-commandArgs(TRUE)
if(length(args)!=1)stop('Usage: Rscript tests/check_results.R PRIVATE_OUTPUT_DIRECTORY')
o<-normalizePath(args[1],mustWork=TRUE)
read<-function(p)read.csv(file.path(o,p),check.names=FALSE)
risks<-read('primary/risks.csv');draws<-read('primary/bootstrap.csv')
stopifnot(nrow(draws)==9000,length(unique(risks$scenario))==9,all(risks$B==1000))
for(s in unique(risks$scenario))for(e in c('control','detailed','RD')){
 r<-risks[risks$scenario==s&risks$estimand==e,]
 ci<-quantile(draws[draws$scenario==s,e],c(.025,.975),na.rm=TRUE)
 stopifnot(max(abs(ci-c(r$lower95,r$upper95)))<1e-12)
}
v<-read('primary/point_validation.csv')
stopifnot(all(v$independent_AJ_error<1e-10))
a<-read('coding_selection/audit_summary.csv')
stopifnot(a$value[a$check=='raw_mapping_mismatches']==0)
p<-read('coding_selection/adjusted_coding_probabilities.csv')
b<-read('coding_selection/selection_bootstrap.csv')
stopifnot(nrow(b)==1000,all(is.finite(p$estimate)),all(p$lower95<=p$upper95))
stopifnot(all(is.na(b$error)|b$error==''),all(is.na(b$warnings)|b$warnings==''))
verify<-read('coding_selection/independent_validation.csv')
stopifnot(all(verify$value[1:2]<1e-10))
cat('PASS: primary percentile intervals, independent AJ, mapping audit and coding-model verification.\n')
