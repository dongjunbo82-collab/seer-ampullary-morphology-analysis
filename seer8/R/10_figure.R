# Code-only distribution. Use run_analysis.R; outputs may contain sensitive records.
# Evidence contract: modern coding coverage is concentrated in recorded-surgery cases.
# Quantitative grid: observed denominators, adjusted coding probabilities, prognosis.
# TIFF-only output generated locally; no image is distributed with this code.
O<-'coding_selection'
cov<-read.csv(file.path(O,'modern_surgery_era_coverage.csv'))
pr<-read.csv(file.path(O,'adjusted_coding_probabilities.csv'))
rr<-read.csv('primary/risks.csv')
draw<-function(){
 layout(matrix(1:3,3,1),heights=c(1,1.25,1.1))
 par(family='sans',mgp=c(2.1,.6,0),tcl=-.25,cex=.9)
 par(mar=c(3.2,4.5,2,1))
 a<-cov[cov$level!='Unknown',];a<-a[match(c('2010-2017No','2010-2017Yes','2018-2022No','2018-2022Yes'),paste0(a$era,a$level)),]
 stopifnot(nrow(a)==4,!anyNA(a$coverage),identical(as.integer(a$detailed),c(3L,71L,13L,202L)))
 bx<-barplot(100*a$coverage,names.arg=c('No surgery','Surgery','No surgery','Surgery'),ylim=c(0,46),
 col=rep(c('#777777','#287E9C'),2),border=NA,ylab='Detailed-code coverage (%)',cex.names=.85)
 text(bx,100*a$coverage+2.5,sprintf('%d / %d',a$detailed,a$n),cex=.83)
 mtext('2010-2017',1,line=2,at=mean(bx[1:2]),cex=.85);mtext('2018-2022',1,line=2,at=mean(bx[3:4]),cex=.85)
 mtext('A  Modern coding coverage by surgery and era',3,adj=0,line=.6,font=2)
 b<-pr[pr$estimand=='probability',];yy<-c(9,8,6.5,5.5,4,3,2)
 par(mar=c(3,8.6,2,7.8))
 plot(100*b$estimate,yy,xlim=c(0,36),ylim=c(1.5,10),yaxt='n',ylab='',xlab='Adjusted probability of detailed coding (%)',pch=16,col='#287E9C',bty='l',cex.axis=.85)
 segments(100*b$lower95,yy,100*b$upper95,yy,col='#287E9C',lwd=1.4)
 axis(2,yy,c('No recorded surgery','Recorded surgery','2010-2017','2018-2022','Localized summary stage','Regional summary stage','Distant summary stage'),las=1,tick=FALSE,cex.axis=.8)
 text(37.5,yy,sprintf('%.1f (%.1f to %.1f)',100*b$estimate,100*b$lower95,100*b$upper95),adj=0,xpd=NA,cex=.77)
 text(37.5,9.9,'Probability (95% CI)',adj=0,xpd=NA,cex=.77,font=2)
 mtext('B  Standardized coding probabilities',3,adj=0,line=.6,font=2)
 nm<-c('supported3_all','supported3_surgery','narrow3_all','narrow3_surgery')
 z<-rr[rr$estimand=='RD',];z<-z[match(nm,z$scenario),]
 par(mar=c(5,8.6,2,7.8))
 plot(100*z$estimate,4:1,xlim=c(-23,15),ylim=c(.5,4.7),xlab='3-year cancer-death risk difference (percentage points)',ylab='',yaxt='n',pch=16,col='#287E9C',bty='l',cex.axis=.85)
 abline(v=0,lty=2,col='#999999');segments(100*z$lower95,4:1,100*z$upper95,4:1,col='#287E9C',lwd=1.4)
 axis(2,4:1,c('Detailed vs NOS, overall','Detailed vs NOS, surgery','Exclude 8255, overall','Exclude 8255, surgery'),las=1,tick=FALSE,cex.axis=.8)
 text(16.5,4:1,sprintf('%.1f (%.1f to %.1f)',100*z$estimate,100*z$lower95,100*z$upper95),adj=0,xpd=NA,cex=.77)
 text(16.5,4.7,'Difference (95% CI)',adj=0,xpd=NA,cex=.77,font=2)
 mtext('C  Prognostic comparisons in 2010-2019 diagnoses',3,adj=0,line=.7,font=2)
 mtext('Each risk comparison targets its own overlap population.',1,line=3.1,cex=.75)
 mtext('Within-surgery intervals remain wide; they do not establish equivalence.',1,line=4,cex=.75)
}
# Portable base-R sans font; no fonts or packages installed.
dir.create('figure',showWarnings=FALSE)
tiff('figure/Figure_1.tiff',width=183,height=235,units='mm',res=600,compression='lzw',pointsize=9);draw();dev.off()
