# Rscript figure1.R SEER17_RESULT_DIR CONFIRMATION_RESULT_DIR FIGURE_OUTPUT_DIR
args<-commandArgs(TRUE);stopifnot(length(args)==3)
A<-args[1];C<-args[2];P<-args[3];dir.create(P,recursive=TRUE,showWarnings=FALSE)
r<-read.csv(file.path(A,'08_adjusted_coding_probabilities.csv'))
extra<-read.csv(file.path(C,'20_confirmation_sensitivity.csv'))
r<-rbind(r,extra)
sc<-c('broad','narrow','omit_stage','unknown_stage_category','unrestricted_broad','unrestricted_narrow')
z<-r[r$estimand=='Yes minus No',];z<-z[match(sc,z$scenario),]
comp<-read.csv(file.path(A,'06_composition_all_reference.csv'),check.names=FALSE)
co<-comp[match(c('surgery_pct','localized_pct','regional_pct','distant_pct','unknown_stage_pct'),comp$characteristic),c('Detailed','NOS','Reference')]
cov<-read.csv(file.path(A,'04_observed_coverage.csv'));cv<-cov[cov$variable=='surgery' & cov$level %in% c('Yes','No'),];cv<-cv[match(c('Yes','No'),cv$level),]
blue<-'#246B87';gray<-'#777777';black<-'#222222'
if(.Platform$OS.type=='windows')windowsFonts(Arial=windowsFont('Arial'))
draw<-function(){
 layout(matrix(1:3,3,1),heights=c(.8,1.18,1.05));par(family='Arial',cex=.96,mgp=c(2,.55,0),tcl=-.22,lwd=.8)
 par(mar=c(3.5,4.5,2.5,1.1))
 bx<-barplot(cv$broad_pct,names.arg=c('Recorded surgery','No recorded surgery'),col=c(blue,gray),border=NA,ylim=c(0,39),ylab='Coverage (%)',cex.names=.88)
 text(bx,cv$broad_pct+2.2,sprintf('%s / %s',format(cv$broad,big.mark=','),format(cv$N,big.mark=',')),cex=.86)
 mtext('a  SEER 17 observed coding coverage in 2018-2022',3,adj=0,line=.8,font=2,cex=.96)
 par(mar=c(4,10.2,2.5,6.6));yy<-6:1
 plot(100*z$estimate,yy,xlim=c(0,36),ylim=c(.5,6.8),yaxt='n',ylab='',xlab='Adjusted surgery-minus-no-surgery difference (pp)',pch=16,col=blue,bty='l',cex.lab=.85,cex.axis=.85)
 segments(100*z$lower95,yy,100*z$upper95,yy,col=blue,lwd=1.5)
 axis(2,yy,c('Broad: primary','8144/8163 only','Stage omitted','Unknown stage retained','No confirmation restriction','No confirmation; 8144/8163'),las=1,tick=FALSE,cex.axis=.73)
 text(37,yy,sprintf('%.1f (%.1f-%.1f)',100*z$estimate,100*z$lower95,100*z$upper95),adj=0,xpd=NA,cex=.73)
 text(37,6.65,'Difference (95% CI)',adj=0,xpd=NA,cex=.69,font=2)
 mtext('b  SEER 17 standardized coding differences',3,adj=0,line=.8,font=2,cex=.96)
 par(mar=c(4,8.5,3.6,8.2));yy<-5:1
 plot(co$Detailed,yy+.12,xlim=c(0,100),ylim=c(.45,5.6),yaxt='n',ylab='',xlab='Observed proportion (%)',pch=16,col=blue,bty='l',cex.axis=.87,cex.lab=.87)
 points(co$NOS,yy,pch=17,col=gray);points(co$Reference,yy-.12,pch=5,col=black)
 axis(2,yy,c('Recorded surgery','Localized stage','Regional stage','Distant stage','Unknown stage'),las=1,tick=FALSE,cex.axis=.82)
 text(103,yy,sprintf('%.1f / %.1f / %.1f',co$Detailed,co$NOS,co$Reference),adj=0,xpd=NA,cex=.71)
 text(103,5.65,'Detailed / NOS / Ref. (%)',adj=0,xpd=NA,cex=.68,font=2)
 legend('top',legend=c('Detailed (n = 579)','NOS (n = 1,951)','Reference (n = 2,675)'),pch=c(16,17,5),col=c(blue,gray,black),bty='n',horiz=TRUE,cex=.69,inset=c(-.025,-.02),xpd=NA)
 mtext('c  SEER 17 clinical composition',3,adj=0,line=1.9,font=2,cex=.96)
 mtext('The full reference cohort includes detailed-code and NOS records.',1,line=2.9,cex=.72)
}
tiff(file.path(P,'Fig1.tiff'),width=7.2,height=8.8,units='in',res=600,compression='lzw',type='cairo');draw();dev.off()
png(file.path(P,'Fig1_review.png'),width=1440,height=1760,res=200,type='cairo');draw();dev.off()
