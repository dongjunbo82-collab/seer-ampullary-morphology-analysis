# Rscript figure2.R SEER8.txt SEER17.txt SEER17_RESULT_DIR FIGURE_OUTPUT_DIR
args<-commandArgs(TRUE);stopifnot(length(args)==4)
A<-args[3];P<-args[4];dir.create(P,recursive=TRUE,showWarnings=FALSE)
annual<-read.csv(file.path(A,'18_annual_surgery_coverage.csv'))
codes<-c(8140,8144,8163,8210,8211,8255,8260:8263,8480,8481,8490)
micro<-c('Positive histology','Positive exfoliative cytology, no positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies','Positive microscopic confirm, method not specified')
flow<-do.call(rbind,lapply(1:2,function(i){
 z<-read.delim(args[i],colClasses='character',check.names=FALSE,quote='"',na.strings=c('','NA'),fileEncoding='UTF-8')
 h<-as.integer(z[['Histologic Type ICD-O-3']]);ref<-h %in% codes & z[['Diagnostic Confirmation']] %in% micro
 typed<-ref & h %in% c(8144,8163,8255)
 data.frame(dataset=c('SEER 8, 1975-2023','SEER 17, 2018-2022')[i],source=nrow(z),reference=sum(ref),excluded=sum(!ref),detailed=sum(typed),other_reference=sum(ref & !typed))
}))
stopifnot(identical(flow$source,c(5695L,3112L)),identical(flow$reference,c(4811L,2675L)))
blue<-'#246B87';gray<-'#777777';black<-'#222222'
if(.Platform$OS.type=='windows')windowsFonts(Arial=windowsFont('Arial'))
box<-function(x,y,w,h,label,fill='white',col=black,cex=.86){rect(x-w/2,y-h/2,x+w/2,y+h/2,col=fill,border=col,lwd=.9);text(x,y,label,cex=cex,col=col)}
draw<-function(){
 layout(matrix(c(1,2,3),3,1),heights=c(1.13,1,.25))
 par(family='Arial',mgp=c(2.3,.6,0),tcl=-.25,lwd=.8)
 par(mar=c(.5,.8,2.6,.8));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
 mtext('a  Source and reference populations',3,adj=0,line=.8,font=2,cex=1)
 for(j in 1:2){
  x<-c(.245,.755)[j];f<-flow[j,]
  text(x,.96,f$dataset,font=2,cex=.94)
  box(x,.81,.43,.17,sprintf('Malignant ampullary tumors\nAge 20 years or older\nn = %s',format(f$source,big.mark=',')))
  arrows(x,.72,x,.59,length=.07)
  text(x+.02,.65,sprintf('Excluded: %s',format(f$excluded,big.mark=',')),adj=0,cex=.77)
  box(x,.48,.43,.20,sprintf('Reference adenocarcinoma\nMicroscopic confirmation\nn = %s',format(f$reference,big.mark=',')))
  segments(x,.38,x,.31);segments(x-.115,.31,x+.115,.31)
  arrows(c(x-.115,x+.115),.31,c(x-.115,x+.115),.26,length=.06)
  box(x-.115,.16,.205,.19,sprintf('Detailed codes\nn = %s\n%.1f%%',format(f$detailed,big.mark=','),100*f$detailed/f$reference),fill='#EAF2F5',col=blue,cex=.79)
  box(x+.115,.16,.205,.19,sprintf('Other reference\nn = %s\n%.1f%%',format(f$other_reference,big.mark=','),100*f$other_reference/f$reference),fill='#F3F3F3',cex=.79)
 }
 text(.5,.03,'Counts are tumor records; SEER 8 and SEER 17 areas overlap.',cex=.76)
 par(mar=c(2,4.4,3.2,1))
 plot(2018:2022,rep(0,5),type='n',xlim=c(2017.8,2022.2),ylim=c(0,48),xaxt='n',xlab='',ylab='Observed detailed-code coverage (%)',bty='l',cex.lab=.92)
 axis(1,2018:2022,2018:2022)
 for(g in c('All','Yes','No')){
  z<-annual[annual$surgery==g,];col<-switch(g,Yes=blue,No=gray,All=black);pch<-switch(g,Yes=16,No=17,All=1);lty<-if(g=='All')2 else 1
  if(g!='All'){arrows(z$year,z$lower,z$year,z$upper,angle=90,code=3,length=.035,col=adjustcolor(col,.65),lwd=.8)}
  lines(z$year,z$pct,col=col,lwd=1.6,lty=lty);points(z$year,z$pct,pch=pch,col=col,cex=.9)
 }
 mtext('b  SEER 17 annual coverage by recorded surgery',3,adj=0,line=1.5,font=2,cex=1)
 legend('topleft',legend=c('Surgery','No surgery','Full reference'),col=c(blue,gray,black),pch=c(16,17,1),lty=c(1,1,2),bty='n',horiz=TRUE,cex=.87,inset=c(0,-.01))
 par(mar=c(.2,4.4,0,1));plot.new();plot.window(xlim=c(2017.8,2022.2),ylim=c(0,1))
 text(2017.75,.94,'Detailed / N',adj=1,xpd=NA,cex=.69,font=2)
 for(i in 1:3){g<-c('Yes','No','Unknown')[i];z<-annual[annual$surgery==g,];y<-c(.7,.4,.1)[i]
  text(2017.75,y,c('Surgery','No surgery','Unknown')[i],adj=1,xpd=NA,cex=.73,col=c(blue,gray,black)[i])
  text(z$year,y,ifelse(z$N>0,paste(z$detailed,z$N,sep='/'),'-'),cex=.74,col=c(blue,gray,black)[i])
 }
}
tiff(file.path(P,'Fig2.tiff'),width=7.2,height=7.5,units='in',res=600,compression='lzw',type='cairo');draw();dev.off()
png(file.path(P,'Fig2_review.png'),width=1728,height=1800,res=240,type='cairo');draw();dev.off()
cat('Figure 2 and aggregate sources verified.\n')
