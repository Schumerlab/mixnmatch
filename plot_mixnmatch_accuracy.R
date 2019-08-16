arrArgs<-commandArgs(trailingOnly = TRUE);

infile<-as.character(arrArgs[1])
data<-read.csv(file=infile,sep="\t",head=FALSE)

accuracy<-{}
errors_false_het<-{}
errors_false_homopar1<-{}
errors_false_homopar2<-{}
errors_by_tract<-{}
for(x in 1:length(unique(data$V1))){
focal<-subset(data,data$V1==paste("indiv",x,sep=""))
accuracy<-c(accuracy,sum(focal$V8)/(sum(focal$V8)+sum(focal$V9)))
errors_false_het<-c(errors_false_het,sum(subset(focal$V4,focal$V7!= "par1par2" & focal$V7!="par2par1")))
errors_false_homopar1<-c(errors_false_homopar1,sum(subset(focal$V6,focal$V7!= "par2par2")))
errors_false_homopar2<-c(errors_false_homopar2,sum(subset(focal$V5,focal$V7!= "par1par1")))
errors_by_tract<-rbind(errors_by_tract,cbind(focal$V3-focal$V2,focal$V9))
}#for all individuals

per_bp_error<-errors_by_tract[,2]/errors_by_tract[,1]
false_het_ratio<-errors_false_het/(errors_false_het+errors_false_homopar1+errors_false_homopar2)
false_homo1_ratio<-errors_false_homopar1/(errors_false_het+errors_false_homopar1+errors_false_homopar2)
false_homo2_ratio<-errors_false_homopar2/(errors_false_het+errors_false_homopar1+errors_false_homopar2)
allerrors<-(errors_false_het+errors_false_homopar1+errors_false_homopar2)

pdf(paste(infile,"_accuracy.pdf",sep=""),width=3.5,height=4)
hist(accuracy,xlab="Accuracy per individual",main="",col="lightblue",xlim=c(min(accuracy,na.rm=TRUE)-0.05,max(accuracy,na.rm=TRUE)+0.05),cex.lab=1.15,cex.axis=1.1)
dev.off()

pdf(paste(infile,"_errors.pdf",sep=""),width=8,height=4)
par(mfrow=c(1,4))
hist(allerrors,xlab="Number of errors\nper indiv",main="",col="lightblue",xlim=c(min(allerrors,na.rm=TRUE)-50,max(allerrors,na.rm=TRUE)+50),cex.lab=1.15,cex.axis=1.1)
hist(false_het_ratio,xlab="Proportion of errors\nfalse heterozygotes",main="",col="lightblue",xlim=c(min(false_het_ratio,na.rm=TRUE)-0.05,max(false_het_ratio,na.rm=TRUE)+0.05),cex.lab=1.15,cex.axis=1.1)
hist(false_homo1_ratio,xlab="Proportion of errors\nfalse parent 1",main="",col="lightblue",xlim=c(min(false_homo1_ratio,na.rm=TRUE)-0.05,max(false_homo1_ratio,na.rm=TRUE)+0.05),cex.lab=1.15,cex.axis=1.1)
hist(false_homo2_ratio,xlab="Proportion of errors\nfalse parent 2",main="",col="lightblue",xlim=c(min(false_homo2_ratio,na.rm=TRUE)-0.05,max(false_homo2_ratio,na.rm=TRUE)+0.05),cex.lab=1.15,cex.axis=1.1)
dev.off()

pdf(paste(infile,"_accuracy_vs_tract_length.pdf",sep=""),width=3.5,height=3.5)
plot(errors_by_tract[,1]/1000,per_bp_error,pch=20,xlab="Tract length (kb)",ylab="Error rate per basepair",col=rgb(120/255,132/255,150/255,alpha=0.1))
dev.off()
