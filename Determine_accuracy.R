arrArgs <- commandArgs(trailingOnly = TRUE);

if(length(arrArgs)<4){
stop("usage is: Rscript Determine_accuracy.R true_ancestry.bed individ genotypes_file posterior1 posterior2 path\n") 
}

binsfile<-as.character(arrArgs[1])
indiv<-as.character(arrArgs[2])
infile<-as.character(arrArgs[3])
pp1<-as.character(arrArgs[4])
pp2<-as.character(arrArgs[5])
path<-as.character(arrArgs[6])

bins<-read.csv(file=binsfile,sep="\t",head=FALSE)


bins<-bins[order(bins$V6),]
#command 0 set

command0=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file",sep="")

system(command0)

command0b=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file_pp1",sep="")
command0c=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file_pp2",sep="")

system(command0b)
system(command0c)

#command 1 set
command1=paste("grep ",indiv,"_ ",infile," ",">> accuracy_",indiv,"_","genotypes_file",sep="")

system(command1)

command1b=paste("grep ",indiv,"_ ",pp1," ",">> accuracy_",indiv,"_","genotypes_file_pp1",sep="")
command1c=paste("grep ",indiv,"_ ",pp2," ",">> accuracy_",indiv,"_","genotypes_file_pp2",sep="")

system(command1b)
system(command1c)

command2=paste("perl ",path,"/transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file",sep="")
system(command2)

command2b=paste("perl ",path,"/transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file_pp1",sep="")
command2c=paste("perl ",path,"/transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file_pp2",sep="")

system(command2b)
system(command2c)

#command 2 set
data<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_transposed",sep=""),sep="\t",head=TRUE)
pp1<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_pp1_transposed",sep=""),sep="\t",head=TRUE)
pp2<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_pp2_transposed",sep=""),sep="\t",head=TRUE)

options(scipen=999)
whole_genome<-{}
start=1
stop=bins$V9[1]
counts_het=0
counts_par1=0
counts_par2=0
accurate_counts=0
inaccurate_counts=0
mean_pos<-{}
for (x in 1:length(bins$V1)){

focal<-subset(data,data$pos>=start & data$pos<=stop)
focalpp1<-subset(pp1,pp1$pos>=start & pp1$pos<=stop)
focalpp2<-subset(pp2,pp2$pos>=start & pp2$pos<=stop)

counts_het=length(subset(focal[,1],focal[,3] ==1))
counts_par1=length(subset(focal[,1],focal[,3] ==0))
counts_par2=length(subset(focal[,1],focal[,3] ==2))

if((bins$V4[x] == "par1" & bins$V8[x] == "par2") | (bins$V4[x] == "par2" & bins$V8[x] == "par1")){
accurate_counts=accurate_counts+counts_het
inaccurate_counts=inaccurate_counts+counts_par1+counts_par2
mean_pos=c(mean_pos,(focalpp1[,3]+focalpp2[,3]))
}
if(bins$V4[x] == "par1" & bins$V8[x] == "par1"){
accurate_counts=accurate_counts+counts_par1
inaccurate_counts=inaccurate_counts+counts_het+counts_par2
mean_pos=c(mean_pos,(1-focalpp1[,3]))
}
if(bins$V4[x] == "par2" & bins$V8[x] == "par2"){
accurate_counts=accurate_counts+counts_par2
inaccurate_counts=inaccurate_counts+counts_het+counts_par1
mean_pos=c(mean_pos,(1-focalpp2[,3]))
}
#write.table(cbind(mean_pos))

whole_genome<-rbind(whole_genome,cbind(indiv,start,stop,counts_het,counts_par1,counts_par2,paste(bins$V4[x],bins$V8[x],sep=""),accurate_counts,inaccurate_counts,mean(mean_pos)))

counts_het=0
counts_par1=0
counts_par2=0

start=stop+1
stop=start+bins$V9[x+1]
accurate_counts=0
inaccurate_counts=0
}

write.table(whole_genome,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
