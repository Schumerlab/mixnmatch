arrArgs <- commandArgs(trailingOnly = TRUE);

if(length(arrArgs)<3){
stop("usage is: Rscript Determine_accuracy.R true_ancestry.bed individ genotypes_file\n") 
}

binsfile<-as.character(arrArgs[1])
indiv<-as.character(arrArgs[2])
infile<-as.character(arrArgs[3])

bins<-read.csv(file=binsfile,sep="\t",head=FALSE)

bins<-bins[order(bins$V6),]
#head(bins)

command0=paste("head -n 1 ",infile," > accuracy_",indiv,"_","genotypes_file",sep="")

system(command0)

command1=paste("grep ",indiv,"_ ",infile," ",">> accuracy_",indiv,"_","genotypes_file",sep="")
#command1
system(command1)

command2=paste("perl transpose_genotypes_tsv.pl ","accuracy_",indiv,"_","genotypes_file",sep="")
system(command2)
#command2
data<-read.csv(file=paste("accuracy_",indiv,"_","genotypes_file_transposed",sep=""),sep="\t",head=TRUE)

options(scipen=999)
whole_genome<-{}
start=1
stop=bins$V9[1]
counts_het=0
counts_par1=0
counts_par2=0
accurate_counts=0
inaccurate_counts=0
for (x in 1:length(bins$V1)){

focal<-subset(data,data$pos>=start & data$pos<=stop)
counts_het=length(subset(focal[,1],focal[,3] ==1))
counts_par1=length(subset(focal[,1],focal[,3] ==0))
counts_par2=length(subset(focal[,1],focal[,3] ==2))

if((bins$V4[x] == "par1" & bins$V8[x] == "par2") | (bins$V4[x] == "par2" & bins$V8[x] == "par1")){
accurate_counts=accurate_counts+counts_het
inaccurate_counts=inaccurate_counts+counts_par1+counts_par2
}
if(bins$V4[x] == "par1" & bins$V8[x] == "par1"){
accurate_counts=accurate_counts+counts_par1
inaccurate_counts=inaccurate_counts+counts_het+counts_par2
}
if(bins$V4[x] == "par2" & bins$V8[x] == "par2"){
accurate_counts=accurate_counts+counts_par2
inaccurate_counts=inaccurate_counts+counts_het+counts_par1
}

whole_genome<-rbind(whole_genome,cbind(indiv,start,stop,counts_het,counts_par1,counts_par2,paste(bins$V4[x],bins$V8[x],sep=""),accurate_counts,inaccurate_counts))

counts_het=0
counts_par1=0
counts_par2=0

start=stop+1
stop=start+bins$V9[x+1]
accurate_counts=0
inaccurate_counts=0
}

write.table(whole_genome,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
