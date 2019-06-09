arrArgs<-commandArgs(trailingOnly = TRUE);

options(scipen=999)

infile<-as.character(arrArgs[1])
data<-read.csv(file=infile,sep="\t",head=FALSE)

start=as.numeric(arrArgs[2])
stop=as.numeric(arrArgs[3])

base_rate=as.numeric(arrArgs[4])
base_rate=(base_rate/1000)
chromlength=as.numeric(arrArgs[5])
Mlength=chromlength*base_rate

data$V1<-round(data$V1*chromlength)
data$V2<-round(data$V2*chromlength)
data<-subset(data,data$V1!=data$V2) #0 bp intervals not allowed
cM_per_bp<-(((data$V2-data$V1)*data$V3)/sum((data$V2-data$V1)*data$V3)*Mlength)

data$V4<-cM_per_bp/(data$V2-data$V1)

rate_index_start=which.min(c(abs(data$V1-start),abs(data$V2-start)))
rate_index_stop=which.min(c(abs(data$V1-stop),abs(data$V2-stop)))

if(rate_index_start == rate_index_stop){
total_rate=(stop-start)*data$V4[rate_index_start]
} else{
start_edge=(data[rate_index_start,]$V2-start)*data[rate_index_start,]$V4
stop_edge=(stop-data[rate_index_start,]$V1)*data[rate_index_stop,]$V4
intermediate<-subset(data,as.numeric(rownames(data)) >rate_index_start & as.numeric(rownames(data))<rate_index_stop)
intermediate_rate<-sum((intermediate$V2-intermediate$V1+1)*intermediate$V4)
total_rate=start_edge+stop_edge+intermediate_rate
}

write.table(cbind(total_rate),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")