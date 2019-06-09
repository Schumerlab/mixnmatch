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

last_interval=as.numeric(arrArgs[6])

data$V1<-round(data$V1*chromlength)
data$V2<-round(data$V2*chromlength)
data<-subset(data,data$V1!=data$V2) #0 bp intervals not allowed
cM_per_bp<-(((data$V2-data$V1)*data$V3)/sum((data$V2-data$V1)*data$V3)*Mlength)

data$V4<-cM_per_bp/(data$V2-data$V1)

current_length=as.numeric(arrArgs[7])
counter=0

bp=last_interval

x=0
y=0

track=0
rate_index_last=1
while(current_length<stop){

rate_index=which.min(c(abs(data$V1-bp),abs(data$V2-bp)))
total_heat=(data$V2[rate_index]-data$V1[rate_index])*data$V4[rate_index]

if(rate_index <= length(data$V2)){

if(counter==0){
if((current_length+total_heat)<start){
current_length=current_length+total_heat
bp=data$V2[rate_index]+1
} else{
current_length=current_length+data$V4[rate_index]
bp=bp+1
}#add whole or partial window
}#if still within start interval

if(counter==1){
if((current_length+total_heat)<stop){
current_length=current_length+total_heat
bp=data$V2[rate_index]+1
} else{
current_length=current_length+data$V4[rate_index]
bp=bp+1
}#add whole or partial window
}#if still within stop interval

if((current_length>=start) & (counter==0)){
x=last_interval #deal w/windows that ended in the middle from previous interval
counter=1
}

if(current_length>=stop){
y=bp
}

rate_index_last=rate_index

}else{current_length=stop}#don't try to run through the end of the file

}#keep churning through

if(x>0 & y>0){
write.table(cbind(x,y,0),sep="\t",col.names=FALSE,row.names=FALSE)
} else{
write.table(cbind(last_interval,chromlength,1),sep="\t",col.names=FALSE,row.names=FALSE)
}
