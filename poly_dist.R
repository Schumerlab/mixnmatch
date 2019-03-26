arrArgs<-commandArgs(trailingOnly = TRUE);

shared_poly<-as.numeric(arrArgs[1])
poly1<-as.numeric(arrArgs[2])
poly2<-as.numeric(arrArgs[3])
length<-as.numeric(arrArgs[4])
num_indiv<-as.numeric(arrArgs[5])
aims<-as.character(arrArgs[6])
error<-as.numeric(arrArgs[7])

aims_dat<-read.csv(file=aims,sep="\t",head=FALSE)

x<-1:num_indiv
y<-1/x
total<-sum(y)
SFS<-cbind(x/num_indiv,1/(x*total))

distribution_par1<-{}
distribution_par2<-{}

num_sites_par1<-round(poly1*length) - shared_poly
num_sites_par2<-round(poly2*length) - shared_poly

for (j in 1:length(SFS[,1])){
distribution_par1<-c(distribution_par1, rep(SFS[j,1],round(SFS[j,2]*num_sites_par1)))
}

for (k in 1:length(SFS[,1])){
distribution_par2<-c(distribution_par2, rep(SFS[k,1],round(SFS[k,2]*num_sites_par2)))
}

write.table(sample(distribution_par1),file="polymorphism_distribution_par1",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(sample(distribution_par2),file="polymorphism_distribution_par2",sep="\t",row.names=FALSE,col.names=FALSE)

#determine frequency of aims and shared polymorphism
freq<-{}
for (m in 1:length(aims_dat[,1])){

shared<-rbinom(1,1,shared_poly)

if(shared == 1){
freq<-rbind(freq,cbind(0.5,0.5))
} else{
curr_error<-runif(1,0,error)
freq<-rbind(freq,cbind(1-curr_error,curr_error))
}#site is shared or not

}#all sites

write.table(freq,file="shared_polymorphism_distribution",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
