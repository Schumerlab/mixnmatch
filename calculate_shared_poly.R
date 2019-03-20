args<-commandArgs(TRUE)
num_indiv<-as.numeric(args[1])
num_sites_par1<-as.numeric(args[2])
num_sites_par2<-as.numeric(args[3])

#determine SFS

x<-1:num_indiv
y<-1/x
total<-sum(y)
SFS<-cbind(x/num_indiv,1/(x*total))
write.table(cbind(x/num_indiv,1/(x*total)),file="SFS",sep="\t",row.names=FALSE,col.names=FALSE)

distribution_par1<-{}
distribution_par2<-{}

for (j in 1:length(SFS[,1])){
distribution_par1<-c(distribution_par1, rep(SFS[j,1],round(SFS[j,2]*num_sites_par1)))
}

for (k in 1:length(SFS[,1])){
distribution_par2<-c(distribution_par2, rep(SFS[k,1],round(SFS[k,2]*num_sites_par2)))
}

write.table(distribution_par1,file="polymorphism_distribution_par1",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(distribution_par2,file="polymorphism_distribution_par2",sep="\t",row.names=FALSE,col.names=FALSE)
