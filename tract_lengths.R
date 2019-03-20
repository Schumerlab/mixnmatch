args<-commandArgs(TRUE)
ma<-as.numeric(args[1]) #mix prop
ta<-as.numeric(args[2]) #admixture time
xa<-(as.numeric(args[3])) #recombination rate
chr<-(as.numeric(args[4])) #length

#Define the exponential distribution function
exponential<-function(m,t,x){
M<-rexp(1,rate=(m*t))
round(M*(1/x)*100*1000)
}

#Draw from the function until the chromosome length is generated

ancestry_tracts<-{}
total_length=0

while(total_length < chr){
par1_tract<-exponential(ma,ta,xa)
par2_tract<-exponential(1-ma,ta,xa)

total_length<-total_length+par1_tract+par2_tract

	if(total_length>chr){
		diff=total_length-chr
	if(par2_tract>diff){
		par2_tract=par2_tract-diff;
	} else{
	if(par2_tract<diff){
		diff=diff-par2_tract
		#return(diff)
		par2_tract=0
		par1_tract=par1_tract-diff
		}#take off both
	     }#only if required to remove some of the previous tract
	}#correct hanging tract
	
ancestry_tracts<-c(ancestry_tracts,par1_tract,par2_tract)

}

write.table(cbind(ancestry_tracts),row.names=FALSE,col.names=FALSE,sep="\t")
