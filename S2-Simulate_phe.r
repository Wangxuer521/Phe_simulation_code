###Phenotypic simulation parameters: h2 = 0.8，QTL_number = 100。

library("data.table")

#Read genotype file
geno = as.data.frame(fread("./S1-Select_SNPs/Output_file/selected_SNPs.raw"))
cleang<-geno[,-c(1,3,4,5,6)]
rownames(cleang) = cleang[,1]
num_SNPs = ncol(geno)-6
colname<-c("ID",paste0("SNP", c(1:num_SNPs)))
colnames(cleang)<- colname
cleang1<-cleang[,-1]

#Randomly select 100 QTLs from the genotypes
data<-c(1:num_SNPs)
locinum<-sample(data,100,replace=FALSE)
qtlpos<-locinum[order(locinum)]

#QTL genotypes matrix
QTLloci<-cleang1[,qtlpos]

#Calculate the frequency of the given alleles for each QTL
p<-apply(QTLloci, 2, mean)/2

#calculate the variance of each QTL genotype
mult <- 2*p*(1-p)
qtl<-as.matrix(QTLloci)

#Centralize each column
qtlcen<-sweep(qtl, 2, colMeans(qtl))

#Simulation of allele substitution effects of QTL
eff<-sqrt((0.8*(mult^(-1)))/(length(mult)))

#Half of the QTL allele substitution effects were set as negative and half as positive, so that the mean breeding value of the population is zero
ef<-cbind(1:length(eff),eff)
rows<-sample(1:length(eff),length(eff)/2,replace=FALSE)
minus<-ef[rows,]
minus[,2]<-minus[,2]*-1
effe<-rbind(ef[-rows,],minus)
effec<-effe[order(effe[,1]),2]

#Simulated the real breeding values of individuals
tbv<-qtlcen%*%effec
snpvar<-mult*(effec^2)
num_IDs = nrow(geno)
res<-rnorm(num_IDs,mean=0,sd=sqrt(0.2))

#Simulate residual effects
re<-sqrt((0.2*var(tbv))/(0.8*var(res)))[1]*(res-rep(mean(res),num_IDs))

#Simulate phenotype
phe<-tbv+re
indinfo<-cbind(cleang[,1],tbv,re,phe)
qtlinfo<-cbind(qtlpos,p,mult,effec,snpvar)
write.table(qtlinfo,"qtlinfo",row.names=FALSE,col.names=FALSE,quote=F)   
#5 columns：SNP, allele frequency (p), expected variance (2pq), allele substitution effect, SNP variation
write.table(indinfo,"indinfo",row.names=FALSE,col.names=FALSE,quote=F)
#4 columns：ID, TBV, residual, phenotype

