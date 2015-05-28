rm(list=ls())
library(gplm)

####
ker="uniform"   # quartic triweight gaussian "uniform" 
k=3
###
# read the Patristic Distance Data
PatricDist=read.csv(file="clostRooted_distance.csv",header=TRUE,row.names=1,check.names=FALSE)
# TestDist=read.csv("testTree_Distance.csv",header=TRUE,row.names=1,check.names=FALSE)
lz =read.table("LZprofile.txt",header=TRUE,row.names=1)
Bacnames=read.csv("NewTreeNamesID.txt",sep=",",header=FALSE)
# colnames(Bacnames)=c('name','ID')
# as.character  %in%
# only matched ID distance
subPDM=PatricDist[as.character(Bacnames[,2])]
subPDM=subPDM[colnames(subPDM),]
#write.csv(subPDM,file="subPDM.csv")
# t = table(Bacnames$name)
# t = as.data.frame(t)
# t[t$Freq>1,]
# a = as.character(Bacnames$ID)
# a = as.data.frame(table(a))
# a[a$Freq>1,] # delete begin with 4
# droplevels
# which(table(Bacnames$ID)>1) find the reduplicated items ID has >1 Bac
#lzprofile2 = lz[,colnames(lz) %in% Bacnames[,1]]
lzprofile = lz[,colnames(lz) %in% Bacnames[,1]]

df=data.frame(t(Bacnames[,2]))
colnames(df)=Bacnames[,1]
colnames(lzprofile) = df[1,colnames(lzprofile)] # replace the bac name in lzprofile with feature id
#write.csv(lzprofile,file="lzprofile.csv")
# lzprofile = as.matrix(lzprofile)
lzprofile[lzprofile>0]=1
# bacid=Bacnames
# bacid$ID=as.character(bacid$ID) # convert ID to string
#write.csv(lzprofile,file="lzprofile01.csv")

#############
# cross check remove not exist
df1 = df[1,df[1,] %in% colnames(lzprofile)]
PDM = subPDM[,colnames(subPDM) %in% as.character(df1[1,])]
PDM = PDM[colnames(PDM),]

# remove lzprofile rows which are all 0
# don't need lzp = lzprofile x[!apply(x == 0, 1, all),] x[rowSums(x)>0,]

# reduce lzprofile
lzprofile = lzprofile[rowSums(lzprofile)>=200 & rowSums(lzprofile)<=400,]  ### too much

ngene = nrow(lzprofile)  # 540/6505
nbac = ncol(lzprofile)   # 662

# end of data preprocessing

## initial process
kneighbors = function(bid,distance.matrix,k)
{
     ordermatrix = order(distance.matrix[bid,])
     return(ordermatrix[2:(k+2)])
}
predmatrix = data.frame(matrix(NA, nrow=ngene, ncol=nbac))
rownames(predmatrix)=rownames(lzprofile)
colnames(predmatrix)=colnames(lzprofile)
resimatrix = predmatrix 

###

num = 0
for (bac in colnames(lzprofile)) {
	Idex=kneighbors(bac,PDM,k)
	for (gene in rownames(lzprofile)) {

lneigh = PDM[bac,Idex[1:k]]/PDM[bac,Idex[k+1]]
lneigh = ifelse(lneigh==1,0.999,lneigh)
vneigh = lzprofile[gene,colnames(PDM)[Idex[1:k]]]
wlneigh = kernel.function(as.numeric(lneigh), kernel=ker)
wlneigh = wlneigh/sum(wlneigh)
predmatrix[gene,bac] = ifelse(as.matrix(vneigh) %*% as.matrix(wlneigh)>0.5,1,0)
resimatrix[gene,bac] = lzprofile[gene,bac] - predmatrix[gene,bac]
}
num=num+1
print(num)
}
err = sum(abs(resimatrix))/(ngene*nbac)
print(err)
write.csv(resimatrix, paste0("ResidualMatrix",k,ker,".csv"))
write.table(t(c(k,ker,err)),"results.txt",col.names=F,row.names=F,append=T)

