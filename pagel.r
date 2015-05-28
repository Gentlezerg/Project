## for prepare the data (20 genes in test)
rm(list=ls())
lz =read.table("LZprofile.txt",header=TRUE,row.names=1)
clostreeid=read.csv("ClosTreeID.txt",sep=",",header=FALSE)
treegis=read.csv("TreeGIs.txt",sep=",",header=FALSE)
Bacnames=read.csv("NewTreeNamesID.txt",sep=",",header=FALSE)
# as.character()
lzcol = colnames(lz)

#write.table(mydata, "c:/mydata.txt", sep="\t")

#dim(clostreeid) 695 1
# treegis[,2]   719 2
lzprofile = lz[,colnames(lz) %in% Bacnames[,1]]
# dim(lz) 6505 860
# dim(lzprofile) 6505 662

i = 0
vector=c()
for (id in clostreeid[,1]) {
	if (!(id %in% treegis[,2])){
		print(id)
		i=i+1
		vector=c(vector,id)
	}

}
# m = matrix(sample(0:1,695*2, replace=TRUE),695,2)
# testprofile = cbind(clostreeid[,1],m)
# write.table(testprofile, "testClostree.txt", sep="\t",row.names=FALSE,col.names=FALSE)

df=data.frame(t(Bacnames[,2]))
colnames(df)=Bacnames[,1]
colnames(lzprofile) = df[1,colnames(lzprofile)]

lzcol = colnames(lzprofile)
i = 0
vector=c()
for (id in clostreeid[,1]) {
	if (!(id %in% lzcol)){
		print(id)
		i=i+1
		vector=c(vector,id)
	}

}

################## for one pair ##########
m = matrix(0,33,2)
missprofile = cbind(vector,m)
dfmissprofile=missprofile[,2:3]
rownames(dfmissprofile)=missprofile[,1]
# the bacs on the tree can't be found in profile
# gene loop 6000 * (6000-1)
gpair=c("496549815", "511538055")
# only 0 1
lzprofile[lzprofile>0]=1
pairprofile=lzprofile[gpair,]
pairprofile=t(pairprofile)

pairprofile = rbind(pairprofile,dfmissprofile)
write.table(pairprofile, "pairprofile.txt", sep="\t",row.names=TRUE,col.names=FALSE,quote=FALSE)
write.csv(pairprofile, "pairprofile.csv")

############################

################# for loop profile ################ test using 20 genes
m = matrix(0,33,20)
rownames(m) = vector
lzprofile[lzprofile>0]=1
loopprofile = lzprofile[1:20,]
loopprofile = t(loopprofile)
loopprofile = rbind(loopprofile,m)
write.table(loopprofile, "loopprofile.txt", sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.csv(loopprofile, "loopprofile.csv")

########

with open("pairprofile.txt.log.txt",'r') as f :
 for line in f:
     pass
line.split()[1]




