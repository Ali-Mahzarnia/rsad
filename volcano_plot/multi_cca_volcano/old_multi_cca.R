
#run this right after running Step A matlab that would produce the connectivity data by reading and re-arranging
#install.packages("R.matlab")
library(R.matlab)
path="/Users/ali/Desktop/may/rnahuman/resultsconnectivity_all_ADDecode_Dipy.mat"
path2="/Users/ali/Desktop/may/rnahuman/resultsresponse_array.mat"


data=readMat(path)
connectivity=data$connectivity

# noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab

#temp=connectivity[-noreadcsf,-noreadcsf,1]
temp=connectivity[,,1]
indexlower=lower.tri(temp, diag=FALSE)
indexlowertrue=which(indexlower==TRUE)
temp=temp[indexlower]
len=sum(indexlower)  



data2=readMat(path2, fixNames=TRUE)
response=data2$response.array 
# 'MRI_Exam', 'sex', 'age', 'Weight', 'risk_for_ad', 'genotype'

#riskfactors=matrix(NA,  dim(response)[1], (dim(response)[2]-1))
riskfactors=matrix(NA,  (dim(response)[1]), (dim(response)[2]-1)) #
# 'sex', 'age', 'Weight', 'risk_for_ad', 'genotype'
#sum(riskfactors[,2]==3)


subjnameofconnectivity=data$subjlist



for (i in 1:dim(riskfactors)[1]) {
  ind=which(response[i,10]==subjnameofconnectivity)
  if (i!=ind) cat("here", i)
  #riskfactors[ind,]=response[ind, 2:(dim(response)[2]) ]
  temp=response[ind, 1:(dim(response)[2]) ];
  temp=temp[-c(10)];
  riskfactors[ind,]= temp;
}


###covar names as well:
pathnames0="/Users/ali/Desktop/may/sccapapr/resultsresponse_tablename.mat"
temp=readMat(pathnames0, fixNames=T)
temp=unlist(temp$varnames)
temp=temp[-c(9)]
temp=c("id"  ,temp)
colnames(riskfactors)=temp
#riskfactors=riskfactors[,-c(3)] # no weight
riskfactors=riskfactors[,c(1:11)] # phisio
riskfactorsorig=riskfactors
riskfactors=riskfactors[,c(2,3,10,11)] #only sex age gene and family
#riskfactors=riskfactors[,-c(3,6,8)] # no bmi if weight is there high corrs also no sys and pulse
#also no height between height and sex
rownames(riskfactors)=t(riskfactorsorig[,1])




######## pull all ad as mci
famindex=which(colnames(riskfactorsorig)=="risk_for_ad")
tempaaa=riskfactorsorig[,famindex]
tempaaa[tempaaa==3]=2
riskfactorsorig[,famindex]=tempaaa
#########


#riskfactorind=riskfactors>0
#sum(riskfactorind)
#riskfactors=riskfactors[riskfactorind,] # removing riskfactor 2,3

rna=read.delim('/Users/ali/Desktop/may/rnahuman/humanwithsymbol.txt',sep="\t", header = T)
sum(is.na(rna$sym))
j=1 # replace some nonfound rna sym with ids (282 of them)
for (i in 1:length(rna$sym)) {
  if (is.na(rna$sym[i])) { cat(j, '   ', i, '\n');rna$sym[i]=rna$gene_id[i]; j=j+1}
}
sum(is.na(rna$sym)) # check to see all replaced
rna$sym[25697]


colnames(rna) # we only have about 56 subjects rna seq
subjrnalen=length(which(colnames(rna) !='gene_id'  & colnames(rna) !='sym'  ))
len=dim(rna)[1]

rnadata=matrix(NA,  dim(riskfactors)[1], len) # -6 becasue of cfs removal

colomnames=colnames(rna);
colomnames=colomnames[1:(dim(rna)[2])]
colomnames=as.numeric(gsub(".*?([0-9]+).*", "\\1", colomnames))      

 for (i in 1:dim(riskfactors)[1]){
indexofsubj=which(as.numeric(rownames(riskfactors)[i]) == colomnames)
if (length(indexofsubj)>0){
  if (length(indexofsubj)>1){indexofsubj=indexofsubj[1] }
 rnadata[i,]=rna[,indexofsubj]
}
 }
dim(rnadata)
temprowremoved=which(is.na(rnadata),arr.ind=T)[,1]
temprowremoved=unique(temprowremoved)
rnadata=rnadata[-temprowremoved,]
dim(rnadata)
sum(is.na(rnadata))

riskfactors=riskfactors[-temprowremoved,]
riskfactorsorig=riskfactorsorig[-temprowremoved,]
#rnadata[,1]

#rnadata=rnadata[riskfactorind,]

#recordzerocols # these are zero cols that we remove and add at the edn 
# we rmove now because cca needs to standardize and sd of them are zero
inddrna=0
for (i in 1:dim(rnadata)[2]) if(sd(rnadata[,i])==0 ) {inddrna=rbind(inddrna,i);  cat ( i , sd(rnadata[,i]), "\n" );}
if (length(inddrna)>1){
inddrna=inddrna[2:dim(inddrna)[1]]
rnadata=rnadata[,-inddrna] }


inddz=0
for (i in 1:dim(riskfactors)[2]) if(sd(riskfactors[,i])==0 ) {inddz=rbind(indd,i);  cat ( i , sd(riskfactors[,i]), "\n" );}
if (length(inddz)>1){
inddz=inddz[2:dim(inddz)[1]]
riskfactors=riskfactors[,-inddz]
}



ageind=which(colnames(riskfactors)=="age")
medianage=median(riskfactors[,ageind])
agecat=riskfactors[,ageind];agecat[agecat<=medianage]=1;agecat[agecat>medianage]=2   #agecat
#riskfactors=riskfactors[,c(5,9)] # NO WEIGHT it is sex, age, diet, gene
#riskfactors[,2]=agecat;


# penalty=c(sqrt(dim(riskfactors)[2]),sqrt(dim(riskfactors)[2])*2.5,4*sqrt(dim(riskfactors)[2]) )
# out <- MultiCCA(xlist, type=c("standard", "standard", "standard"),
#                 penalty=penalty, ncomponents=1, ws=perm.out$ws.init)

#temp=connectivity[-noreadcsf,-noreadcsf,1]
temp=connectivity[,,1]
indexlower=lower.tri(temp, diag=FALSE)
indexlowertrue=which(indexlower==TRUE)
temp=temp[indexlower]
len=sum(indexlower)  

image=matrix(NA,  dim(connectivity)[3], len) # -6 becasue of cfs removal

for (i in 1:dim(connectivity)[3]){
  #temp=connectivity[-noreadcsf,-noreadcsf,i]
  temp=connectivity[,,i]
  indexlower=lower.tri(temp, diag=FALSE)
  temp=temp[indexlower]
  image[i,]=temp
}
dim(image)

image=image[-temprowremoved,]

indd=0
for (i in 1:dim(image)[2]) if(sd(image[,i])==0 ) {indd=rbind(indd,i);  cat ( i , sd(image[,i]), "\n" );}
if (length(indd)>1){
  indd=indd[2:dim(indd)[1]]
  image=image[,-indd] }









#lets run
## Not run:
#install.packages("PMA")
#install.packages("https://gitlab.oit.duke.edu/am983/PMA2/-/archive/master/PMA2-master.tar.gz", repos = NULL, type="source")
library(PMA2)
xlist= list (riskfactors, image, rnadata)
#perm.out <- MultiCCA.permute(xlist, nperm=1000)


#perm.out$penalties
###or load the data after this
#penalty=perm.out$bestpenalties
#penalty[2]=sqrt(dim(image)[2])*2.5
#penalty[2]=sqrt(dim(riskfactors)[2])*2.5

penalty=c(sqrt(dim(riskfactors)[2]),sqrt(dim(riskfactors)[2])*1.5,3*sqrt(dim(riskfactors)[2]) )
out <- MultiCCA(xlist, type=c("standard", "standard", "standard"),
                penalty=penalty, ncomponents=1, ws=perm.out$ws.init)
out$penalty

ws=out$ws
mean(ws[[3]]==0)
sum(ws[[3]]!=0)
print(out)

perm.out2 <- MultiCCA.permute(xlist, nperm=100, penalties = as.matrix(penalty) )

rnaw=matrix(NA, length(ws[[3]])+length(inddrna),1 )
#put those zeros back
rnaw[inddrna]=0
rnaw[-inddrna]=ws[[3]]


symnames=rna[,which(colnames(rna)=='sym')]
table=cbind(symnames[rnaw!=0],rnaw[rnaw!=0] )
colnames(table)=c('sym','multicca')
library(xlsx)
write.xlsx2(table, "60table.xlsx")




u=out$ws[[2]]

sum(u==0)
#len=length(u)
sum(u!=0)
u[u!=0]
sum(u==0)/len #sparsity 