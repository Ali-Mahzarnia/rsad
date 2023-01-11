
library(readxl)
library(dplyr)

#path_master='/Users/ali/Desktop/Jul/apoe/MasterSheet_Experiments2021.xlsx'
#data2=read_xlsx(path_master, sheet = '18ABB11_readable02.22.22_BJ_Cor' )
path_master='/Users/ali/Desktop/nov/RSAD/human/AD_DECODE_data050922selected.csv'
data=read.csv(path_master)
datatemp=data%>%dplyr::select(genotype,Weight, sex, age, Subject,risk_for_ad)#subselect
#nchar(datatemp[111,1])
datatemp=na.omit(datatemp)
#datatemp[nchar(datatemp$Subject)==1,]=matrix(NA,1,dim(datatemp)[2])
#datatemp=na.omit(datatemp)
#datatemp[substr(datatemp$Subject,1,1)!="N",]=matrix(NA,1,dim(datatemp)[2])
#datatemp=na.omit(datatemp) ## ommit all na and zero character Subject and died durring
#datatemp$Subject=as.numeric(substr(datatemp$Subject,2,6)) # make Subject numeric
#datatemp=datatemp[datatemp$Genotype=="APOE22",]
# datatemp=datatemp[datatemp$Genotype %in% c("APOE22", "APOE33", "APOE44" , "APOE44HN"),]

# ####
# path_connec="/Users/ali/Desktop/nov/RSAD/apoe234_connectomes/"
# file_list=list.files(path_connec)
# temp_conn= read_xlsx( paste0(path_connec,file_list[1]) )
# temp_conn=temp_conn[,2: dim(temp_conn)[2]]
# connectivity=array( NA ,dim=c(dim(temp_conn)[1],dim(temp_conn)[2],dim(datatemp)[1]))
# dim(connectivity)
# 
# notfound=0
# ##read connec
# for (i in 1:dim(connectivity)[3]) {
#   
#   temp_index=which(datatemp$DWI[i]==as.numeric(substr(file_list,2,6)))
#   if (length(temp_index)>0) 
#   {
#   temp_connec=read_xlsx( paste0(path_connec,file_list[temp_index]) )
#   temp_connec=temp_connec[,2:dim(temp_connec)[2]]
#   colnames(temp_connec)=NA
#   temp_connec=as.matrix(temp_connec) ##
#   diag(temp_connec)=0### diag is 0
#   temp_connec=temp_connec/sum(temp_connec) ### divide by sum of connectviity
#   temp_connec=100*temp_connec
#   connectivity[,,i]=as.matrix(temp_connec)
#   }
#   else
#     notfound=c(notfound, datatemp$DWI[i])
#   
# }
# 
# notfound=notfound[2:length(notfound)]
# not_found_index=which( datatemp$DWI  %in%  notfound )
# 
# datatemp=datatemp[-not_found_index,]
# connectivity=connectivity[,,-not_found_index]
# sum(is.na(connectivity))

response=datatemp
#setwd(system("pwd", intern = T) )


RNA=read.delim('/Users/ali/Desktop/nov/RSAD/human/humanwithsymbol.txt',sep="\t", header = T)
RNA=t(RNA)
RNA_rows=rownames(RNA)

RNA_rows= as.numeric(gsub("\\D", "", RNA_rows))
#RNA_rows=RNA_rows[2:(length(RNA_rows)-1)]

RNA_data=RNA[2:(dim(RNA)[1]-1),]
colnames(RNA_data)=RNA[dim(RNA)[1],]
rownames(RNA_data)=RNA_rows[2:(length(RNA_rows)-1)]
#rownames(RNA_data)=NULL
sum(is.na(RNA_data))
#RNA_data[!is.na(RNA_data)]=NA
sum(is.na(RNA_data))


exist_response_ID=sub("\\s*\\:.*$", "",  response$Subject)
response$Subject=as.numeric(gsub("\\D", "", exist_response_ID )) 
#as.numeric(rownames(RNA_data))


found=0
notfound=0
for (jj in 1:dim(RNA_data)[1]) {
  
  index=which( as.numeric(row.names( RNA_data)   [jj])== as.numeric(response$Subject)      )
  if (length(index)>1) {cat(jj,",", index, ",")}
  if (length(index)>0){
        found=c(found,index)
  }
  
  else notfound=c(notfound,jj)
    
  }
found=found[2:length(found)]

response=response[found,]
dim(response)

#connectivity=connectivity[,,found]
#dim(connectivity)

#notfound=notfound[2:length(notfound)]
#RNA_data=RNA_data[-notfound,]

dim(RNA_data)

rownames(RNA_data)=NULL

RNA_data = as.data.frame(RNA_data)
response =as.data.frame(response)
response$RSAD1= RNA_data$RSAD1

cor(as.numeric(response$age),as.numeric(response$RSAD1) )
plot(as.numeric(response$age),as.numeric(response$RSAD1))

############################# Rsad1

response2 = response
#response2= response2 [ as.numeric(response2$RSAD1)>0,] 
#response2= response2 [ as.numeric(response2$RSAD1)<200,] 
response2$genotype[response2$genotype =="APOE23"] = "APOE22"
response2$genotype[response2$genotype =="APOE22"] = "APOE33"
response2$genotype[response2$genotype =="APOE34"] = "APOE44"
response2$risk_for_ad[response2$risk_for_ad ==3] = 2

# lm=lm( RSAD1~age*as.factor(genotype), data=response2)
# summary(aov(lm))
# response2$RSAD1=lm$resid


dodge <- position_dodge(width = 1)

plot1a<-ggplot(response2, aes(x=as.factor(sex), y=as.numeric(RSAD1), fill = as.factor(risk_for_ad))) +
  geom_violin(inherit.aes=TRUE,position=dodge) +
  #scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  #scale_fill_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  facet_grid(. ~ as.factor(risk_for_ad))  +
  facet_wrap(~as.factor(risk_for_ad)) +
  scale_alpha_discrete(range = c(0.4,0.8)) +
  geom_boxplot(color="black", outlier.color="black", width=0.2, alpha=.8, position=dodge) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2, alpha=0.6, position=dodge)+
  #geom_jitter(size = 0.1, height = 0, width = 0.1, aes(color = Sex)) +
  #labs(title = "VOL")+
  theme_minimal()+
  #background_grid(major = 'xy', minor = "none") + # add thin horizontal lines
  #panel_border() +
  theme_bw()+
  #labs(x = "genotype", y = "Diastolic_LV_Volume", title = "Diastolic_LV_Volume ") +
  #stat_summary(fun.y=median, geom="point", size=2, color="black") +
  theme(legend.position="bottom")+
  theme_bw()
plot1a
ggsave("human_Rsad1.png", dpi = 300)

summary(aov(RSAD1 ~ risk_for_ad *sex*genotype, data = response2))

save(response, file="response.rda")
save(connectivity, file="connectivity.rda")
save(RNA_data, file="RNA_data.rda")