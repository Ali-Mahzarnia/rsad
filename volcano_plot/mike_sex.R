#data=read.csv('/Users/ali/Desktop/aug/volcano_plot/mikes_sex.csv')%>%select(Gene.Symbol, log2.FC., ProbChiSq)
data=read.csv('/Users/ali/Desktop/aug/volcano_plot/mikes_sex.csv')%>%select(Gene.Symbol, log2.FC., fdr_p)

data[2,1]
len_na=length(data[which(data=="",arr.ind = T)])
data[which(data=="",arr.ind = T)]=matrix(NA,len_na, 1)
data[2,1]
data=na.omit(data)
data[2,1]
sum(data=="")
data$log2.FC.=as.numeric(data$log2.FC.)
#data$ProbChiSq=as.numeric(data$ProbChiSq)
data$fdr_p=as.numeric(data$fdr_p)
data=na.omit(data)

colnames(data)=c('gene_symbol', 'log2FoldChange','pvalue')



sum(abs(data$log2FoldChange)>10)

de=data[abs(data$log2FoldChange)<9,]
sum(abs(de$log2FoldChange)>10)





library(ggplot2)
# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=de, aes(x=log2FoldChange, y=pvalue)) + geom_point()
# Convert directly in the aes()
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
plot(p)
# Add more simple "theme"
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
plot(p)




# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
plot(p2)




# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.6 & de$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.6 & de$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2


## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
p3

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3


# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()


# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
final=ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size=8) +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  theme(axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),  
        axis.title.x = element_text(color = "black", size = 20),
        axis.title.y = element_text(color = "black", size = 20))

# theme(axis.text = element_text(size = 8, colour="black")) 

ggsave("mike_sex.png", plot = final, 
       device='png', scale=1, width=12, height=10, unit=c("in"), dpi=200)


