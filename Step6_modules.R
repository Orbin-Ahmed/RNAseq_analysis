library(tidyverse)
library(limma)
library(RColorBrewer)
library(gplots)
library(gameofthrones)
library(heatmaply)
library(d3heatmap)

myheatcolors1 <- bluered(75)
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)

myheatcolors2 <- colorRampPalette(colors=c("blue","white","red"))(100)
myheatcolors3 <- brewer.pal(name="RdBu", n=11)
myheatcolors3 <- c("#fed976", "#268f9c")
got_palette <- got(75, option = "Arya")

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
dim(diffGenes)
diffGenes
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")

module.assign <- cutree(clustRows, k=2)

module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 

heatmaply(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module.color,
          #showticklabels=c(FALSE,FALSE),
          scale='row')

d3heatmap(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          row_side_colors = module.color,
          scale='row')

colnames(diffGenes) <- targets$treatment

diffGenes.AVG <- avearrays(diffGenes)

names(module.color) <- names(module.assign) 

module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df,
                          cols = 1:203,
                          names_to = "geneID",
                          values_to = "module")

module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))


ggplot(module.assign.pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

modulePick <- 2
myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA, 
          labRow = NA,
          col=rev(myheatcolors3), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20)) 

moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"module_upRegulated.tsv")

mySelectedGenes <- mydata.filter
mySelectedGenes.matrix <- as.matrix(mySelectedGenes)
mySelectedGenes
hr <- hclust(as.dist(1-cor(t(mySelectedGenes), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(mySelectedGenes, method="spearman")), method="average")

heatmap.2(mySelectedGenes.matrix, 
          Rowv=NA, Colv=NA, 
          col=myheatcolors1, 
          scale="row", density.info="none", 
          trace="none", labCol=NA, 
          cexRow=1.5, cexCol=1, margins=c(8,20), key = F)
