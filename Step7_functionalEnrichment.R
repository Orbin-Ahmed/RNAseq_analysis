library(tidyverse)
library(limma)
library(gplots)
library(DT)
library(GSEABase)
library(Biobase)
library(GSVA)
library(gprofiler2)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
gost.res <- gost(rownames(myTopHits), organism = "mmusculus", correction_method = "fdr")
mygostplot = gostplot(gost.res, interactive = F, capped = T)

publish_gostplot(
  mygostplot,
  highlight_terms = c("GO:0034987"),
  filename = NULL,
  width = NA,
  height = NA)

publish_gosttable(
  gost.res,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = NULL,
  ggplot=TRUE)

msigdbr_species()
hs_gsea <- msigdbr(species = "Mus musculus")
hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)


hs_gsea_c2 <- msigdbr(species = "Mus musculus",
                      category = "C2") %>%
  dplyr::select(gs_name, gene_symbol)


mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)

mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)

gseaplot2(myGSEA.res, 
          geneSetID = c(6, 47, 262),
          pvalue_table = FALSE,
          title = myGSEA.res$Description[47])

myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "Dexamethasone",
    NES < 0 ~ "vehicle"))

ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

mySig <- toupper(rownames(myTopHits))
mySig2 <- sample((rownames(v.DEGList.filtered.norm$E)), size = 50, replace = FALSE)
collection <- list(real = mySig, fake = mySig2)
camera.res <- camera(v.DEGList.filtered.norm$E, collection, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

mroast(v.DEGList.filtered.norm$E, collection, design, contrast=1)

broadSet.C2.ALL <- getGmt("m2.all.v2023.1.Mm.symbols.gmt", geneIdType=SymbolIdentifier())

broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)
camera.res <- camera(v.DEGList.filtered.norm$E, broadSet.C2.ALL, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# filter based on FDR and display as interactive table
camera.df <- filter(camera.df, FDR<=0.01)

datatable(camera.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)


camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "Dexamethasone",
    Direction == "Down" ~ "vehicle"))


camera.df.sub <- camera.df %>%
  dplyr::filter(str_detect(setName, "CDC25|CYTOTOX"))


ggplot(camera.df[1:25,], aes(x = phenotype, y = setName)) +
  geom_point(aes(size = NGenes, color = Direction, alpha = -log10(FDR))) +
  scale_color_manual(values = c("up" = "green", "down" = "red")) +  # Adjust colors
  theme_bw()


v.DEGList.filtered.norm$E
GSVA.res.C2CP <- gsva(v.DEGList.filtered.norm$E,
                      broadSet.C2.ALL,
                      min.sz=5, max.sz=500,
                      mx.diff=FALSE,
                      method="gsva")

fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)
summary(res.C2CP)

diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
head(diffSets.C2CP)
dim(diffSets.C2CP)
 
hr.C2CP <- hclust(as.dist(1-cor(t(diffSets.C2CP), method="pearson")), method="complete")
hc.C2CP <- hclust(as.dist(1-cor(diffSets.C2CP, method="spearman")), method="complete")
mycl.C2CP <- cutree(hr.C2CP, k=2)
mycolhc.C2CP <- rainbow(length(unique(mycl.C2CP)), start=0.1, end=0.9) 
mycolhc.C2CP <- mycolhc.C2CP[as.vector(mycl.C2CP)] 


myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)

heatmap.2(diffSets.C2CP, 
          Rowv=as.dendrogram(hr.C2CP), 
          Colv=as.dendrogram(hc.C2CP), 
          col=myheatcol, scale="row",
          density.info="none", trace="none", 
          cexRow=0.9, cexCol=1, margins=c(10,14))