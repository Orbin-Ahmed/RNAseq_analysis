
library(tidyverse)
library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)

group <- factor(targets$treatment)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
fit <- lmFit(v.DEGList.filtered.norm, design)

contrast.matrix <- makeContrasts(infection = Vehicle - Dexamethasone,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)

ebFit <- eBayes(fits)

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")

myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)

vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

head(results)
summary(results)
vennDiagram(results, include="up")

head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

write_tsv(diffGenes.df,"DiffGenes.txt")

library(IsoformSwitchAnalyzeR)

targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = group) %>%
  dplyr::select(sampleID, condition)

Txi_trans <- importIsoformExpression(sampleVector = path)

colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels)
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels)

mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  ignoreAfterPeriod=TRUE,
  isoformExonAnnoation = "Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta       = "Homo_sapiens.GRCh38.cdna.all.fa",
  showProgress = TRUE)

mySwitchList <- isoformSwitchAnalysisCombined(
  switchAnalyzeRlist   = mySwitchList,
  pathToOutput = 'isoform_output')

extractSwitchSummary(mySwitchList)

extractTopSwitches(
  mySwitchList,
  filterForConsequences = TRUE,
  n = 50,
  sortByQvals = FALSE)
switchPlot(
  mySwitchList,
  gene='FCGR3B',
  condition1 = 'disease',
  condition2 = 'healthy',
  localTheme = theme_bw())