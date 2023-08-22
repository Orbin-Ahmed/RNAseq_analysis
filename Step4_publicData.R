library(tidyverse)
library(rhdf5)
library(edgeR)

archs4.mouse <- "archs4_gene_mouse_v2.1.2.h5"
h5ls(archs4.mouse)

all.samples.mouse <- h5read(archs4.mouse, name="meta/samples/geo_accession")
dim(all.samples.mouse)

mySamples <- c("GSM2310941",
               "GSM2310942",
               "GSM2310943",
               "GSM2310944",
               "GSM2310945",
               "GSM2310946",
               "GSM2310947",
               "GSM2310948",
               "GSM2310949",
               "GSM2310950",
               "GSM2310951",
               "GSM2310952")

my.sample.locations <- which(all.samples.mouse %in% mySamples) # first time you've seen the %in% operator.

genes <- h5read(archs4.mouse, "meta/genes/gene_symbol")

expression <- h5read(archs4.mouse, "data/expression",
                     index=list(my.sample.locations, 1:length(genes)))

expression <- t(expression)

rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
colSums(expression)
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)
colSums(archs4.cpm)

table(rowSums(archs4.dgelist$counts==0)==12)
keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

sample_source_name <- h5read(archs4.mouse, "meta/samples/source_name_ch1")

sample_title <- h5read(archs4.mouse, name="meta/samples/title")

sample_characteristics<- h5read(archs4.mouse, name="meta/samples/characteristics_ch1")

studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                      treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

genotype <- factor(studyDesign$genotype)
treatment <- factor(studyDesign$treatment)
sampleName <- studyDesign$Sample_title

pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
ls(pca.res)
summary(pca.res)
pca.res$rotation
pca.res$x
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2) +
  geom_point(size=4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(treatment) %>%
  add_column(genotype) %>%
  add_column(sampleName)


pca.pivot <- pivot_longer(pca.res.df,
                          cols = PC1:PC4,
                          names_to = "PC",
                          values_to = "loadings")

ggplot(pca.pivot) +
  aes(x=sampleName, y=loadings, fill=treatment) +
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()