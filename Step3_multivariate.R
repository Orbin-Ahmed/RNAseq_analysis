library(tidyverse)
library(DT)
library(plotly)
library(gt)

targets
group <- targets$treatment
group <- factor(group)

log2.cpm.filtered.norm.df

distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "average")
plot(clusters, labels=sampleLabels)

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)

ls(pca.res)
summary(pca.res)
pca.res$rotation
pca.res$x

screeplot(pca.res)
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) +
  geom_point(size=4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)
  
pca.pivot <- pivot_longer(pca.res.df,
                          cols = PC1:PC4,
                          names_to = "PC",
                          values_to = "loadings")

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(Dexamethasone.AVG = (SCR1.1 + SCR1.2 + SCR2.1 + SCR2.2 + SCR3.1 + SCR3.2 + SCR4.1 + SCR4.2 + SCR5.1 + SCR5.1)/10,
         Vehicle.AVG = (SiP1.1 + SiP1.2 + SiP2.1 + SiP2.2 + SiP3.1 + SiP3.2 + SiP4.1 + SiP4.2 + SiP5.1 + SiP5.2)/10,
         LogFC = (Vehicle.AVG - Dexamethasone.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

mydata.df
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)
print(mydata.df)

mydata.filter <- mydata.df %>%
  dplyr::filter(geneID=="0610007P14Rik" | geneID=="0610009B22Rik" | geneID=="0610009D07Rik" | geneID=="0610009O20Rik" | geneID=="0610010F05Rik"
                | geneID=="0610011F06Rik" | geneID=="0610030E20Rik" | geneID=="0610031J06Rik" | geneID=="0610037L13Rik" ) %>%
  dplyr::select(geneID, Dexamethasone.AVG, Vehicle.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneID)) %>%
  dplyr::select(geneID, Dexamethasone.AVG, Vehicle.AVG, LogFC) %>%
  dplyr::arrange(desc(geneID))

gt(mydata.filter)
mydata.filter %>%
  gt() %>%
  fmt_number(columns=2:4, decimals = 1) %>%
  tab_header(title = md("**Regulators of skin pathogenesis**"),
             subtitle = md("*during cutaneous leishmaniasis*")) %>%
  tab_footnote(
    footnote = "Deletion or blockaid ameliorates Vehicle in mice",
    locations = cells_body(
      columns = geneID,
      rows = c(6, 7))) %>% 
  tab_footnote(
    footnote = "Associated with treatment failure in multiple studies",
    locations = cells_body(
      columns = geneID,
      rows = c(2:9))) %>%
  tab_footnote(
    footnote = "Implicated in parasite control",
    locations = cells_body(
      columns = geneID,
      rows = c(2))) %>%
  tab_source_note(
    source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))


datatable(mydata.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

myplot <- ggplot(mydata.df) + 
  aes(x=Dexamethasone.AVG, y=Vehicle.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("Vehicle vs. Dexamethasone") +
  theme_bw()

ggplotly(myplot)

myplot <- ggplot(mydata.df) +
  aes(x=Dexamethasone.AVG, y=Vehicle.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("Vehicle vs. Dexamethasone") +
  theme_bw()

ggplotly(myplot)