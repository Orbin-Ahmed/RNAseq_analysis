library(rhdf5)
library(tidyverse)
library(tximport)
library(ensembldb)
library(EnsDb.Mmusculus.v75)
library(biomaRt)
library(datapasta)

targets <- read_tsv("studydesign.txt")


path <- file.path(targets$sample, "abundance.tsv")

all(file.exists(path)) 

Tx <- transcripts(EnsDb.Mmusculus.v75, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

library(biomaRt) 

listMarts()
myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
available.datasets <- listDatasets(myMart)
dog.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "clfamiliaris_gene_ensembl")
dog.attributes <- listAttributes(dog.anno)

Tx.dog <- getBM(attributes=c('ensembl_transcript_id_version',
                             'external_gene_name'),
                mart = dog.anno)

Tx.dog <- as_tibble(Tx.dog)
Tx.dog <- dplyr::rename(Tx.dog, target_id = ensembl_transcript_id_version, 
                        gene_name = external_gene_name)


Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)