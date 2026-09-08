library(data.table)
library(dplyr)

# LOAD STUFF IN
pu1 <- fread(
  "HL60_SPI1_ENCFF797QGP_GRCh38_optimalIDR.bed.gz",
  header = FALSE
)

colnames(pu1) <- c(
  "chr",
  "start",
  "end",
  "name",
  "score",
  "strand",
  "signalValue",
  "pValue",
  "qValue",
  "peak"
)

hs_genes <- fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_relevant_genes.csv") %>% 
	.$hgnc_symbol %>% unique()

# GET PROMOTER REGIONS
Sys.setenv(TZ = "America/New_York")
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
library(AnnotationDbi) 

BiocManager::install(c(
  "GenomicRanges",
  "GenomicFeatures",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "rtracklayer"
))