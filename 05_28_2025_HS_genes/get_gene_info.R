library(dplyr)
library(biomaRt)
library(data.table)

enzymes <- c("Hs2st1", "Hs3st1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Ndst3", 
"Ndst2", "Ndst1", "Glce", "Ext2", "Extl2", "Extl1", "Extl3", 
"Ext1", "Gpc4", "Tgfbr3", "Agrn")

# Your list of gene symbols (example list)
gene_list <- c("GATA2", "SPI1", "RUNX1", "CEBPA", "TAL1")

# Set up the Ensembl connection (human genes)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Query biomaRt for gene information
gene_info <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "external_gene_name", "description", "gene_biotype"),
  filters = "hgnc_symbol",
  values = toupper(enzymes),
  mart = ensembl
)

# View results
print(gene_info)

# Optionally save to CSV
write.csv(gene_info, file = "HS_relevant_genes.csv", row.names = FALSE)