library(biomaRt)
library(dplyr)
library(readr)
library(data.table)

# Define gene list
# hs_genes <- fread("inputs/HS_relevant_genes.csv")

species <- "human"

genes <- fread("inputs/metabolic_pathway_genes_human.csv")
genes <- genes$cholesterol
genes <- genes[genes != ""]

bed <- genes_to_bed(genes, species)

write_tsv(bed, "inputs/cholesterol_promoters_human.bed", col_names = FALSE)

genes_to_bed <- function(genes, species) {
  # Connect to Ensembl
  ensembl <- if (species == "mouse") {
    useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  } else {
    useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }

  # Retrieve canonical transcript TSS and info
  gene_coords <- getBM(
    attributes = c("external_gene_name", "chromosome_name", "strand", 
                   "transcript_start", "transcript_end", "transcription_start_site", 
                   "transcript_is_canonical"),
    filters = "external_gene_name",
    values = genes,
    mart = ensembl
  ) %>%
    filter(transcript_is_canonical == 1) %>%
    filter(grepl("^\\d+$|^X$|^Y$", chromosome_name))  # keep standard chromosomes

  # Format as BED: 1000 bp upstream to 500 bp downstream (strand-aware)
  promoter_bed <- gene_coords %>%
    mutate(
      chrom = paste0("chr", chromosome_name),
      start = ifelse(strand == 1,
                     pmax(transcription_start_site - 1000, 0),
                     transcription_start_site - 500),
      end = ifelse(strand == 1,
                   transcription_start_site + 500,
                   transcription_start_site + 1000),
      name = external_gene_name
    ) %>%
    dplyr::select(chrom, start, end, name)

  return(promoter_bed)
}