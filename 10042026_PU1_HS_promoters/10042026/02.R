Sys.setenv(TZ = "UTC")

library(dplyr)
library(rtracklayer)

genes <- import("genes.gtf")

tss <- promoters(
    genes,
    upstream = 0,
    downstream = 1
)

tss$name <- paste(
    tss$gene_name,
    tss$gene_id,
    tss$transcript_id,
    sep = "|"
)

length(tss)
table(width(tss))

tss[tss$gene_name == "Xkr4"]

tss_bed <- tss

mcols(tss_bed) <- DataFrame(
    name = paste(
        tss$gene_name,
        tss$gene_id,
        tss$transcript_id,
        sep = "|"
    ),
    score = rep(0L, length(tss))
)

head(mcols(tss_bed))
anyNA(tss_bed$score)

export(tss_bed, "genes_tss.bed", format = "BED")




