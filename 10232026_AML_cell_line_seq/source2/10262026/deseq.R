library(DESeq2)

fc <- read.delim(
    "counts/gene_counts.txt",
    comment.char = "#",
    row.names = 1,
    check.names = FALSE
)

# Remove annotation columns: Chr, Start, End, Strand, Length.
raw_counts <- as.matrix(fc[, -(1:5), drop = FALSE])

colnames(raw_counts) <- sub(
    "_ERR[0-9]+$",
    "",
    sub(
        "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
        "",
        basename(colnames(raw_counts))
    )
)

metadata <- data.frame(row.names = colnames(raw_counts))

dds <- DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = metadata,
    design = ~ 1
)

dds <- estimateSizeFactors(dds)
normalized <- counts(dds, normalized = TRUE)

save_table <- function(x, filename) {
    write.table(
        data.frame(gene_id = rownames(x), x, check.names = FALSE),
        filename,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}

save_table(raw_counts, "gene_raw_counts.tsv")
save_table(normalized, "gene_normalized_counts.tsv")
save_table(log2(normalized + 1), "gene_log2_normalized_counts.tsv")
