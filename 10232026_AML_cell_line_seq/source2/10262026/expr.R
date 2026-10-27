library(data.table)
library(dplyr)
library(rtracklayer)
library(ggplot2)

counts <- fread("gene_normalized_counts.tsv")

gtf_path <- "/gs/gsfs0/users/rscharf/Tools/star_files/index/gencode.v44.annotation.gtf"

# Read only feature type (column 3) and attributes (column 9).
gtf <- read.delim(
    gtf_path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE,
    colClasses = c(
        "NULL", "NULL", "character",
        rep("NULL", 5), "character"
    )
)

# Extract identifiers from gene records.
attributes <- gtf$V9[gtf$V3 == "gene"]

gene_map <- unique(data.table(
    gene_id = sub('.*gene_id "([^"]+)".*', '\\1', attributes),
    gene_name = sub('.*gene_name "([^"]+)".*', '\\1', attributes)
))

# Match without changing row order.
setDT(counts)
counts[, gene_name := gene_map$gene_name[
    match(gene_id, gene_map$gene_id)
]]

setcolorder(
    counts,
    c("gene_id", "gene_name",
      setdiff(names(counts), c("gene_id", "gene_name")))
)

counts[1:4, 1:6]


ext1 <- melt(
    counts[gene_name == "EXT1"],
    id.vars = c("gene_id", "gene_name"),
    variable.name = "cell_line",
    value.name = "normalized_count"
)

setorder(ext1, -normalized_count)

pdf("out1.pdf", height = 10)

	ggplot(ext1, aes(
	    x = reorder(cell_line, normalized_count),
	    y = log2(normalized_count + 1)
	)) +
	    geom_col(fill = "steelblue") +
	    coord_flip() +
	    labs(
	        title = "EXT1 expression across cell lines",
	        x = NULL,
	        y = "log2(normalized count + 1)"
	    ) +
	    theme_classic() +
	    theme(axis.text.y = element_text(size = 7))

dev.off()


hs_genes <- fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_genes_10182026.csv")

synthesis_genes <- unique(hs_genes[Category == "Synthesis_1", Gene])

synthesis_long <- melt(
    counts[gene_name %in% synthesis_genes],
    id.vars = c("gene_id", "gene_name"),
    variable.name = "cell_line",
    value.name = "normalized_count"
)

plots <- setNames(
    lapply(synthesis_genes, function(gene) {
        ggplot(
            synthesis_long[gene_name == gene],
            aes(
                x = reorder(cell_line, normalized_count),
                y = log2(normalized_count + 1)
            )
        ) +
            geom_col(fill = "steelblue") +
            coord_flip() +
            labs(
                title = paste(gene, "expression across cell lines"),
                x = NULL,
                y = "log2(normalized count + 1)"
            ) +
            theme_classic() +
            theme(axis.text.y = element_text(size = 7))
    }),
    synthesis_genes
)

pdf("out2.pdf", width = 8, height = 16)
	for (p in plots) print(p)
dev.off()



targets <- c("HL-60", "NB-4")

long <- melt(
    counts[gene_name %in% hs_genes$Gene],
    id.vars = c("gene_id", "gene_name"),
    variable.name = "cell_line",
    value.name = "normalized_count"
)

long[, `:=`(
    expression_rank = rank(
        -normalized_count,
        ties.method = "average",
        na.last = "keep"
    ),
    n_cell_lines = sum(!is.na(normalized_count))
), by = .(gene_id, gene_name)]

selected <- copy(long[cell_line %in% targets])

selected[, cell_line := fifelse(cell_line == "HL-60", "HL60", "NB4")]

rank_table <- dcast(
    selected,
    gene_id + gene_name + n_cell_lines ~ cell_line,
    value.var = c("normalized_count", "expression_rank")
)

rank_table <- merge(
    hs_genes,
    rank_table,
    by.x = "Gene",
    by.y = "gene_name",
    all.x = TRUE,
    sort = FALSE
)

rank_table[Category == "Glycoprotein"]
rank_table[Category == "Synthesis_1"]
rank_table[Category == "Synthesis_2"]
rank_table[Category == "Degradation"]
rank_table[Category == "Remodeling"]
rank_table[Category == "Sulfation"]
rank_table[Category == "Other"]
Degradation




comparison <- merge(
    unique(hs_genes[, .(Gene, Category)]),
    counts[, .(
        Gene = gene_name,
        gene_id,
        NB4 = `NB-4`,
        THP1 = `THP-1`
    )],
    by = "Gene",
    all.x = TRUE,
    sort = FALSE
)

comparison[, log2_ratio := log2((NB4 + 1) / (THP1 + 1))]

setorder(comparison, log2_ratio)

comparison[]


comparison[log2_ratio <= -2 & THP1 >= 20]

comparison[
    Gene %in% c(paste0("SDC", 1:4), paste0("GPC", 1:6)),
    .(Gene, NB4, THP1, log2_ratio)
][order(-NB4)]

comparison_long <- melt(
    comparison,
    id.vars = c("Gene", "Category", "gene_id"),
    measure.vars = c("NB4", "THP1"),
    variable.name = "cell_line",
    value.name = "normalized_count"
)

pdf("out3.pdf")

ggplot(
    comparison_long,
    aes(
        x = Gene,
        y = log2(normalized_count + 1),
        fill = cell_line
    )
) +
    geom_col(position = "dodge") +
    coord_flip() +
    facet_wrap(~ Category, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(
        NB4 = "#D55E00",
        THP1 = "#0072B2"
    )) +
    labs(
        title = "HS pathway expression: NB-4 versus THP-1",
        x = NULL,
        y = "log2(normalized count + 1)",
        fill = NULL
    ) +
    theme_classic()

dev.off()




