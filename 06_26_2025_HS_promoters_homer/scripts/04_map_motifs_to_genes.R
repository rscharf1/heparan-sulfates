library(data.table)
library(dplyr)


# Load file
hits <- read.delim("outputs/motif_hits.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 1. Extract gene identifier
gene_col <- "Gene.Name"
motif_cols <- setdiff(names(hits), names(hits)[1:21])  # First 21 columns are metadata

# 2. Convert motif columns to binary (1 = motif hit present, 0 = not present)
# Each motif cell has a string like: 145(ACGTGTGT,+,0.00),...
# So we check whether each cell is non-empty
motif_matrix <- as.data.frame(sapply(hits[, motif_cols], function(x) as.integer(nchar(x) > 0)))

# 3. Add gene names as rownames or first column
motif_matrix$Gene <- hits[[1]]

# 4. Optional: Reorder so Gene is first column
motif_matrix <- motif_matrix[, c("Gene", setdiff(names(motif_matrix), "Gene"))]

# 5. (Optional) View top entries
colnames(motif_matrix) <- gsub("..Homer.Distance.From.Peak.sequence.strand.conservation.", "", colnames(motif_matrix))
motif_matrix[1:5,1:5]

clean_names <- function(col) {
  if (col == "Gene") return(col)
  strsplit(col, "\\.")[[1]][1]  # Take just the first part before first "."
}

colnames(motif_matrix) <- sapply(colnames(motif_matrix), clean_names)

