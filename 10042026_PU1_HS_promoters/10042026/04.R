Sys.setenv(TZ = "UTC")

library(Biostrings)
library(TFBSTools)
library(motifmatchr)
library(data.table)
library(dplyr)

promoters <- readDNAStringSet(
    "mouse_protein_coding_promoters.fa"
)

length(promoters)
head(names(promoters))

motif_lines <- readLines("../PU1_motif.jaspar")

motif_rows <- motif_lines[
    grepl("^[ACGT][[:space:]]", motif_lines)
]

pfm <- do.call(
    rbind,
    lapply(motif_rows, function(line) {
        numbers_only <- gsub("[^0-9]+", " ", line)
        scan(text = numbers_only, quiet = TRUE)
    })
)

rownames(pfm) <- substr(trimws(motif_rows), 1, 1)
colnames(pfm) <- seq_len(ncol(pfm))

pu1_pfm <- PFMatrix(
    ID = "MA0080.3",
    name = "Spi1",
    matrixClass = "ETS",
    strand = "+",
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    profileMatrix = pfm
)

pu1_scan <- matchMotifs(
    pu1_pfm,
    promoters,
    bg = "even",
    p.cutoff = 5e-5,
    out = "scores"
)

has_pu1 <- as.vector(motifMatches(pu1_scan)[, 1])
pu1_count <- as.vector(motifCounts(pu1_scan)[, 1])
pu1_best_score <- as.vector(motifScores(pu1_scan)[, 1])

scan_results <- data.table(
    promoter = names(promoters),
    has_pu1 = has_pu1,
    motif_count = pu1_count,
    best_score = pu1_best_score
)

scan_results$has_pu1 %>% table()
scan_results$motif_count %>% table()
scan_results[order(-motif_count)]

scan_results[grepl("Csf", promoter)]
scan_results[grepl("Epc", promoter)]

header <- sub(
    "\\([+-]\\)$",
    "",
    scan_results$promoter
)

header_parts <- strsplit(
    header,
    "|",
    fixed = TRUE
)

scan_results$gene_name <- vapply(
    header_parts,
    `[`,
    character(1),
    1
)

scan_results$gene_id <- vapply(
    header_parts,
    `[`,
    character(1),
    2
)

scan_results$transcript_id <- vapply(
    header_parts,
    `[`,
    character(1),
    3
)

HS_genes <- read.csv(
    "../HS_relevant_genes.csv",
    check.names = FALSE
) %>% .$hgnc_symbol %>% unique()


table(tolower(HS_genes) %in% tolower(scan_results$gene_name))

scan_results[tolower(gene_name) %in% tolower(HS_genes)]

# Background 
gc_counts <- letterFrequency(
    promoters,
    letters = c("G", "C"),
    as.prob = TRUE
)

scan_results$gc_fraction <- rowSums(gc_counts)

scan_results$is_target <-
    tolower(scan_results$gene_name) %in% tolower(HS_genes)

target_results <- scan_results[
    scan_results$is_target,
]

nrow(target_results)

target_results[
    ,
    c(
        "gene_name",
        "gene_id",
        "has_pu1",
        "motif_count",
        "best_score"
    )
]



targets <- as.data.frame(target_results)

background_pool <- as.data.frame(
    scan_results[!scan_results$is_target, ]
)

gene_tests <- vector(
    "list",
    nrow(targets)
)

for (i in seq_len(nrow(targets))) {

    gene <- targets[i, ]

    gc_difference <- abs(
        background_pool$gc_fraction -
        gene$gc_fraction
    )

    # Prefer controls within 1 percentage point of GC.
    control_indices <- which(gc_difference <= 0.01)

    # Ensure enough controls for a useful empirical p-value.
    if (length(control_indices) < 500) {
        control_indices <- order(gc_difference)[
            seq_len(min(500, nrow(background_pool)))
        ]
    }

    controls <- background_pool[control_indices, ]

    gene_tests[[i]] <- data.frame(
        gene_name = gene$gene_name,
        gene_id = gene$gene_id,
        gc_fraction = gene$gc_fraction,

        motif_count = gene$motif_count,
        best_score = gene$best_score,

        n_controls = nrow(controls),
        maximum_gc_difference =
            max(abs(controls$gc_fraction - gene$gc_fraction)),

        control_median_count =
            median(controls$motif_count),

        count_percentile =
            100 * mean(controls$motif_count < gene$motif_count),

        count_p_value =
            (1 + sum(controls$motif_count >= gene$motif_count)) /
            (1 + nrow(controls)),

        score_percentile =
            100 * mean(controls$best_score < gene$best_score),

        score_p_value =
            (1 + sum(controls$best_score >= gene$best_score)) /
            (1 + nrow(controls))
    )
}

gene_tests <- do.call(
    rbind,
    gene_tests
)
