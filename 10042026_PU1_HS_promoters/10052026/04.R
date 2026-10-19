Sys.setenv(TZ = "UTC")

library(Biostrings)
library(TFBSTools)
library(motifmatchr)
library(data.table)
library(dplyr)

# Scan for PU.1 in all promoters 
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
    p.cutoff = 1e-4,
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

gc_counts <- letterFrequency(
    promoters,
    letters = c("G", "C"),
    as.prob = TRUE
)

scan_results$gc_fraction <- rowSums(gc_counts)

# Controls
pos_ctrls <- c("Csf1r", "Itgam", "Spi1", "Cybb", "Il1b", "Fcgr1", "Fcgr2b", "Cd14", "Cd68", "Lyz2", "Ncf2", "Mpo", "Elane", "Ctsg", "Tlr4", "Fcer1g", "Tyrobp")
neg_ctrls <- c("Alb", "Ttr", "Myh6", "Tnnt2", "Myod1", "Myog", "Ins1", "Ins2", "Pdx1", "Krt5", "Krt14", "Camk2a")

scan_results[gene_name %in% pos_ctrls]
scan_results[gene_name %in% neg_ctrls]


which(scan_results$gene_name == "Epcam")
promoters[19208] %>% as.character()

# Stats on results 

HS_genes <- read.csv(
    "../HS_relevant_genes.csv",
    check.names = FALSE
) %>% .$hgnc_symbol %>% unique()

test <- c("Ext1", "B3gat1")

lapply(test, function(gene) {
    print(gene)

    this_gc <- scan_results[gene_name == gene]$gc_fraction

    control_pool <- scan_results[!(tolower(gene_name) %in% tolower(HS_genes))][abs(gc_fraction - this_gc) < 0.01]

    
    n_controls <- nrow(control_pool)
    n_hits     <- sum(control_pool$has_pu1)
    p_hit      <- mean(control_pool$has_pu1)

    data.table(
        target_gene       = gene,
        target_has_pu1    = scan_results[gene_name == gene]$has_pu1,
        control_promoters = n_controls,
        control_hits      = n_hits,
        control_hit_rate  = p_hit
    )

    hit_rate_test <- binom.test(n_hits, n_controls)

    hit_rate_test$estimate
    hit_rate_test$conf.int
})


#########
# My own crude motif finder 

scan_results[tolower(gene_name) %in% tolower(HS_genes)][motif_count > 0]









