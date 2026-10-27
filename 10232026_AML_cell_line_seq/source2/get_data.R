library(dplyr)
library(tidyr)

metadata_url <- paste0(
  "https://www.ebi.ac.uk/ena/portal/api/filereport?",
  "accession=PRJEB30312",
  "&result=read_run",
  "&fields=run_accession,sample_accession,sample_alias,",
  "sample_title,experiment_title,fastq_ftp",
  "&format=tsv"
)

metadata <- read.delim(
  metadata_url,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Display names alongside run IDs
metadata[, c(
  "run_accession", "sample_accession",
  "sample_alias", "sample_title", "experiment_title"
)]

download_fastqs <- function(metadata, outdir = "fastq") {
  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (i in seq_len(nrow(metadata))) {
    paths <- metadata$fastq_ftp[i]
    if (is.na(paths) || !nzchar(paths)) next

    cell_line <- gsub("[^A-Za-z0-9_-]", "_", metadata$sample_title[i])
    sample_dir <- file.path(outdir, cell_line)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

    for (path in strsplit(paths, ";", fixed = TRUE)[[1]]) {
      path <- trimws(path)
      if (!nzchar(path)) next

      url <- paste0("https://", sub("^(ftp|https?)://", "", path))
      destination <- file.path(sample_dir, basename(path))
      partial <- paste0(destination, ".part")

      if (file.exists(destination)) {
        message("Already exists: ", destination)
        next
      }

      message("Downloading: ", cell_line, " / ", basename(path))

      status <- download.file(
        url,
        destfile = partial,
        method = "libcurl",
        mode = "wb"
      )

      if (status != 0L) {
        stop("Download failed: ", url)
      }

      if (!file.rename(partial, destination)) {
        stop("Could not rename downloaded file: ", partial)
      }
    }
  }

  invisible(NULL)
}

selected <- metadata %>%
  filter(sample_title %in% c("HL-60", "MOLM-13", "THP-1"))

download_fastqs(metadata, outdir = "fastq")





download_table <- metadata %>%
  select(sample_title, fastq_ftp) %>%
  separate_rows(fastq_ftp, sep = ";") %>%
  mutate(
    fastq_ftp = trimws(fastq_ftp),
    cell_line = gsub("[^A-Za-z0-9_-]", "_", sample_title),
    url = paste0("https://", sub("^(ftp|https?)://", "", fastq_ftp)),
    filename = basename(fastq_ftp)
  ) %>%
  filter(!is.na(fastq_ftp), nzchar(fastq_ftp)) %>%
  distinct(cell_line, url, filename)

write.table(
  download_table,
  "fastq_downloads.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)