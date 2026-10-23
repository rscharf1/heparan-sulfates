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
