#!/usr/bin/env Rscript

extensions <- "\\.(R|py|sh|tsv|csv|fa|txt|md)$"
max_size_bytes <- 0.5 * 1024^2

# Detect files reached through symbolic-link directories.
is_inside_symlink <- function(path) {
  path <- sub("^\\./", "", path)
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]

  prefixes <- vapply(
    seq_along(parts),
    function(i) do.call(file.path, as.list(parts[seq_len(i)])),
    character(1)
  )

  any(nzchar(Sys.readlink(prefixes)))
}

# Stage deletions and modifications to already tracked files.
stage_tracked_status <- system2("git", c("add", "-u"))

if (stage_tracked_status != 0) {
  stop("Failed to stage tracked changes.")
}

# Find existing files with allowed extensions.
files <- list.files(
  path = ".",
  pattern = extensions,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

# Exclude files reached through symbolic links.
files <- files[
  !vapply(files, is_inside_symlink, logical(1))
]

# Exclude missing files and files at or above the size limit.
sizes <- file.info(files)$size
files <- files[
  !is.na(sizes) & sizes < max_size_bytes
]

# Stage each qualifying file.
for (file in files) {
  stage_status <- system2(
    "git",
    c("add", "--", shQuote(file))
  )

  if (stage_status != 0) {
    warning(sprintf("Could not stage: %s", file))
  }
}

# Show exactly what will be committed without opening a pager.
system2(
  "git",
  c("--no-pager", "status", "--short")
)

system2(
  "git",
  c("--no-pager", "diff", "--cached", "--stat")
)

# Request confirmation.
args <- commandArgs(trailingOnly = TRUE)

if (!"--yes" %in% args) {
  message(
    "\nPreview only: nothing was committed.",
    "\nIf everything above looks correct, run:",
    "\n\nRscript sync2.R --yes"
  )
  quit(status = 0)
}

# Commit with the current date and time.
commit_message <- format(
  Sys.time(),
  "%Y-%m-%d_%H-%M-%S"
)

commit_status <- system2(
  "git",
  c("commit", "-m", shQuote(commit_message))
)

if (commit_status != 0) {
  stop("Commit failed; nothing was pushed.")
}

# Integrate remote changes without creating a merge commit.
pull_status <- system2(
  "git",
  c("pull", "--rebase", "origin", "main")
)

if (pull_status != 0) {
  stop(
    paste(
      "Pull/rebase failed.",
      "Resolve the problem before pushing."
    )
  )
}

# Push the updated main branch.
push_status <- system2(
  "git",
  c("push", "-u", "origin", "main")
)

if (push_status != 0) {
  stop("Push failed.")
}

message("Commit and push completed successfully.")
