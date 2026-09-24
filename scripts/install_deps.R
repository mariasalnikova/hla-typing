#!/usr/bin/env Rscript
# Install every R package the pipeline needs (CRAN + Bioconductor).
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L))
cran <- c("igraph", "stringdist", "stringr", "data.table", "dplyr", "tidyr", "jsonlite", "BiocManager")
need <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need)
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
ok <- vapply(c(cran, "Biostrings"), requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("Some packages failed to install: ", paste(names(ok)[!ok], collapse = ", "))
