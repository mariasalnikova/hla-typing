#!/usr/bin/env Rscript
# Download hla_nuc.fasta from ANHIG/IMGTHLA into db/ (or a given directory).
#   Rscript scripts/download_hla_db.R [--dest db] [--release Latest]
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
HLA_PIPELINE_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(HLA_PIPELINE_ROOT, "R", "utils.R"))
source(file.path(HLA_PIPELINE_ROOT, "R", "reference_db.R"))
source(file.path(HLA_PIPELINE_ROOT, "R", "cli.R"))
spec <- list(
  dest    = list(default = file.path(HLA_PIPELINE_ROOT, "db"), type = "character", help = "destination directory"),
  release = list(default = "Latest", type = "character", help = "ANHIG/IMGTHLA git ref, e.g. Latest or v3.58.0-alpha")
)
opt <- parse_cli_args(commandArgs(trailingOnly = TRUE), spec, "Usage: Rscript scripts/download_hla_db.R [--dest DIR] [--release REF]\n")
download_hla_db(opt$dest, opt$release)
