# ---------------------------------------------------------------------------
# IPD-IMGT/HLA reference database: download and load `hla_nuc.fasta`.
#
# The pipeline maps assembled amplicons onto the nucleotide (cDNA/CDS)
# sequences of every known HLA allele. The database is distributed by the
# ANHIG/IMGTHLA GitHub repository; a specific release can be pinned by
# passing e.g. `release = "v3.58.0-alpha"` (see the repository tags).
# ---------------------------------------------------------------------------

IMGTHLA_RAW_URL <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/%s/%s"

#' Download `hla_nuc.fasta` (and the matching release_version.txt) into `dest_dir`.
#'
#' @param dest_dir Directory where the files are stored.
#' @param release  Git ref in ANHIG/IMGTHLA: "Latest" or a tag such as "v3.58.0-alpha".
#' @return Path to the downloaded FASTA file (invisibly).
download_hla_db <- function(dest_dir = "db", release = "Latest") {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  fasta <- file.path(dest_dir, "hla_nuc.fasta")
  log_msg("Downloading hla_nuc.fasta (", release, ") to ", fasta)
  utils::download.file(sprintf(IMGTHLA_RAW_URL, release, "hla_nuc.fasta"), fasta, mode = "wb", quiet = TRUE)
  utils::download.file(sprintf(IMGTHLA_RAW_URL, release, "release_version.txt"),
                       file.path(dest_dir, "release_version.txt"), mode = "wb", quiet = TRUE)
  log_msg("Database version: ", hla_db_version(dest_dir))
  invisible(fasta)
}

#' Read the IPD-IMGT/HLA version string stored next to the FASTA file, if any.
hla_db_version <- function(db_dir) {
  f <- file.path(db_dir, "release_version.txt")
  if (!file.exists(f)) return(NA_character_)
  v <- grep("^#\\s*version:", readLines(f, warn = FALSE), value = TRUE)
  if (length(v) == 0) NA_character_ else trimws(sub("^#\\s*version:\\s*", "", v[1]))
}

#' Load `hla_nuc.fasta` into the allele table used by the pipeline.
#'
#' FASTA headers look like `>HLA:HLA00005 A*01:01:01:01 1098 bp`.
#'
#' @return data.frame with columns Fasta_Id, Allele, length, Sequence, HLA_class
#'   (HLA_class is the locus part of the allele name, e.g. "A", "DRB1", "DQB1").
load_hla_db <- function(fasta_path) {
  if (!file.exists(fasta_path)) {
    stop("HLA database not found: ", fasta_path,
         "\nRun `Rscript scripts/download_hla_db.R` first.")
  }
  ss  <- Biostrings::readDNAStringSet(fasta_path)
  hdr <- stringr::str_split(names(ss), stringr::fixed(" "), simplify = TRUE)
  db <- data.frame(
    Fasta_Id  = hdr[, 1],
    Allele    = hdr[, 2],
    length    = hdr[, 3],
    Sequence  = as.character(ss),
    HLA_class = stringr::str_split(hdr[, 2], stringr::fixed("*"), simplify = TRUE)[, 1],
    stringsAsFactors = FALSE
  )
  rownames(db) <- NULL
  attr(db, "version") <- hla_db_version(dirname(fasta_path))
  attr(db, "dna_set") <- ss   # DNAStringSet of the same sequences, reused by stage 1
  log_msg("Loaded HLA database: ", nrow(db), " alleles",
          if (!is.na(attr(db, "version"))) paste0(" (", attr(db, "version"), ")") else "")
  db
}
