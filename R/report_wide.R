# ---------------------------------------------------------------------------
# STUDY-SPECIFIC post-processing (kept from the original script `hla_wide`).
#
# Turns a manually curated typing table (at most two allele groups per locus
# per donor) into one row per donor with columns <Locus>_1 / <Locus>_2, and
# flags the DR3 (DRB1*03:01-DQB1*02:01-DQA1*05:01) and DR4
# (DRB1*04-DQB1*03:02/03:01/02-DQA1*03:01) haplotypes used in the type 1
# diabetes study. It also contains hard-coded allele renamings that only
# make sense for that data set. Review before using on other projects.
# ---------------------------------------------------------------------------

#' @param path  tab-separated typing table (output of write_typing_table()),
#'              filtered by hand to <= 2 rows per locus and donor.
#' @return list(wide = one row per donor, long = one row per allele)
typing_table_to_wide <- function(path) {
  d <- utils::read.delim(path, stringsAsFactors = FALSE)
  d$Locus        <- vapply(strsplit(d$tidyAllele, "*", fixed = TRUE), `[[`, character(1), 1)
  d$FirstAllele  <- vapply(strsplit(d$tidyAllele, " ", fixed = TRUE), `[[`, character(1), 1)
  d$SimpleAllele <- stringr::str_sub(d$FirstAllele, 1, -4)
  long <- d[, c("donor", "Locus", "FirstAllele", "SimpleAllele", "tidyAllele")]

  wide <- d[, c("donor", "Locus", "FirstAllele")]
  # Data-set specific corrections of the original study
  wide$FirstAllele <- gsub("DPB1*107:01", "DPB1*13:01:01", wide$FirstAllele, fixed = TRUE)
  wide$FirstAllele <- gsub("B*07:100 ",   "B*15:01:01",    wide$FirstAllele, fixed = TRUE)
  wide$FirstAllele <- gsub("DPB1*138:01", "DQB1*23:01:01", wide$FirstAllele, fixed = TRUE)

  wide <- dplyr::mutate(dplyr::group_by(wide, donor), n = as.integer(duplicated(Locus)) + 1)
  wide$Locus <- paste0(wide$Locus, "_", wide$n)
  wide <- tidyr::spread(dplyr::ungroup(wide)[, c("donor", "Locus", "FirstAllele")], Locus, FirstAllele)
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)

  wide$genotype <- apply(wide[, -1, drop = FALSE], 1, paste0, collapse = "_")
  has <- function(pattern) grepl(pattern, wide$genotype)
  wide$DR3 <- has("DRB1[*]03:01") & has("DQB1[*]02:01") & has("DQA1[*]05:01")
  wide$DR4 <- has("DRB1[*]04:(01|02|04|05|08)") & (has("DQB1[*]03:(01|05)") | has("DQB1[*]02")) & has("DQA1[*]03:01")
  wide$same_genotype <- vapply(wide$genotype, function(g) {
    paste(wide$donor[grepl(g, wide$genotype, fixed = TRUE)], collapse = "_")
  }, character(1), USE.NAMES = FALSE)
  list(wide = wide, long = long)
}
