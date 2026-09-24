# ---------------------------------------------------------------------------
# Assembling an amplicon from read 1 and the read-2 consensus.
#
# All functions take/return the per-amplicon "cluster table" produced by
# summarise_clusters() (columns read1, read2cons, readnumber) and add:
#   Overlap_max_aligned  length of the longest identical run in the overlap
#                        (or the assumed overlap length for fixed merging)
#   Overlap_start        1-based position in read 1 where read 2 (reverse
#                        complemented) starts
#   assembled            the merged amplicon sequence. Inside the overlap the
#                        read-1 bases are used because read 1 is of higher
#                        quality.
# ---------------------------------------------------------------------------

#' Locate where sequence `b` starts inside sequence `a` (one-sided overlap).
#'
#' For every offset i the suffix a[i..] is compared position by position with
#' the prefix of b and the longest run of identical nucleotides is recorded.
#' The offset with the longest run wins.
#' @return c(longest_identical_run, offset)
find_read_overlap <- function(a, b) {
  a <- NT_CODE[strsplit(a, "", fixed = TRUE)[[1]]]
  b <- NT_CODE[strsplit(b, "", fixed = TRUE)[[1]]]
  n <- length(a)
  res <- a
  for (i in seq_len(n)) {
    runs   <- rle(a[i:n] - b[1:(n - i + 1)])
    res[i] <- suppressWarnings(max(runs$lengths[runs$values == 0], na.rm = TRUE))
  }
  c(max(res), which(res == max(res))[1])
}

# read 1 + the part of revcomp(read 2) that extends beyond the overlap
.join_reads <- function(tab) {
  rc2 <- revcomp(tab$read2cons)
  paste0(tab$read1, substr(rc2, nchar(tab$read1) - tab$Overlap_start + 2, nchar(tab$read2cons)))
}

#' Merge read 1 and read 2 using the overlap found by find_read_overlap().
merge_by_overlap <- function(tab, ...) {
  if (nrow(tab) == 0) return(tab)
  ov <- t(vapply(seq_len(nrow(tab)),
                 function(i) find_read_overlap(tab$read1[i], revcomp(tab$read2cons[i])),
                 numeric(2)))
  tab$Overlap_max_aligned <- ov[, 1]
  tab$Overlap_start       <- ov[, 2]
  tab$assembled           <- .join_reads(tab)
  tab
}

#' Merge read 1 and read 2 assuming a fixed overlap of `shift` nucleotides
#' between trimmed reads of length `rlength`.
merge_fixed_overlap <- function(tab, shift = 1, rlength = 231) {
  if (nrow(tab) == 0) return(tab)
  tab$Overlap_max_aligned <- shift
  tab$Overlap_start       <- rlength - shift + 1
  tab$assembled           <- .join_reads(tab)
  tab
}

#' Join non-overlapping reads with `insert` 'N' characters in between.
merge_without_overlap <- function(tab, insert = 1) {
  if (nrow(tab) == 0) return(tab)
  tab$Overlap_max_aligned <- insert
  tab$Overlap_start       <- nchar(tab$read1[1]) + insert
  tab$assembled           <- paste0(tab$read1, strrep("N", insert), revcomp(tab$read2cons))
  tab
}

#' Apply the merging rule of one amplicon (see amplicon_table()).
merge_amplicon <- function(tab, spec, merge_method = c("overlap", "fixed"), rlength = 231) {
  merge_method <- match.arg(merge_method)
  switch(spec$merge_mode,
    standard = if (merge_method == "overlap") merge_by_overlap(tab)
               else merge_fixed_overlap(tab, shift = spec$shift, rlength = rlength),
    overlap_plus_fixed = {
      shifts <- as.numeric(strsplit(spec$fix_shifts, ",", fixed = TRUE)[[1]])
      do.call(rbind, c(list(merge_by_overlap(tab)),
                       lapply(shifts, function(s) merge_fixed_overlap(tab, shift = s, rlength = rlength))))
    },
    nonoverlap = merge_without_overlap(tab, insert = 1),
    stop("Unknown merge_mode: ", spec$merge_mode)
  )
}

#' Exact lookup of every assembled sequence in the allele database.
#'
#' 'N' in the assembly matches any nucleotide. Adds the columns
#'   Exact      alleles whose sequence contains the assembly (space separated)
#'   Exact_rev  the same for the reverse complement of the assembly
#'
#' Implementation: Biostrings::vcountPattern with IUPAC matching, which is
#' several times faster than the regular-expression search of the original
#' script and gives identical results as long as the database contains only
#' A/C/G/T (checked; otherwise the regex search is used).
match_assemblies_exact <- function(tab, db) {
  if (nrow(tab) == 0) return(tab)
  subject <- attr(db, "dna_set")
  if (is.null(subject)) subject <- Biostrings::DNAStringSet(db$Sequence)
  use_regex <- any(grepl("[^ACGT]", db$Sequence))
  lookup <- function(seqs) vapply(seqs, function(s) {
    hit <- if (use_regex) {
      grepl(gsub("N", "[ATGC]", s, fixed = TRUE), db$Sequence)
    } else {
      Biostrings::vcountPattern(s, subject, fixed = FALSE) != 0
    }
    paste0(db$Allele[hit], collapse = " ")
  }, character(1), USE.NAMES = FALSE)
  tab$Exact     <- lookup(tab$assembled)
  tab$Exact_rev <- lookup(revcomp(tab$assembled))
  tab
}

#' For inverse amplicons: report the assembly in the forward orientation.
straighten_inverse <- function(tab) {
  if (nrow(tab) == 0) return(tab)
  tab$Exact     <- tab$Exact_rev
  tab$assembled <- revcomp(tab$assembled)
  tab
}
