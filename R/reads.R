# ---------------------------------------------------------------------------
# Reading FASTQ files, clustering identical reads and building consensuses.
# ---------------------------------------------------------------------------

#' Read a pair of (gzipped) FASTQ files into a two-column data.frame.
#'
#' Only the sequence lines are kept; qualities are ignored. If the two files
#' have a different number of lines (e.g. after adapter trimming) the longer
#' one is truncated to the length of the shorter one.
read_fastq_pair <- function(r1_path, r2_path) {
  for (p in c(r1_path, r2_path)) if (!file.exists(p)) stop("FASTQ file not found: ", p)
  l1 <- readLines(r1_path)
  l2 <- readLines(r2_path)
  n  <- min(length(l1), length(l2))
  if (n < 4) stop("Empty FASTQ input: ", r1_path)
  seq_lines <- seq(2, n, by = 4)
  data.frame(read1 = l1[seq_lines], read2 = l2[seq_lines], stringsAsFactors = FALSE)
}

#' Trim reads and cluster read pairs with an identical read 1.
#'
#' Reads are cut to positions `primer_length`..`read_length` (1-based,
#' inclusive, i.e. `read_length - primer_length + 1` nt remain). Only read-1
#' sequences observed at least twice are kept, and among those the `threshold`
#' most abundant clusters with at least 3 stored copies.
#'
#' NOTE (kept from the original implementation): `duplicated()` drops the
#' first occurrence of every sequence, so the reported `readnumber` of a
#' cluster is one less than the true number of read pairs in it.
#'
#' @return Named list of data.frames (one per cluster, name = read-1 sequence),
#'   sorted alphabetically by read 1.
cluster_identical_reads <- function(reads, threshold = 100, read_length = 250, primer_length = 20) {
  reads$read1 <- substr(reads$read1, primer_length, read_length)
  reads$read2 <- substr(reads$read2, primer_length, read_length)
  reads  <- reads[duplicated(reads$read1), , drop = FALSE]
  groups <- split(reads, f = reads$read1)
  if (length(groups) > 0) {
    sizes  <- vapply(groups, nrow, integer(1))
    groups <- groups[(rank(-sizes, ties.method = "first") <= threshold) & (sizes > 2)]
  }
  groups
}

#' Majority-vote consensus of equal-length sequences.
#' A position gets the most frequent nucleotide if its share exceeds `thres`,
#' otherwise 'N'.
consensus_sequence <- function(seqs, thres = 0.7) {
  m <- do.call(rbind, strsplit(seqs, "", fixed = TRUE))
  cons <- vapply(seq_len(ncol(m)), function(j) {
    counts <- table(m[, j])
    if (max(prop.table(counts)) > thres) names(counts)[counts == max(counts)][1] else "N"
  }, character(1))
  paste0(cons, collapse = "")
}

#' Collapse each read-1 cluster into one row: read 1, consensus of its read 2s,
#' and the number of stored read pairs.
summarise_clusters <- function(groups) {
  data.frame(
    read1      = names(groups),
    read2cons  = vapply(groups, function(g) consensus_sequence(g$read2), character(1), USE.NAMES = FALSE),
    readnumber = vapply(groups, nrow, integer(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}
