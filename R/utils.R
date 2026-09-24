# ---------------------------------------------------------------------------
# Small helpers shared by both pipeline stages.
# ---------------------------------------------------------------------------

#' Print a timestamped progress message to the console.
log_msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), paste0(...), "\n", sep = " ")
}

#' Reverse-complement a character vector of DNA sequences (A/C/G/T/N only).
revcomp <- function(seqs) {
  rc <- c(A = "T", T = "A", C = "G", G = "C", N = "N")
  vapply(strsplit(seqs, "", fixed = TRUE),
         function(chars) paste0(rev(rc[chars]), collapse = ""),
         character(1), USE.NAMES = FALSE)
}

#' Number of 'N' characters in a single sequence.
n_content <- function(seq) {
  sum(strsplit(seq, "", fixed = TRUE)[[1]] == "N")
}

#' Numeric code for each nucleotide, used by the overlap search.
NT_CODE <- c(A = 1, G = 2, C = 3, T = 4, N = 5)
