# ---------------------------------------------------------------------------
# STAGE 1: from raw FASTQ to assembled, error-filtered amplicon sequences.
#
# assemble_amplicons() is the equivalent of HLA_amplicones_full() in the
# original script. Its result is a list:
#   safety1     one data.frame per amplicon with ALL read-1 clusters that
#               passed clustering (candidate sequences with graph statistics)
#   safety2     the same tables restricted to the trusted ("parent")
#               sequences that are passed on to typing
#   statistics  read counts at every step, per amplicon
# ---------------------------------------------------------------------------

#' Run stage 1 for one sample.
#'
#' @param r1_path,r2_path  FASTQ(.gz) files with read 1 and read 2.
#' @param db               allele table from load_hla_db().
#' @param threshold        keep at most this many read-1 clusters per amplicon.
#' @param read_length      reads are cut at this position (250 for 2x250 or
#'                         2x300 runs; the fixed overlap lengths in
#'                         amplicon_table() assume this value).
#' @param primer_length    reads are cut from this position, removing the primer.
#' @param min_read_length  read pairs with a shorter read are discarded.
#' @param downsample       if not NA, randomly keep this many read pairs.
#' @param merge_method     "overlap": locate the read-1/read-2 overlap by
#'                         sequence search; "fixed": use the amplicon's `shift`.
#' @param amplicons        amplicon definitions, see amplicon_table().
#' @param min_reads        clusters with fewer stored reads are not trusted.
#' @param min_freq_secondary  a sequence ranked 2..4 in its graph component
#'                         is trusted only if its share of reads exceeds this.
#' @param seed             random seed used for downsampling.
assemble_amplicons <- function(r1_path, r2_path, db,
                               threshold = 100, read_length = 250, primer_length = 20,
                               min_read_length = 200, downsample = NA,
                               merge_method = c("overlap", "fixed"),
                               amplicons = amplicon_table(),
                               min_reads = 10, min_freq_secondary = 0.05,
                               seed = NULL) {
  merge_method <- match.arg(merge_method)
  rlength <- read_length - primer_length + 1   # length of a trimmed read

  # 1. read the FASTQ pair --------------------------------------------------
  log_msg("Stage 1: reading ", basename(r1_path))
  reads <- read_fastq_pair(r1_path, r2_path)
  if (!is.na(downsample) && downsample > 0 && downsample < nrow(reads)) {
    if (!is.null(seed)) set.seed(seed)
    reads <- reads[sample(seq_len(nrow(reads)), size = downsample, replace = FALSE), ]
  }
  reads <- reads[nchar(reads$read1) > min_read_length & nchar(reads$read2) > min_read_length, ]
  log_msg("  read pairs after length filter: ", nrow(reads))

  # 2. assign read pairs to amplicons by primer -----------------------------
  by_amp    <- demultiplex_reads(reads, amplicons)
  amp_reads <- vapply(by_amp, nrow, integer(1))
  all_reads <- sum(amp_reads)
  log_msg("  read pairs assigned to amplicons: ", all_reads)
  rm(reads)

  # 3. cluster identical read 1, consensus of read 2 ------------------------
  log_msg("  clustering identical reads")
  by_amp <- lapply(by_amp, function(x) {
    summarise_clusters(cluster_identical_reads(x, threshold = threshold,
                                               read_length = read_length,
                                               primer_length = primer_length))
  })
  n_champ     <- vapply(by_amp, nrow, integer(1))
  champ_reads <- vapply(by_amp, function(x) sum(x$readnumber), numeric(1))

  # 4. merge read 1 with read 2 into the amplicon sequence ------------------
  log_msg("  assembling read pairs")
  for (nm in names(by_amp)) {
    by_amp[[nm]] <- merge_amplicon(by_amp[[nm]], amplicons[nm, ], merge_method, rlength)
  }

  # 5. exact database lookup (both orientations) ----------------------------
  log_msg("  exact matching of assemblies against ", nrow(db), " alleles")
  for (nm in names(by_amp)) by_amp[[nm]] <- match_assemblies_exact(by_amp[[nm]], db)

  # 6. graph of 1-mismatch neighbours: which sequences are "parents" --------
  log_msg("  building mismatch graphs")
  by_amp <- lapply(by_amp, rank_parents)
  inverse <- amplicons$name[amplicons$inverse]
  by_amp[inverse] <- lapply(by_amp[inverse], straighten_inverse)

  # 7. keep trusted sequences only -----------------------------------------
  trusted <- lapply(by_amp, function(x) {
    if (nrow(x) == 0) return(x)
    keep <- x$readnumber > (min_reads - 1) &
      (x$parents == 1 | (x$parents %in% 2:4 & x$freq > min_freq_secondary))
    x[keep, , drop = FALSE]
  })
  n_trusted     <- vapply(trusted, nrow, integer(1))
  trusted_reads <- vapply(trusted, function(x) sum(x$readnumber), numeric(1))
  used          <- vapply(trusted, function(x) sum(x$readnumber[x$Exact != ""]), numeric(1))

  statistics <- list(
    n_reads = all_reads, amps_reads = amp_reads,
    n_champs = n_champ, champ_reads = champ_reads,
    n_intersected = n_trusted, reads_intersected = trusted_reads,
    used = used, non_used = trusted_reads - used
  )
  log_msg("Stage 1 done: ", sum(n_trusted), " trusted sequences in ",
          sum(n_trusted > 0), "/", length(trusted), " amplicons")
  list(safety1 = by_amp, safety2 = trusted, statistics = statistics)
}

#' Run stage 1 for several samples, optionally in parallel (forking; not on Windows).
#'
#' @param samples data.frame with columns `sample`, `r1`, `r2`.
#' @param cores   number of parallel processes.
#' @param ...     passed to assemble_amplicons().
#' @return named list of stage-1 results.
assemble_amplicons_batch <- function(samples, db, cores = 1, ...) {
  run_one <- function(i) try(assemble_amplicons(samples$r1[i], samples$r2[i], db = db, ...))
  idx <- seq_len(nrow(samples))
  res <- if (cores > 1) parallel::mclapply(idx, run_one, mc.cores = cores) else lapply(idx, run_one)
  names(res) <- samples$sample
  failed <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(failed)) warning("Stage 1 failed for: ", paste(names(res)[failed], collapse = ", "))
  res
}
