# ---------------------------------------------------------------------------
# STAGE 2: mapping trusted amplicon sequences onto the allele database and
# choosing the alleles that explain the most amplicons.
#
# type_hla() is the equivalent of get_HLA_result() in the original script and
# returns list(safety3, safety4, safety5):
#   safety3  per locus group: allele x amplicon tables with (a) the summed
#            read share of matching sequences and (b) the indices of matching
#            sequences; plus the read counts per amplicon
#   safety4  candidate allele groups (alleles that are indistinguishable by
#            the amplicons are merged) with a score and dominance graph
#   safety5  final call table: Allele, Score, no_amps, sumreads
# ---------------------------------------------------------------------------

# Locus groups: which alleles of the database are typed by which amplicons.
# Table/graph names (`table`, `graph`) are kept from the original output.
LOCUS_GROUPS <- list(
  Iclass = list(class_pattern = "^(A|B|C)$", table = "DT1", graph = "gr1",
                allele_key = "HLA_allele_Iclass", index_key = "HLA_Iclass_amps", reads_key = "n_readsI",
                amplicons = c("Iamp1", "Iamp1_inv", "Iamp2", "Iamp2_inv",
                              "Iamp1alt", "Iamp1alt_inv", "Iamp2alt", "Iamp2alt_inv")),
  DQB = list(class_pattern = "^DQB", table = "DQB", graph = "grDQB",
             allele_key = "HLA_allele_DQB", index_key = "HLA_DQB_amps", reads_key = "n_readsDQB",
             amplicons = c("IIamp1_DQB", "IIamp1_DQB_inv", "IIamp2", "IIamp2_inv",
                           "IIamp1_DQBalt", "IIamp1_DQBalt_inv", "IIamp2_DQBalt", "IIamp2_DQBalt_inv")),
  DRB = list(class_pattern = "^DRB", table = "DRB", graph = "grDRB",
             allele_key = "HLA_allele_DRB", index_key = "HLA_DRB_amps", reads_key = "n_readsDRB",
             amplicons = c("IIamp1_Others", "IIamp1_Others_inv", "IIamp2", "IIamp2_inv",
                           "IIamp1_DRBalt", "IIamp1_DRBalt_inv", "IIamp2_DRBalt", "IIamp2_DRBalt_inv")),
  DPB = list(class_pattern = "^DPB", table = "DPB", graph = "grDPB",
             allele_key = "HLA_allele_DPB", index_key = "HLA_DPB_amps", reads_key = "n_readsDPB",
             amplicons = c("IIamp1_DPB", "IIamp1_DPB_inv", "IIamp2_DPB", "IIamp2_DPB_inv")),
  DQA = list(class_pattern = "^DQA", table = "DQA", graph = "grDQA",
             allele_key = "HLA_allele_DQA", index_key = "HLA_DQA_amps", reads_key = "n_readsDQA",
             amplicons = c("IIamp_DQA", "IIamp_DQA_inv")),
  DPA = list(class_pattern = "^DPA", table = "DPA", graph = "grDPA",
             allele_key = "HLA_allele_DPA", index_key = "HLA_DPA_amps", reads_key = "n_readsDPA",
             amplicons = c("IIamp_DPA", "IIamp_DPA_inv"))
)

# Loci whose scores are normalised separately in the final table.
LOCUS_SCORE_PATTERNS <- list(
  list(pattern = "A*", fixed = TRUE), list(pattern = "B*", fixed = TRUE), list(pattern = "C*", fixed = TRUE),
  list(pattern = "DQB", fixed = TRUE), list(pattern = "DRB1", fixed = TRUE), list(pattern = "DRB[2-9]", fixed = FALSE),
  list(pattern = "DPB", fixed = FALSE), list(pattern = "DQA", fixed = FALSE), list(pattern = "DPA", fixed = FALSE)
)

#' Run stage 2 (typing) on a stage-1 result.
#'
#' @param stage1     output of assemble_amplicons().
#' @param db         allele table from load_hla_db().
#' @param min_freq   trusted sequences must hold at least this share of the
#'                   amplicon's reads (on top of the stage-1 filter).
#' @param pad_n      number of 'N' added to both ends of every database
#'                   sequence, so that amplicons overhanging short (partial)
#'                   allele sequences still match.
#' @param p_miss     probability that an amplicon of a true allele is missed;
#'                   Score = prod over amplicons of (1 - p_miss) if seen, p_miss if not.
#' @return list(safety3, safety4, safety5), or NULL if no class I amplicon has
#'   trusted sequences.
type_hla <- function(stage1, db, min_freq = 0.01, pad_n = 150, p_miss = 0.05) {
  class_I <- LOCUS_GROUPS$Iclass$amplicons
  if (sum(vapply(stage1$safety2[class_I], nrow, integer(1))) == 0) {
    log_msg("Stage 2: no trusted class I sequences, sample skipped")
    return(NULL)
  }
  trusted <- lapply(stage1$safety2, function(x) {
    if (nrow(x) == 0) return(x)
    x[x$parents %in% 1:4 & x$freq > min_freq & !grepl("NULL", x$assembled, fixed = TRUE), , drop = FALSE]
  })
  log_msg("Stage 2: mapping trusted sequences onto the allele database")
  safety3 <- map_sequences_to_alleles(trusted, db, pad_n = pad_n)
  safety4 <- score_allele_groups(safety3, p_miss = p_miss)
  safety5 <- call_alleles(safety4)
  log_msg("Stage 2 done: ", nrow(safety5), " allele calls")
  list(safety3 = safety3, safety4 = safety4, safety5 = safety5)
}

#' safety3: for every locus group, which alleles are matched by which trusted
#' sequences.
#'
#' For each allele and amplicon two things are accumulated over the trusted
#' sequences whose (N-tolerant) exact match includes the allele:
#'   * the summed read share `freq` of those sequences  -> HLA_allele_<group>
#'   * the row indices of those sequences, e.g. "0,1,3"  -> HLA_<group>_amps
#'     (the leading "0" is the initial value and is kept for compatibility)
map_sequences_to_alleles <- function(trusted, db, pad_n = 150) {
  pad <- strrep("N", pad_n)
  db$Sequence <- paste0(pad, db$Sequence, pad)
  out <- list()
  for (g in names(LOCUS_GROUPS)) {
    grp <- LOCUS_GROUPS[[g]]
    alleles <- db[grepl(grp$class_pattern, db$HLA_class), c("Allele", "Sequence")]
    rownames(alleles) <- NULL
    subject <- Biostrings::DNAStringSet(alleles$Sequence)
    freq_tab <- alleles
    idx_tab  <- alleles["Allele"]
    n_reads  <- setNames(rep(0, length(grp$amplicons)), grp$amplicons)
    for (amp in grp$amplicons) {
      freq_tab[[amp]] <- 0
      idx_tab[[amp]]  <- "0"
      seqs <- trusted[[amp]]
      if (is.null(seqs) || nrow(seqs) == 0) next
      for (i in seq_len(nrow(seqs))) {
        hit <- Biostrings::vcountPattern(seqs$assembled[i], subject,
                                         max.mismatch = 0, with.indels = FALSE, fixed = FALSE) != 0
        freq_tab[[amp]][hit] <- freq_tab[[amp]][hit] + seqs$freq[i]
        idx_tab[[amp]][hit]  <- paste0(idx_tab[[amp]][hit], ",", i)
      }
      n_reads[amp] <- sum(seqs$readnumber)
    }
    out[[grp$allele_key]] <- freq_tab
    out[[grp$index_key]]  <- idx_tab
    out[[grp$reads_key]]  <- n_reads
  }
  out
}

# TRUE when every amplicon index set of `child` is contained in that of `parent`
.dominates <- function(parent, child) {
  all(mapply(function(p, c) length(setdiff(c, p)) == 0, parent, child))
}

#' safety4: merge alleles with identical amplicon evidence, score them and
#' build the dominance graph.
#'
#' Per locus group, alleles sharing the same index signature `amps` are
#' collapsed into one row (`Allele` lists them all). Then:
#'   sumreads  sum over amplicons of (read share x amplicon read count)
#'   Score     prod over amplicons of (1 - p_miss) if matched else p_miss
#'   no_amps   amplicons in which the allele group was not seen
#'   degree    number of OTHER groups whose evidence is a superset of this
#'             group's evidence (in-degree in the dominance graph). 0 means
#'             nothing explains the data better.
score_allele_groups <- function(safety3, p_miss = 0.05) {
  out <- list()
  for (g in names(LOCUS_GROUPS)) {
    grp  <- LOCUS_GROUPS[[g]]
    amps <- grp$amplicons
    freq_tab <- safety3[[grp$allele_key]]
    idx_tab  <- safety3[[grp$index_key]]
    n_reads  <- safety3[[grp$reads_key]]
    freq_tab$amps    <- apply(idx_tab[amps], 1, paste0, collapse = "_")
    freq_tab$no_amps <- apply(idx_tab[amps], 1, function(r) paste0(amps[r == "0"], collapse = "   "))

    dt <- data.table::as.data.table(freq_tab)
    zero_sig <- paste0(rep("0", length(amps)), collapse = "_")
    tab <- dt[amps != zero_sig,
              c(list(Allele = paste0(Allele, collapse = " ")),
                lapply(.SD, unique),
                list(no_amps = unique(no_amps))),
              by = amps, .SDcols = amps]
    tab <- as.data.frame(tab, stringsAsFactors = FALSE)

    freq_mat <- as.matrix(tab[amps])
    tab$sumreads <- as.numeric(freq_mat %*% n_reads[amps])
    tab[amps] <- ifelse(freq_mat > 0, 1 - p_miss, p_miss)
    tab$Score <- apply(tab[amps], 1, prod)

    # dominance graph: edge i -> j when group i's evidence covers group j's
    sig <- lapply(strsplit(tab$amps, "_", fixed = TRUE),
                  function(s) strsplit(s, ",", fixed = TRUE))
    n <- length(sig)
    adj <- matrix(0, n, n)
    for (i in seq_len(n)) for (j in seq_len(n)) adj[i, j] <- .dominates(sig[[i]], sig[[j]])
    graph <- igraph::simplify(igraph::graph_from_adjacency_matrix(adj))
    tab$degree <- igraph::degree(graph, mode = "in")

    out[[grp$table]] <- tab
    out[[grp$graph]] <- graph
  }
  out
}

#' safety5: final calls. A group is reported when it is not dominated by any
#' other group (degree 0) or when it was seen in every amplicon. Scores are
#' normalised to sum to 1 within each locus.
call_alleles <- function(safety4) {
  parts <- lapply(LOCUS_GROUPS, function(grp) {
    tab <- safety4[[grp$table]]
    tab[tab$degree == 0 | tab$no_amps == "", c("Allele", "Score", "no_amps", "sumreads"), drop = FALSE]
  })
  res <- do.call(rbind, parts)
  rownames(res) <- NULL
  for (p in LOCUS_SCORE_PATTERNS) {
    sel <- grepl(p$pattern, res$Allele, fixed = p$fixed)
    if (any(sel)) res$Score[sel] <- prop.table(res$Score[sel])
  }
  res
}
