# ---------------------------------------------------------------------------
# Amplicon definitions.
#
# Every sequencing read pair is assigned to an amplicon by the primer sequence
# found at the very start of read 1 (and, for two amplicons, of read 2). Each
# amplicon also carries the rule used to merge read 1 with read 2.
#
# Naming convention (kept from the original scripts):
#   I*        class I amplicons (HLA-A/B/C);  II*  class II amplicons
#   amp1/amp2 first/second overlapping amplicon of a locus
#   *alt      alternative primer system for the same region
#   *_inv     the same amplicon read in the opposite orientation
#             (read 1 starts from the reverse primer). Such assemblies are
#             reverse-complemented after the database lookup so that all
#             amplicons are reported in the same orientation.
#   IIamp2 and IIamp2_DPB share the same read-1 primer; the reads are simply
#   processed twice, once for DRB/DQB and once for DPB typing.
#
# merge_mode:
#   "standard"           merge read 1 with the reverse complement of the read-2
#                        consensus. The overlap is found by sequence search
#                        (merge_method = "overlap") or set to a fixed length
#                        `shift` (merge_method = "fixed").
#   "overlap_plus_fixed" produce the sequence-search assembly AND one fixed-
#                        overlap assembly for every value in `fix_shifts`
#                        (rows are stacked; used for DQA/DPA/DPB whose
#                        overlaps are hard to locate).
#   "nonoverlap"         reads do not overlap: join them with a single 'N'.
#
# `shift` is the expected overlap length between read 1 and read 2 (in nt)
# for the trimmed 231-nt reads. It is only used by fixed-overlap merging.
# ---------------------------------------------------------------------------

amplicon_table <- function() {
  a <- function(name, group, read1_pattern, read1_prefix, inverse = FALSE,
                read2_pattern = NA_character_, read2_prefix = NA_integer_,
                merge_mode = "standard", shift = NA_real_, fix_shifts = "") {
    data.frame(name = name, group = group,
               read1_pattern = read1_pattern, read1_prefix = read1_prefix,
               read2_pattern = read2_pattern, read2_prefix = read2_prefix,
               inverse = inverse, merge_mode = merge_mode, shift = shift,
               fix_shifts = fix_shifts, stringsAsFactors = FALSE)
  }
  tab <- rbind(
    # ---- forward orientation --------------------------------------------
    a("Iamp1",          "I",   "CCCTGACC[GC]AGACCTG",                                      20, shift = 29),
    a("Iamp2",          "I",   "CGACGGCAA[AG]GATTAC",                                      20, shift = 9),
    a("IIamp1_DQB",     "DQB", "AG[GT]CTTTGCGGATCCC",                                      20, shift = 71),
    a("IIamp1_Others",  "DRB", "CTGAGCTCCC[GC]ACTGG",                                      20, shift = 130),
    a("IIamp2",         "DRB/DQB", "GGAACAGCCAGAAGGA",                                    20, shift = 155),
    a("Iamp1alt",       "I",   "TC[CT]CACTCCATGAGGTATTTC|TCCCACTCCATGAAGTATTTC",           22, merge_mode = "nonoverlap"),
    a("Iamp2alt",       "I",   "GGCAA[AG]GATTACATCGCC|GGCAAGGATTACATCGCT",                 20, shift = 40),
    a("IIamp1_DQBalt",  "DQB", "TGAGGGCAGAGAC[CT]CTCC",                                    22, shift = 32),
    a("IIamp1_DRBalt",  "DRB", "TGACAGTGACACTGATGG|TGACAGTGACATTGACGG",                    22, shift = 35),
    a("IIamp2_DRBalt",  "DRB", "GAGAGCTTCAC[AG]GTGCAG",                                    22, shift = 166),
    a("IIamp2_DQBalt",  "DQB", "ACCATCTCCCCATCCAG",                                        22, shift = 205),
    # ---- inverse orientation (read 1 starts from the other end) ----------
    a("Iamp1_inv",      "I",   "GGGCCGCCTCC[AC]ACTTG|GGGCCGTCTCCCACTTG|GGACCGCCTCCCACTTG", 20, TRUE, shift = 29),
    a("Iamp2_inv",      "I",   "[CT]GGTGG[AG]CTGGGAAGA",                                   20, TRUE, shift = 9),
    a("IIamp1_DQB_inv", "DQB", "[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG",               20, TRUE, "AG[GT]CTTTGCGGATCCC", 20, shift = 71),
    a("IIamp1_Others_inv", "DRB", "[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG",            20, TRUE, "CTGAGCTCCC[GC]ACTGG", 20, shift = 130),
    a("IIamp2_inv",     "DRB/DQB", "CCAC[GT]TGGCAGGTGTA|CCACTTGGCAAGTGTA",                 20, TRUE, shift = 155),
    a("Iamp1alt_inv",   "I",   "GAGC[GC]ACTCCACGCAC|GAGCCCGTCCACGCAC",                     22, TRUE, merge_mode = "nonoverlap"),
    a("Iamp2alt_inv",   "I",   "TCAGGGTGAGGGGCT|TCAGGGTGCAGGGCT",                          20, TRUE, shift = 40),
    a("IIamp1_DQBalt_inv", "DQB", "GTCCAGTCACC[AG]TTCCTA",                                 22, TRUE, shift = 32),
    a("IIamp1_DRBalt_inv", "DRB", "CAG[CT]CTTCTCTTCCTGGC",                                 22, TRUE, shift = 35),
    a("IIamp2_DRBalt_inv", "DRB", "TGCTCTGTGCAGATTCAG",                                    22, TRUE, shift = 166),
    a("IIamp2_DQBalt_inv", "DQB", "TG[CT]TCTGGGCAGATTCAG",                                 22, TRUE, shift = 205),
    # ---- DQA1, DPB1, DPA1 -------------------------------------------------
    a("IIamp_DQA",      "DQA", "ACAAAGCTCTG[AC]TGCTGGG",                                   24, merge_mode = "overlap_plus_fixed", shift = 33, fix_shifts = "33,36"),
    a("IIamp_DQA_inv",  "DQA", "AGAAACA[GC]CTTCTGTGACTG",                                  24, TRUE, merge_mode = "overlap_plus_fixed", shift = 33, fix_shifts = "33,36"),
    a("IIamp1_DPB",     "DPB", "GCGTTACTGATGGTGCTGC",                                      24, merge_mode = "overlap_plus_fixed", shift = 104, fix_shifts = "104"),
    a("IIamp1_DPB_inv", "DPB", "ATC[CT]GTCACGTGGCAGAC",                                    23, TRUE, merge_mode = "overlap_plus_fixed", shift = 104, fix_shifts = "104"),
    a("IIamp2_DPB",     "DPB", "GGAACAGCCAGAAGGA",                                         20, merge_mode = "overlap_plus_fixed", shift = 155, fix_shifts = "155"),
    a("IIamp2_DPB_inv", "DPB", "CCAC[GT]TGGCA[GA][GA]TGTA",                                24, TRUE, merge_mode = "overlap_plus_fixed", shift = 155, fix_shifts = "155"),
    a("IIamp_DPA",      "DPA", "CGGACCATGTGTCAACTTAT",                                     20, merge_mode = "overlap_plus_fixed", shift = 33, fix_shifts = "33,36"),
    a("IIamp_DPA_inv",  "DPA", "CTCTGCTGAGGGCACA",                                         24, TRUE, merge_mode = "overlap_plus_fixed", shift = 33, fix_shifts = "33,36")
  )
  rownames(tab) <- tab$name
  tab
}

#' Split a table of read pairs into one table per amplicon (see amplicon_table()).
#' A read pair is assigned to an amplicon when the first `read1_prefix`
#' nucleotides of read 1 match `read1_pattern` (and likewise for read 2 when a
#' read-2 pattern is given). A pair may fall into several amplicons.
demultiplex_reads <- function(reads, amplicons = amplicon_table()) {
  res <- lapply(seq_len(nrow(amplicons)), function(i) {
    spec <- amplicons[i, ]
    keep <- grepl(spec$read1_pattern, substr(reads$read1, 1, spec$read1_prefix))
    if (!is.na(spec$read2_pattern)) {
      keep <- keep & grepl(spec$read2_pattern, substr(reads$read2, 1, spec$read2_prefix))
    }
    reads[keep, , drop = FALSE]
  })
  names(res) <- amplicons$name
  res
}
