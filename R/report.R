# ---------------------------------------------------------------------------
# Result tables: one row per called allele group, across samples.
# ---------------------------------------------------------------------------

#' Shorten an allele name to three fields (e.g. A*02:01:01:01 -> A*02:01:01).
allele_3fields <- function(allele) {
  f <- strsplit(allele, ":", fixed = TRUE)[[1]]
  paste0(f[seq_len(min(3, length(f)))], collapse = ":")
}

#' Collapse a space-separated list of alleles to the unique 3-field names.
tidy_allele_names <- function(allele_list) {
  paste(unique(vapply(strsplit(allele_list, " ", fixed = TRUE)[[1]], allele_3fields, character(1))),
        collapse = " ")
}

#' Combine stage-2 results of several samples into one table
#' (the "mega table" of the original script).
#'
#' @param results named list (sample -> type_hla() result); NULL entries are skipped.
#' @return data.frame with columns donor, tidyAllele, Allele, Score, n_noamps,
#'   sumreads, no_amps.
combine_typing_results <- function(results) {
  rows <- lapply(names(results), function(s) {
    calls <- results[[s]]$safety5
    if (is.null(calls) || nrow(calls) == 0) return(NULL)
    cbind(donor = s, calls, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab)) {
    return(data.frame(donor = character(), tidyAllele = character(), Allele = character(),
                      Score = numeric(), n_noamps = integer(), sumreads = numeric(),
                      no_amps = character(), stringsAsFactors = FALSE))
  }
  tab$tidyAllele <- vapply(tab$Allele, tidy_allele_names, character(1), USE.NAMES = FALSE)
  tab$n_noamps   <- lengths(strsplit(tab$no_amps, "   ", fixed = TRUE))
  rownames(tab) <- NULL
  tab[, c("donor", "tidyAllele", "Allele", "Score", "n_noamps", "sumreads", "no_amps")]
}

#' Write the typing table as tab-separated text.
write_typing_table <- function(tab, path) {
  utils::write.table(tab, path, sep = "\t", quote = TRUE, row.names = FALSE)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Human-readable summary of one sample's results (printed by the regression
# check so that one can see what the pipeline produces at every stage).
# ---------------------------------------------------------------------------

# "A*02:01:01:01 A*02:01:01:02L ..." -> "A*02:01 [372 alleles]"
.allele_brief <- function(alleles, max_groups = 2) {
  a <- strsplit(alleles, " ", fixed = TRUE)[[1]]
  a <- a[nzchar(a)]
  if (length(a) == 0) return("-")
  g <- unique(sub("^([^:]+:[^:]+).*", "\\1", a))
  paste0(paste(head(g, max_groups), collapse = ", "),
         if (length(g) > max_groups) sprintf(" (+%d)", length(g) - max_groups) else "",
         sprintf(" [%d allele%s]", length(a), if (length(a) == 1) "" else "s"))
}

.print_table <- function(df, indent = "  ") {
  old <- options(width = 10000)   # never wrap table columns
  on.exit(options(old))
  out <- utils::capture.output(print(df, row.names = FALSE, right = FALSE))
  cat(paste0(indent, out), sep = "\n")
}

#' Print what the pipeline produced for one sample: stage-1 read counts and
#' matched alleles per amplicon, stage-2 candidate allele groups per locus,
#' and the final calls.
#'
#' @param stage1,stage2 results of assemble_amplicons() and type_hla().
#' @param top_candidates how many candidate groups to show per locus.
print_results_summary <- function(sample, stage1, stage2, top_candidates = 5) {
  st <- stage1$statistics
  amps <- names(stage1$safety2)
  cat("\n----- Results for ", sample, " -----\n", sep = "")

  # ---- stage 1 --------------------------------------------------------------
  cat("\nStage 1: reads and sequences per amplicon\n",
      "  reads = read pairs with this primer; groups = distinct read-1 sequences kept;\n",
      "  trusted = sequences passed to typing (safety2); matched = their reads with an\n",
      "  exact database match; alleles = allele groups of the trusted sequences, by reads\n", sep = "")
  alleles <- vapply(amps, function(a) {
    x <- stage1$safety2[[a]]
    x <- x[x$Exact != "", , drop = FALSE]
    if (nrow(x) == 0) return("-")
    x <- x[order(-x$readnumber), , drop = FALSE]
    x <- x[!duplicated(x$assembled), , drop = FALSE]
    g <- unique(vapply(x$Exact, function(e) sub("^([^:]+:[^:]+).*", "\\1", strsplit(e, " ", fixed = TRUE)[[1]][1]),
                       character(1), USE.NAMES = FALSE))
    paste0(paste(head(g, 4), collapse = ", "), if (length(g) > 4) sprintf(" (+%d)", length(g) - 4) else "")
  }, character(1))
  .print_table(data.frame(
    amplicon = amps, reads = st$amps_reads, groups = st$n_champs,
    trusted = st$n_intersected, trusted_reads = st$reads_intersected, matched = st$used,
    alleles = alleles, stringsAsFactors = FALSE))
  cat(sprintf("  total read pairs assigned to amplicons: %d\n", as.integer(st$n_reads)))

  # ---- stage 2 --------------------------------------------------------------
  cat("\nStage 2: candidate allele groups per locus (safety4)\n",
      "  seen = amplicons of the locus containing the group; Score before normalisation;\n",
      "  dominated = number of groups explaining the same data and more;\n",
      "  called = reported in the final table (not dominated, or seen in every amplicon)\n", sep = "")
  for (g in names(LOCUS_GROUPS)) {
    grp <- LOCUS_GROUPS[[g]]
    tab <- stage2$safety4[[grp$table]]
    if (is.null(tab) || nrow(tab) == 0) { cat("  ", g, ": no candidates\n", sep = ""); next }
    n_amp  <- length(grp$amplicons)
    missed <- lengths(strsplit(tab$no_amps, "   ", fixed = TRUE))
    called <- tab$degree == 0 | tab$no_amps == ""
    ord <- order(!called, -tab$Score, -tab$sumreads)
    show <- head(ord, max(top_candidates, sum(called)))
    cat(sprintf("  %s: %d candidate group%s, %d called\n", g, nrow(tab), if (nrow(tab) == 1) "" else "s", sum(called)))
    .print_table(data.frame(
      alleles = vapply(tab$Allele[show], .allele_brief, character(1), USE.NAMES = FALSE),
      seen = sprintf("%d/%d", n_amp - missed[show], n_amp),
      Score = signif(tab$Score[show], 3), sumreads = round(tab$sumreads[show], 1),
      dominated = tab$degree[show], called = ifelse(called[show], "yes", "no"),
      stringsAsFactors = FALSE), indent = "    ")
    if (nrow(tab) > length(show)) cat(sprintf("    ... %d more not called\n", nrow(tab) - length(show)))
  }

  # ---- final calls ----------------------------------------------------------
  calls <- stage2$safety5
  cat("\nFinal calls (safety5, the <sample>_typing.tsv table)\n",
      "  Score is normalised within the locus: ~0.5 each for two alleles, ~1 for one;\n",
      "  missing = amplicons of the locus that did not contain the group\n", sep = "")
  .print_table(data.frame(
    locus = sub("\\*.*", "", sub(" .*", "", calls$Allele)),
    alleles = vapply(calls$Allele, .allele_brief, character(1), USE.NAMES = FALSE),
    Score = signif(calls$Score, 3),
    missing = ifelse(calls$no_amps == "", "-", gsub("   ", ", ", calls$no_amps, fixed = TRUE)),
    sumreads = round(calls$sumreads, 1), stringsAsFactors = FALSE))
  invisible(NULL)
}
