#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Regression test: run the pipeline on the local test samples and compare
# every stage with the reference outputs in tests/expected/. Exit status 1 on
# any mismatch.
#
#   Rscript scripts/check_regression.R [--outdir DIR] [--cores N] [--skip-run] [--quiet]
#
# After the comparison of each sample it prints a summary of what the
# pipeline produced (reads per amplicon, candidate alleles, final calls);
# --quiet leaves only the PASS/FAIL lines. The output files of the run stay
# in --outdir (default results/regression) for closer inspection.
#
# A test sample consists of
#   data/<sample>_*_R1_001.fastq.gz, data/<sample>_*_R2_001.fastq.gz   input
#   tests/expected/<sample>_stage1.rds                                  reference
#   tests/expected/<sample>_stage2.rds
#   tests/expected/<sample>_typing.tsv
# The reference files are the pipeline's own output files for that sample
# (see README, section Testing). The test data are donor data and are not
# part of the repository.
#
# The reference outputs are produced with IPD-IMGT/HLA 3.65.0 and allele
# names differ between releases, so the check always runs on that release.
# It is downloaded once into db/test/, independently of db/hla_nuc.fasta
# used for real analyses (which may be at the latest release).
# ---------------------------------------------------------------------------
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
HLA_PIPELINE_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(HLA_PIPELINE_ROOT, "R", "load.R"))

spec <- list(
  outdir     = list(default = file.path(HLA_PIPELINE_ROOT, "results", "regression"), type = "character", help = "where to write the run"),
  cores      = list(default = 1L, type = "integer", help = "samples processed in parallel"),
  `skip-run` = list(default = FALSE, type = "logical", help = "compare an existing run in --outdir"),
  quiet      = list(default = FALSE, type = "logical", help = "do not print the result summaries")
)
opt <- parse_cli_args(commandArgs(trailingOnly = TRUE), spec,
                      "Usage: Rscript scripts/check_regression.R [--outdir DIR] [--cores N] [--skip-run] [--quiet]\n")

TEST_DB_RELEASE <- "v3.65.0-alpha"          # git tag in ANHIG/IMGTHLA
TEST_DB_VERSION <- "IPD-IMGT/HLA 3.65.0"    # as written in release_version.txt
test_db_dir <- file.path(HLA_PIPELINE_ROOT, "db", "test")
test_db     <- file.path(test_db_dir, "hla_nuc.fasta")

data_dir <- file.path(HLA_PIPELINE_ROOT, "data")
exp_dir  <- file.path(HLA_PIPELINE_ROOT, "tests", "expected")
samples  <- sub("_stage1\\.rds$", "", list.files(exp_dir, pattern = "_stage1\\.rds$"))
if (length(samples) == 0) {
  stop("no test data: put <sample>_stage1.rds, _stage2.rds and _typing.tsv into ", exp_dir,
       " and the FASTQ pairs into ", data_dir, " (see README, section Testing)")
}
fastq <- function(s, read) {
  f <- list.files(data_dir, pattern = paste0("^", s, "_.*_R", read, "_001\\.fastq(\\.gz)?$"), full.names = TRUE)
  if (length(f) != 1) stop("expected exactly one R", read, " FASTQ for ", s, " in data/, found ", length(f))
  f
}

if (!opt$`skip-run`) {
  if (!file.exists(test_db) || !identical(hla_db_version(test_db_dir), TEST_DB_VERSION)) {
    download_hla_db(test_db_dir, TEST_DB_RELEASE)
  }
  if (!identical(hla_db_version(test_db_dir), TEST_DB_VERSION)) {
    stop("test database in ", test_db_dir, " is not ", TEST_DB_VERSION)
  }
  sheet <- data.frame(sample = samples, r1 = vapply(samples, fastq, "", read = 1),
                      r2 = vapply(samples, fastq, "", read = 2))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  sheet_path <- file.path(opt$outdir, "samples.tsv")
  utils::write.table(sheet, sheet_path, sep = "\t", quote = FALSE, row.names = FALSE)
  args <- c(file.path(HLA_PIPELINE_ROOT, "scripts", "run_pipeline.R"),
            "--samplesheet", sheet_path, "--cores", opt$cores, "--outdir", opt$outdir,
            "--db", test_db)
  status <- system2("Rscript", args)
  if (status != 0) stop("pipeline run failed")
}

failures <- 0
check <- function(ok, what, detail = "") {
  cat(if (ok) "  PASS " else "  FAIL ", what, if (!ok && nzchar(detail)) paste0(": ", detail) else "", "\n", sep = "")
  if (!ok) failures <<- failures + 1
  invisible(ok)
}
# TRUE, or a short description of the first differences
same <- function(a, b) {
  r <- all.equal(a, b, tolerance = 1e-8, check.attributes = FALSE)
  if (isTRUE(r)) TRUE else paste(head(r, 3), collapse = "; ")
}
check_same <- function(a, b, what) {
  r <- same(a, b)
  check(isTRUE(r), what, if (isTRUE(r)) "" else r)
}

for (sample in samples) {
  cat("\n=== ", sample, " ===\n", sep = "")
  exp_path <- function(suffix) file.path(exp_dir, paste0(sample, suffix))
  new_path <- function(suffix) file.path(opt$outdir, paste0(sample, suffix))

  # ---- stage 1 -----------------------------------------------------------
  cat("Stage 1 (candidate and trusted sequences per amplicon)\n")
  exp1 <- readRDS(exp_path("_stage1.rds"))
  new1 <- readRDS(new_path("_stage1.rds"))
  check(identical(names(exp1$safety2), names(new1$safety2)), "amplicon names and order")
  for (amp in names(exp1$safety2)) {
    check_same(exp1$safety1[[amp]], new1$safety1[[amp]], paste0(amp, ": all candidate sequences (safety1)"))
    check_same(exp1$safety2[[amp]], new1$safety2[[amp]], paste0(amp, ": trusted sequences (safety2)"))
  }
  for (st in names(exp1$statistics)) {
    check_same(exp1$statistics[[st]], new1$statistics[[st]], paste0("statistics$", st))
  }

  # ---- stage 2 -----------------------------------------------------------
  cat("Stage 2 (allele calls)\n")
  exp2 <- readRDS(exp_path("_stage2.rds"))
  new2 <- readRDS(new_path("_stage2.rds"))
  for (tab in setdiff(names(exp2$safety3), "safety2_merged")) {
    check_same(exp2$safety3[[tab]], new2$safety3[[tab]], paste0("safety3$", tab))
  }
  for (tab in names(exp2$safety4)) {
    if (inherits(exp2$safety4[[tab]], "igraph")) next   # graphs are summarised by the `degree` column
    check_same(exp2$safety4[[tab]], new2$safety4[[tab]], paste0("safety4$", tab))
  }
  check_same(exp2$safety5, new2$safety5, "safety5 (final calls)")

  # ---- final table -------------------------------------------------------
  cat("Final typing table\n")
  exp3 <- utils::read.delim(exp_path("_typing.tsv"), stringsAsFactors = FALSE)
  new3 <- utils::read.delim(new_path("_typing.tsv"), stringsAsFactors = FALSE)
  check_same(exp3, new3, "typing table")

  if (!opt$quiet) print_results_summary(sample, new1, new2)
}

cat("\nOutput files of this run: ", normalizePath(opt$outdir), "\n", sep = "")
cat("\n", if (failures == 0) "ALL CHECKS PASSED" else paste(failures, "CHECK(S) FAILED"), "\n")
quit(status = if (failures == 0) 0 else 1)
