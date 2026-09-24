#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# HLA typing pipeline: FASTQ -> assembled amplicons -> allele calls.
#
#   Rscript scripts/run_pipeline.R --r1 S_R1.fastq.gz --r2 S_R2.fastq.gz [--sample S] [--outdir results]
#   Rscript scripts/run_pipeline.R --samplesheet samples.tsv [--cores 4] [--outdir results]
#   Rscript scripts/run_pipeline.R --stage1 results/S_stage1.rds --sample S   # re-run typing only
#
# samples.tsv is tab separated with a header line and columns: sample, r1, r2.
# ---------------------------------------------------------------------------
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
HLA_PIPELINE_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(HLA_PIPELINE_ROOT, "R", "load.R"))

usage <- "Usage:
  Rscript scripts/run_pipeline.R --r1 R1.fastq.gz --r2 R2.fastq.gz [--sample NAME] [options]
  Rscript scripts/run_pipeline.R --samplesheet samples.tsv [--cores N] [options]
  Rscript scripts/run_pipeline.R --stage1 NAME_stage1.rds --sample NAME [options]
"
spec <- list(
  r1          = list(default = NA, type = "character", help = "read 1 FASTQ(.gz)"),
  r2          = list(default = NA, type = "character", help = "read 2 FASTQ(.gz)"),
  sample      = list(default = NA, type = "character", help = "sample name (default: derived from the R1 file name)"),
  samplesheet = list(default = NA, type = "character", help = "TSV with columns sample, r1, r2 (many samples)"),
  stage1      = list(default = NA, type = "character", help = "existing *_stage1.rds/.json: skip stage 1 and only run typing"),
  db          = list(default = file.path(HLA_PIPELINE_ROOT, "db", "hla_nuc.fasta"), type = "character", help = "IPD-IMGT/HLA hla_nuc.fasta"),
  outdir      = list(default = "results", type = "character", help = "output directory"),
  cores       = list(default = 1L, type = "integer", help = "parallel processes for a samplesheet (fork based, not on Windows)"),
  threshold   = list(default = 100, type = "numeric", help = "max read-1 clusters kept per amplicon"),
  `read-length`     = list(default = 250, type = "numeric", help = "last read position kept (reads are cut to primer-length..read-length)"),
  `primer-length`   = list(default = 20, type = "numeric", help = "first read position kept; 20 keeps positions 20..250, i.e. removes 19 nt, exactly as the original scripts"),
  `min-read-length` = list(default = 200, type = "numeric", help = "discard read pairs with a shorter read"),
  downsample  = list(default = NA, type = "numeric", help = "randomly keep only this many read pairs"),
  seed        = list(default = NA, type = "integer", help = "random seed for --downsample"),
  `merge-method`    = list(default = "overlap", type = "character", help = "overlap | fixed (how read 1 and read 2 are joined)"),
  `min-reads`       = list(default = 10, type = "numeric", help = "stage 1: min reads for a trusted sequence"),
  `min-freq`        = list(default = 0.01, type = "numeric", help = "stage 2: min read share for a trusted sequence"),
  `p-miss`          = list(default = 0.05, type = "numeric", help = "stage 2: probability of missing an amplicon (scoring)"),
  `no-json`         = list(default = FALSE, type = "logical", help = "do not write the stage-1 JSON file")
)
opt <- parse_cli_args(commandArgs(trailingOnly = TRUE), spec, usage)

# ---- samples ---------------------------------------------------------------
if (!is.na(opt$samplesheet)) {
  samples <- utils::read.delim(opt$samplesheet, stringsAsFactors = FALSE, colClasses = "character")
  missing <- setdiff(c("sample", "r1", "r2"), names(samples))
  if (length(missing)) stop("samplesheet lacks columns: ", paste(missing, collapse = ", "))
} else if (!is.na(opt$stage1)) {
  if (is.na(opt$sample)) opt$sample <- sub("_stage1\\.(rds|json)$", "", basename(opt$stage1))
  samples <- data.frame(sample = opt$sample, r1 = NA, r2 = NA, stringsAsFactors = FALSE)
} else if (!is.na(opt$r1) && !is.na(opt$r2)) {
  if (is.na(opt$sample)) opt$sample <- sample_name_from_fastq(opt$r1)
  samples <- data.frame(sample = opt$sample, r1 = opt$r1, r2 = opt$r2, stringsAsFactors = FALSE)
} else {
  cat(usage); print_cli_help(spec); quit(status = 1)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
db <- load_hla_db(opt$db)
amplicons <- amplicon_table()

# ---- stage 1 ---------------------------------------------------------------
if (!is.na(opt$stage1)) {
  stage1 <- setNames(list(read_stage1(opt$stage1)), samples$sample)
} else {
  stage1 <- assemble_amplicons_batch(
    samples, db = db, cores = opt$cores,
    threshold = opt$threshold, read_length = opt$`read-length`,
    primer_length = opt$`primer-length`, min_read_length = opt$`min-read-length`,
    downsample = opt$downsample, merge_method = opt$`merge-method`,
    amplicons = amplicons, min_reads = opt$`min-reads`,
    seed = if (is.na(opt$seed)) NULL else opt$seed
  )
  for (s in names(stage1)) {
    if (inherits(stage1[[s]], "try-error")) next
    write_stage1(stage1[[s]], opt$outdir, s, json = !opt$`no-json`)
  }
}

# ---- stage 2 ---------------------------------------------------------------
typing <- list()
for (s in names(stage1)) {
  if (inherits(stage1[[s]], "try-error")) next
  log_msg("Typing sample ", s)
  typing[[s]] <- type_hla(stage1[[s]], db, min_freq = opt$`min-freq`, p_miss = opt$`p-miss`)
  if (is.null(typing[[s]])) next
  write_stage2(typing[[s]], opt$outdir, s)
  write_typing_table(combine_typing_results(typing[s]), file.path(opt$outdir, paste0(s, "_typing.tsv")))
}
summary_tab <- combine_typing_results(typing)
write_typing_table(summary_tab, file.path(opt$outdir, "typing_summary.tsv"))
log_msg("Done. ", nrow(summary_tab), " allele calls for ", length(typing), " sample(s) written to ", opt$outdir)
