# ---------------------------------------------------------------------------
# Saving / loading pipeline results.
#
# Per sample the pipeline writes into the output directory:
#   <sample>_stage1.rds    full stage-1 result (R object, exact)
#   <sample>_stage1.json   the same as JSON
#   <sample>_stage2.rds    stage-2 result (safety3/4/5)
#   <sample>_typing.tsv    final allele calls (tab separated)
# ---------------------------------------------------------------------------

write_stage1 <- function(res, outdir, sample, json = TRUE) {
  saveRDS(res, file.path(outdir, paste0(sample, "_stage1.rds")))
  if (json) {
    jsonlite::write_json(res, file.path(outdir, paste0(sample, "_stage1.json")),
                         digits = NA, pretty = TRUE)
  }
  invisible(res)
}

read_stage1 <- function(path) {
  if (grepl("\\.json$", path)) {
    x <- jsonlite::read_json(path, simplifyVector = TRUE)
    # empty amplicons are read back as empty lists; make them data.frames
    fix <- function(l) lapply(l, function(t) if (is.data.frame(t)) t else data.frame())
    x$safety1 <- fix(x$safety1)
    x$safety2 <- fix(x$safety2)
    x$statistics <- lapply(x$statistics, unlist)
    x
  } else {
    readRDS(path)
  }
}

write_stage2 <- function(res, outdir, sample) {
  saveRDS(res, file.path(outdir, paste0(sample, "_stage2.rds")))
  invisible(res)
}
