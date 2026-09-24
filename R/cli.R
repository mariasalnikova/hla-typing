# ---------------------------------------------------------------------------
# Minimal command-line argument parsing (no extra dependencies).
# `spec` is a named list: option -> list(default, type = "character"|"numeric"|
# "integer"|"logical", help). Logical options are flags (`--flag`).
# ---------------------------------------------------------------------------

parse_cli_args <- function(args, spec, usage = "") {
  vals <- lapply(spec, `[[`, "default")
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) { cat(usage); print_cli_help(spec); quit(status = 0) }
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a, "\n", usage)
    key <- sub("^--", "", a)
    if (!key %in% names(spec)) stop("Unknown option --", key, "\n", usage)
    type <- spec[[key]]$type
    if (identical(type, "logical")) { vals[[key]] <- TRUE; i <- i + 1; next }
    if (i == length(args)) stop("Option --", key, " needs a value")
    v <- args[i + 1]
    vals[[key]] <- switch(type, numeric = as.numeric(v), integer = as.integer(v), v)
    i <- i + 2
  }
  vals
}

print_cli_help <- function(spec) {
  cat("\nOptions:\n")
  for (k in names(spec)) {
    d <- spec[[k]]$default
    cat(sprintf("  --%-18s %s%s\n", k, spec[[k]]$help,
                if (is.null(d) || is.na(d) || identical(d, FALSE)) "" else paste0(" [default: ", d, "]")))
  }
}

#' Sample name from an Illumina FASTQ file name
#' (`X_S65_L001_R1_001.fastq.gz` -> `X`, `X_R1.fastq.gz` -> `X`).
sample_name_from_fastq <- function(path) {
  b <- basename(path)
  b <- sub("\\.f(ast)?q(\\.gz)?$", "", b)
  b <- sub("_S[0-9]+_L[0-9]+_R[12]_001$", "", b)
  b <- sub("_R[12](_.*)?$", "", b)
  b
}
