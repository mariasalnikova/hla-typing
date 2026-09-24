# Source every module of the pipeline. Usage from the repository root:
#   source("R/load.R")
suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(igraph)
  library(stringdist)
  library(stringr)
  library(dplyr)
})
hla_pipeline_dir <- if (exists("HLA_PIPELINE_ROOT")) HLA_PIPELINE_ROOT else {
  # locate this file when it is source()d, so that the modules are found
  # regardless of the working directory
  f <- NULL
  for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) break }
  if (is.null(f)) "." else dirname(dirname(normalizePath(f)))
}
for (m in c("utils", "reference_db", "amplicons", "reads", "assembly", "graph_filter",
            "stage1_assemble", "stage2_typing", "report", "report_wide", "io", "cli")) {
  source(file.path(hla_pipeline_dir, "R", paste0(m, ".R")))
}
