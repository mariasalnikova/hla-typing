# Agent notes for the HLA typing pipeline

This repository is an R pipeline that calls HLA alleles from paired-end
Illumina amplicon sequencing of HLA cDNA (loci A, B, C, DRB1/3/4/5, DQB1,
DQA1, DPB1, DPA1). Read `README.md` for usage and `docs/algorithm.md` for what
the code does step by step before changing anything.

## Layout

- `R/` – pipeline modules, sourced together by `R/load.R` (no R package
  structure; plain `source()`).
  - `reference_db.R` download/load IPD-IMGT/HLA `hla_nuc.fasta`
  - `amplicons.R` the amplicon/primer table (`amplicon_table()`), demultiplexing
  - `reads.R`, `assembly.R`, `graph_filter.R` stage-1 building blocks
  - `stage1_assemble.R` `assemble_amplicons()` – FASTQ to trusted amplicon sequences
  - `stage2_typing.R` `type_hla()` – map sequences to alleles, score, call
  - `report.R` per-sample and multi-sample result tables, and the readable
    summary printed by the regression check; `report_wide.R` is a
    study-specific post-processing kept from the original code
  - `io.R` output files, `cli.R` argument parsing
- `scripts/` – command-line entry points (`run_pipeline.R`,
  `download_hla_db.R`, `install_deps.R`, `check_regression.R`)
- `data/`, `tests/expected/` – test samples (FASTQ) and the pipeline's
  reference outputs for them (produced with IPD-IMGT/HLA 3.65.0). They are
  real donor data: git-ignored, present only locally, never in the
  repository. See the Testing section of `README.md` for the file names.
- `db/` – the reference database (git-ignored, downloaded by a script);
  `db/test/` holds the pinned release 3.65.0 used only by the regression check

## Rules

- The pipeline must keep reproducing `tests/expected/`. The current
  reference was produced by this pipeline after it had been verified to
  reproduce the original analysis scripts on the same samples; the only
  intended difference from them is that `IIamp_DPA_inv` is reverse-complemented
  like every other inverse amplicon. After any change to
  `R/` run `Rscript scripts/check_regression.R --cores 3` (takes several
  minutes; it runs the whole pipeline on every test sample). If the local
  test data are missing, say that the change could not be verified.
- Never commit sequencing data or typing results of real donors (`data/`,
  `tests/expected/`, `results/`, any FASTQ, RDS or typing tables), and never
  put donor genotypes into documentation, comments or commit messages.
- Several quirks of the original algorithm are deliberately preserved because
  the reference outputs depend on them; they are marked with `NOTE` comments
  (e.g. `readnumber` is one less than the true count, index strings start with
  a "0", primers are trimmed with an off-by-one). Do not "fix" them silently;
  add an option and keep the default behaviour, or regenerate the local
  `tests/expected/` (see README, Testing) and say so in the commit message.
- Output objects keep the historical names `safety1..safety5`, `DT1`, `gr1`
  etc. because downstream analyses of users read them.
- The amplicon order in `amplicon_table()` matters: it defines the order of
  the `statistics` vectors and of the output JSON.
- Dependencies are installed with `scripts/install_deps.R`; do not add new
  packages without updating that script and the README.
- Do not commit the reference database, PDFs or any pipeline output.
- Keep documentation in `README.md` (usage) and `docs/algorithm.md`
  (technical description) in sync with the code.

## Running

```bash
Rscript scripts/install_deps.R
Rscript scripts/download_hla_db.R
Rscript scripts/run_pipeline.R --r1 SAMPLE_R1.fastq.gz --r2 SAMPLE_R2.fastq.gz --outdir results
Rscript scripts/check_regression.R --cores 3
```
