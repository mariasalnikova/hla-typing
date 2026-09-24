# HLA typing from cDNA amplicon sequencing

An R pipeline that determines the HLA genotype of a donor from paired-end
Illumina reads of overlapping HLA cDNA amplicons (HLA-A, -B, -C, -DRB1/3/4/5,
-DQB1, -DQA1, -DPB1, -DPA1). It was developed for a custom library
preparation system in which every locus is covered by short overlapping
amplicons read in both orientations, and it types alleles by exact,
error-tolerant matching against the IPD-IMGT/HLA database rather than by
conventional read mapping.

The method is described in section 4.1 of A. A. Minervina's PhD thesis
("Разработка алгоритма анализа данных секвенирования для HLA-типирования").
The original script is located https://github.com/asya-minervina/HLA_old.
This repository is a restructured, documented version of the updated
analysis scripts.

- `docs/algorithm.md` – what every step does, in technical terms

## Requirements

- R >= 4.1 (tested with 4.5)
- CRAN: igraph, stringdist, stringr, data.table, dplyr, tidyr, jsonlite
- Bioconductor: Biostrings

```bash
Rscript scripts/install_deps.R
```

## Reference database

The pipeline needs `hla_nuc.fasta` (nucleotide sequences of all HLA alleles)
from [ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA). Download the latest
release into `db/`:

```bash
Rscript scripts/download_hla_db.R
```

or pin a release: `Rscript scripts/download_hla_db.R --release v3.58.0-alpha`.
The release used is stored in `db/release_version.txt`. Allele names in the
results depend on the database release, so record it with your results.

## Input

Paired-end FASTQ files (gzipped or not) from a MiSeq 2x250 or 2x300 run of the
HLA amplicon libraries. Reads must start with the amplification primer
(untrimmed), because the primer at the start of read 1 is what assigns a read
pair to an amplicon. About 50 000 read pairs per sample are sufficient.

## Running

One sample:

```bash
Rscript scripts/run_pipeline.R \
  --r1 fastq/D01_S1_L001_R1_001.fastq.gz \
  --r2 fastq/D01_S1_L001_R2_001.fastq.gz \
  --outdir results
```

Many samples, four in parallel (Linux/macOS only):

```bash
Rscript scripts/run_pipeline.R --samplesheet samples.tsv --cores 4 --outdir results
```

where `samples.tsv` is tab-separated with a header and three columns:

```
sample	r1	r2
D01	fastq/D01_S1_L001_R1_001.fastq.gz	fastq/D01_S1_L001_R2_001.fastq.gz
D02	fastq/D02_S2_L001_R1_001.fastq.gz	fastq/D02_S2_L001_R2_001.fastq.gz
```

Re-run only the typing stage (for example after updating the database):

```bash
Rscript scripts/run_pipeline.R --stage1 results/D01_stage1.rds --outdir results
```

`Rscript scripts/run_pipeline.R --help` lists all options. The important ones:

| option | default | meaning |
|---|---|---|
| `--db` | `db/hla_nuc.fasta` | reference database |
| `--threshold` | 100 | max number of distinct read-1 sequences kept per amplicon |
| `--read-length` / `--primer-length` | 250 / 20 | reads keep positions 20..250 (the first 19 nt, the primer, are removed; these are the original values, do not change them to reproduce old results) |
| `--min-read-length` | 200 | read pairs with a shorter read are dropped |
| `--downsample N --seed S` | off | analyse a random subset of N read pairs |
| `--merge-method` | overlap | `overlap`: find the read1/read2 overlap by sequence search; `fixed`: use the per-amplicon expected overlap |
| `--min-reads` | 10 | stage 1: a sequence needs this many reads to be trusted |
| `--min-freq` | 0.01 | stage 2: a sequence needs this share of the amplicon's reads |
| `--p-miss` | 0.05 | scoring: probability that an amplicon of a true allele is missing |

A sample with ~50 000 read pairs takes a few minutes on a laptop; almost all
of the time is spent matching sequences against the database.

## Output

Per sample, in `--outdir`:

| file | content |
|---|---|
| `<sample>_typing.tsv` | **final allele calls** (see below) |
| `<sample>_stage2.rds` | R list `safety3`, `safety4`, `safety5` with all intermediate tables of the typing stage |
| `<sample>_stage1.rds` | R list `safety1` (all candidate sequences), `safety2` (trusted sequences), `statistics` (read counts per amplicon) |
| `<sample>_stage1.json` | the same as JSON (skip with `--no-json`) |

plus `typing_summary.tsv`, the calls of all samples in one table.

Columns of the typing table:

| column | meaning |
|---|---|
| `donor` | sample name |
| `tidyAllele` | the called allele group, names shortened to 3 fields (`A*02:01:01`) |
| `Allele` | all database alleles that are indistinguishable with these amplicons (full names) |
| `Score` | support of the group, normalised to sum to 1 within the locus. Two alleles of a heterozygous locus get ~0.5 each, a homozygous one ~1 |
| `n_noamps`, `no_amps` | how many / which amplicons of the locus did not contain the group. 0 means every amplicon supports it |
| `sumreads` | reads supporting the group, summed over amplicons |

A locus is normally reported as one or two rows with high `Score` and
`n_noamps = 0`. Additional rows with a tiny score (e.g. 4e-6) and several
missing amplicons are weakly supported alternatives and should be discarded
by the analyst; the original workflow filtered the summary table by hand to
at most two allele groups per locus before further analysis
(`typing_table_to_wide()` in `R/report_wide.R`).

## Using the functions from R

```r
source("R/load.R")
db  <- load_hla_db("db/hla_nuc.fasta")
s1  <- assemble_amplicons("D01_R1.fastq.gz", "D01_R2.fastq.gz", db)   # stage 1
s2  <- type_hla(s1, db)                                                 # stage 2
s2$safety5                                                              # calls
combine_typing_results(list(D01 = s2))                                  # summary table
```

## Testing

The regression check `scripts/check_regression.R` runs the pipeline on a set
of test samples and compares every stage with reference outputs saved
earlier for the same samples, so any change in behaviour is detected.

The test data are not included in the repository, because they are
sequencing data of real donors. To run the check, put them locally (both
directories are git-ignored):

- `data/` – the FASTQ pair of every test sample, named
  `<sample>_*_R1_001.fastq.gz` and `<sample>_*_R2_001.fastq.gz`
- `tests/expected/` – the reference outputs for each sample: the files
  `<sample>_stage1.rds`, `<sample>_stage2.rds` and `<sample>_typing.tsv`
  written by `run_pipeline.R` (see *Output*)

Then run (`--cores N` processes the samples in parallel):

```bash
Rscript scripts/check_regression.R --cores 3
```

For every sample the check prints PASS/FAIL lines and then a summary of what
the pipeline produced, so the results can be inspected by eye:

- stage 1: read pairs, read-1 groups and trusted sequences per amplicon, and
  the allele groups the trusted sequences match;
- stage 2: the best candidate allele groups per locus with the number of
  amplicons they were seen in, their score and whether they were called;
- the final calls, as in `<sample>_typing.tsv`.

`--quiet` prints only the PASS/FAIL lines. All output files of the run stay
in `results/regression/` (change with `--outdir`).

Allele names change between IPD-IMGT/HLA releases, so the check always uses
release 3.65.0: on the first run it downloads it (~40 MB) into `db/test/`,
separately from `db/hla_nuc.fasta` used for real analyses, which can stay at
the latest release.

To add a test sample, run the pipeline on it with the test database and copy
the three files into `tests/expected/` after checking the calls by hand:

```bash
Rscript scripts/run_pipeline.R --r1 data/D01_S1_L001_R1_001.fastq.gz \
  --r2 data/D01_S1_L001_R2_001.fastq.gz --sample D01 \
  --db db/test/hla_nuc.fasta --outdir /tmp/D01 --no-json
cp /tmp/D01/D01_stage1.rds /tmp/D01/D01_stage2.rds /tmp/D01/D01_typing.tsv tests/expected/
```

The same procedure regenerates the reference when the pipeline's behaviour
is changed on purpose.
