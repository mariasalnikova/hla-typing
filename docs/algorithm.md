# The algorithm, step by step

The method is described in section 4.1 of A. A. Minervina's PhD thesis
("Разработка алгоритма анализа данных секвенирования для HLA-типирования").
The original script is located https://github.com/asya-minervina/HLA_old.
This repository is a restructured, documented version of the updated
analysis scripts.

This document describes what the code does, in the order it does it, with
pointers to the functions in `R/`. Names such as `safety1`, `safety2` are the
historical names of the intermediate objects and are kept in the output
files.

## Background: the library design

Each HLA locus is amplified as short overlapping cDNA amplicons (< 500 bp) so
that a 2x250 MiSeq read pair covers one amplicon: read 1 covers the first
~230 nt from one primer, read 2 the last ~230 nt from the other primer, and
the two reads overlap in the middle. The primers are designed so that the
same amplicon can be sequenced in either orientation. Consequently every
amplicon exists in two flavours in the code: the forward one (e.g. `Iamp1`),
where read 1 starts at the "left" primer, and the inverse one (`Iamp1_inv`),
where read 1 starts at the "right" primer. Both are processed identically and
the inverse assemblies are reverse-complemented at the end so that all
sequences are in the same orientation.

Amplicons per locus group (see `amplicon_table()` in `R/amplicons.R`):

| group | alleles | amplicons (each also as `_inv`) |
|---|---|---|
| class I | A, B, C | `Iamp1`, `Iamp2`, `Iamp1alt`, `Iamp2alt` |
| DRB | DRB1, DRB3, DRB4, DRB5 | `IIamp1_Others`, `IIamp2`, `IIamp1_DRBalt`, `IIamp2_DRBalt` |
| DQB | DQB1 | `IIamp1_DQB`, `IIamp2`, `IIamp1_DQBalt`, `IIamp2_DQBalt` |
| DPB | DPB1 | `IIamp1_DPB`, `IIamp2_DPB` |
| DQA | DQA1 | `IIamp_DQA` |
| DPA | DPA1 | `IIamp_DPA` |

`alt` amplicons come from an alternative primer set covering the same
region. `IIamp2` (shared by DRB and DQB) and `IIamp2_DPB` use the same
read-1 primer, so the same reads are analysed twice.

## Stage 1: from reads to trusted amplicon sequences

`assemble_amplicons()` in `R/stage1_assemble.R`. Input: one FASTQ pair.
Output: `list(safety1, safety2, statistics)`.

### 1.1 Reading and length filter (`read_fastq_pair`)

Only the sequence lines are read; base qualities are never used. Read pairs
where either read is <= 200 nt are dropped (`min_read_length`).

### 1.2 Demultiplexing by primer (`demultiplex_reads`)

A read pair belongs to amplicon X when the first `read1_prefix` (20-24)
nucleotides of read 1 match the regular expression `read1_pattern` of X.
Patterns contain degenerate positions (`[GC]`) and alternatives (`A|B`) to
cover primer variants. For `IIamp1_DQB_inv` and `IIamp1_Others_inv` read 1
starts with the same primer, so read 2 is used to tell them apart.
Reads matching no pattern are ignored; a pair can land in several
amplicons.

### 1.3 Clustering identical read 1 (`cluster_identical_reads`)

Both reads are cut to positions 20..250 (`primer_length`..`read_length`):
this removes the primer and the low-quality tail, leaving 231 nt. Read pairs
are grouped by their (trimmed) read-1 sequence. Only sequences seen more than
once are kept, then the `threshold` (100) largest groups with more than two
stored copies.

Quirk preserved from the original: the grouping starts from
`duplicated(read1)`, which drops the first copy of every sequence, so the
stored `readnumber` of a group is the true count minus one, and the "more
than two" rule really means "at least four reads".

### 1.4 Read-2 consensus (`consensus_sequence`, `summarise_clusters`)

Read 2 has lower quality than read 1. For each read-1 group a majority
consensus of the read 2s is built: a position gets the majority base if it
holds > 70 % of the reads, otherwise `N`. Each group becomes one row:
`read1`, `read2cons`, `readnumber`.

### 1.5 Assembling the amplicon (`R/assembly.R`)

Read 1 is joined with the reverse complement of the read-2 consensus.

- `merge_by_overlap` (default) slides the reverse-complemented read 2 along
  read 1 (`find_read_overlap`) and, for every offset, records the longest run
  of identical nucleotides. The offset with the longest run is the start of
  the overlap. The assembly is read 1 followed by the non-overlapping tail of
  read 2; inside the overlap only read 1 is used.
- `merge_fixed_overlap` assumes the overlap length `shift` given in the
  amplicon table (`--merge-method fixed`).
- For DQA, DPA and DPB the table stacks the overlap-search assembly *and*
  fixed-overlap assemblies (`fix_shifts`), because the overlap is hard to
  locate reliably there. The same read-1 group therefore appears 2-3 times
  with different `assembled` sequences; the wrong ones simply never match the
  database.
- `Iamp1alt` reads do not overlap; they are joined with a single `N`
  (`merge_without_overlap`).

Columns added: `Overlap_max_aligned`, `Overlap_start`, `assembled`.

### 1.6 Exact database lookup (`match_assemblies_exact`)

Every assembly is searched as a substring in all allele sequences of
`hla_nuc.fasta` (an `N` in the assembly matches any base). Both the assembly
(`Exact`) and its reverse complement (`Exact_rev`) are searched. The result is
a space-separated list of allele names. This is the slowest step.

### 1.7 Error graph (`rank_parents`, `R/graph_filter.R`)

Within an amplicon, all read-1 sequences are vertices of a graph; two are
connected when their Hamming distance is exactly 1. PCR and sequencing errors
create sequences one substitution away from the true one, so true sequences
are hubs. Per sequence:

- `neighbours_read1` – degree in the graph
- `clusters` – connected component
- `freq` – share of the amplicon's reads in this sequence
- `parents` – rank inside the component by degree (1 = most connected)

### 1.8 Orientation and trusted set

Inverse amplicons are reverse-complemented (`straighten_inverse`: `Exact`
takes the value of `Exact_rev`, `assembled` is reverse-complemented).

`safety1` is the full annotated table per amplicon. `safety2` keeps only
trusted sequences: `readnumber >= min_reads - 1` (i.e. at least 10 stored
reads) and either rank 1 in its component, or rank 2-4 with `freq > 0.05`.

`statistics` records per amplicon the number of reads, of groups, of trusted
sequences, and how many of their reads found an exact database match
(`used`/`non_used`).

## Stage 2: from trusted sequences to allele calls

`type_hla()` in `R/stage2_typing.R`. Input: a stage-1 result. Output:
`list(safety3, safety4, safety5)`. Returns `NULL` if no class I amplicon has a
trusted sequence.

### 2.1 Second filter

Sequences from `safety2` with `freq > 0.01` (`min_freq`) and rank 1-4.

### 2.2 Mapping onto alleles (`map_sequences_to_alleles`, `safety3`)

Every allele sequence is padded with 150 `N` on both ends so that an
amplicon that overhangs a short (partial, e.g. exon-only) database entry can
still match. For each locus group and each of its amplicons, every trusted
sequence is matched exactly (`Biostrings::vcountPattern`, `N` tolerant)
against the alleles of the group. Two allele-by-amplicon tables result:

- `HLA_allele_<group>`: sum of `freq` of the matching sequences
- `HLA_<group>_amps`: indices of the matching sequences as a string, e.g.
  `"0,1,3"` (the leading `0` is the initial value, preserved because it is
  part of the group signatures)

plus `n_reads<group>`: reads per amplicon.

### 2.3 Allele groups and scores (`score_allele_groups`, `safety4`)

Alleles with identical index strings across all amplicons of the group are
indistinguishable with this design; they are merged into one row whose
`Allele` lists all of them, keyed by the signature `amps`
(`"0,1_0,2_0_..."`). Alleles matched by nothing are removed. Then per row:

- `sumreads` = sum over amplicons of `freq x reads in the amplicon`
- `no_amps` = amplicons where the allele was not seen
- `Score` = product over amplicons of `0.95` if seen, `0.05` (`p_miss`) if
  not. With 8 amplicons an allele seen everywhere scores 0.95^8 = 0.66, one
  missing 4 amplicons scores 0.95^4 x 0.05^4 = 5e-6.
- `degree`: a directed "dominance" graph is built between rows; row i points
  to row j when every sequence supporting j also supports i (j's evidence is
  a subset of i's). `degree` is the number of rows that dominate the row.

### 2.4 Calls (`call_alleles`, `safety5`)

Rows are reported when nothing dominates them (`degree == 0`) or when they
were seen in every amplicon (`no_amps == ""`). Scores are then normalised to
sum to 1 within each locus (A, B, C, DQB1, DRB1, DRB3/4/5, DPB1, DQA1, DPA1).
A heterozygous locus therefore shows two rows with ~0.5, a homozygous or
heterozygous locus with non-identified allele shows a single row with ~1;
weakly supported alternatives keep a tiny score.

## Summary table (`combine_typing_results`, `R/report.R`)

All samples' `safety5` tables stacked, with `donor`, `tidyAllele` (unique
3-field names) and `n_noamps`. The original workflow then filtered this table
manually to <= 2 rows per locus and donor and converted it to one row per
donor with `typing_table_to_wide()` (`R/report_wide.R`), which also flags the
DR3/DR4 haplotypes of the type 1 diabetes study; that function contains
data-set specific corrections and is kept only for reference.

## Parameters worth knowing

| parameter | where | default | effect |
|---|---|---|---|
| `read_length`, `primer_length` | stage 1 | 250, 20 | trimming; the fixed overlap lengths in the amplicon table assume 231-nt trimmed reads |
| `threshold` | stage 1 | 100 | number of read-1 groups kept per amplicon |
| consensus threshold | `consensus_sequence` | 0.7 | majority needed for a base in the read-2 consensus |
| `min_reads`, `min_freq_secondary` | stage 1 | 10, 0.05 | trusted-sequence filter |
| `min_freq` | stage 2 | 0.01 | second filter |
| `pad_n` | stage 2 | 150 | N padding of database sequences |
| `p_miss` | stage 2 | 0.05 | missing-amplicon probability in the score |

## Known limitations

- Base qualities are ignored; the read-2 consensus and the 1-mismatch graph
  are the only error correction.
- Only substitutions are modelled (Hamming distance); indel errors form
  separate sequences.
- The exact-match lookup in stage 1 searches every assembly against the
  whole database and dominates the run time.
- In case of only one allele identified it is impossible to determine
  whether it is homozygous or heterozygous locus with non-identified allele.
