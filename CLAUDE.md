# PRS-MIX

## What this is

PRS-MIX is a focused fork of the [gauss](https://github.com/statsleelab/gauss) R
package (GWAS summary-statistics analysis using a 32,953-sample, 29-population
reference panel), scoped around one specific goal: building **multi-population,
mixed LD reference panels for PRS-CS**, and verifying they reproduce PRS-CS's own
per-population output.

It carries the full gauss package scaffolding (R/, src/, tests/, etc. — see
below for why) plus the PRS-CS mixing work that isn't upstream anywhere else.
This repo split off from the `gauss3` working copy at
`/Users/gkoo/Documents/Research/gauss3` (origin `gracko95/gauss3`) on
2026-09-17; that repo is the base package's active development line — treat
divergence between the two as expected, not as a merge target, unless told
otherwise.

## Why the whole package is here, not just the mixing scripts

The mixing entry point, `generate_prscs_ld_panels()` (`R/prscs_ld_panel.R`),
calls into the package's compiled reference-panel reader
(`src/gauss.cpp`, `ReadReferenceIndex` / `.read_gauss_index`) to pull genotypes
out of the 33KG reference panel. It cannot run standalone — hence the full
`DESCRIPTION`/`NAMESPACE`/`src/`/`configure` scaffolding is included so this
repo is buildable and testable on its own (`devtools::load_all()`,
`devtools::test()`).

## Key components

- **`R/prscs_ld_panel.R`** — `write_prscs_ld_chr()` / `generate_prscs_ld_panels()`.
  Partitions a chromosome into contiguous, non-overlapping, population-agnostic
  blocks (boundaries come from the full SNP `bp` range in the reference index,
  not from any one population's QC'd SNP set) and writes per-block signed LD
  matrices to HDF5 in the exact layout PRS-CS's `parse_ldblk`/`parse_ldblk_mix`
  expect (`blk_1..blk_N`, each with `ldblk` + `snplist` datasets). Shared block
  boundaries across populations are what let `parse_ldblk_mix` align two
  panels by block index — see the file's roxygen docs for the invariant in
  detail.
- **`inst/prscs-mix/`** — vendored, modified copies of PRS-CS's own scripts
  (`PRScs_mix.py`, `v2parse_genet_mix.py`) that compute posterior SNP effects
  from a weighted blend of multiple reference LD panels (e.g. EUR + ASN)
  instead of a single population. Not part of upstream PRS-CS. Must be run
  alongside an unmodified [PRS-CS](https://github.com/getian107/PRScs)
  checkout so `mcmc_gtb.py`/`gigrnd.py` are importable. See
  `inst/prscs-mix/README.md` for the SNP-intersection bug these scripts fix
  relative to earlier, undocumented versions (phantom zeroed SNPs biasing the
  sampler's shared noise-scale parameter when panels had mismatched SNP
  universes).
- **`inst/scripts/`** — worked examples and validation:
  - `generate_prscs_ld_chr_example.R` — example driver for
    `generate_prscs_ld_panels()`.
  - `verify_panel_vs_prscs_reference.R` — checks a generated panel against
    PRS-CS's published reference LD panel for a population (chr22 EUR result
    checked into `inst/scripts/results/`).
  - `compute_pop_ld_window_example.R` / `compute_superpop_ld_window_example.R`
    — single-population vs. superpopulation LD window examples.
  - `compare_maf_cutoff_snp_overlap.R` — compares SNP overlap under different
    MAF cutoffs across panels being mixed.
- **`tests/testthat/test-prscs-ld-panel.R`** — unit tests for the panel
  generator (block boundaries, HDF5 layout, argument validation).

## Recent history worth knowing (carried over from gauss3)

- `ReadReferenceIndex` caches the parsed reference index to avoid repeated
  file reads; a subsequent fix closed a connection leak in
  `.read_gauss_index` and lowered the default `block_size`.
- `PRScs_mix.py` was patched so `main()` intersects every `--ref_dir` panel's
  SNP universe up front, fixing silently-zeroed SNPs from smaller panels
  (e.g. ASN) that previously biased the whole sampler.

## Cluster

Job submission / environment details for Miami University's RedHawk cluster
are in [`cluster/redhawk.yml`](cluster/redhawk.yml). Read that file before
generating or editing any SLURM submission scripts — it's the source of truth
for account/partition/conda-env bindings, not something to re-derive or
hardcode elsewhere.

## Working conventions

- This is an R package (Rcpp/RcppEigen backend, `hdf5r` for LD panel I/O).
  Standard `devtools`/`testthat` workflow applies.
- Treat `inst/prscs-mix/*.py` as vendored + patched upstream code: keep
  changes minimal and documented (see the README's bug writeup as the model
  for how to describe a deviation from upstream PRS-CS).
