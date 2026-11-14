# Curated paper evidence

Results selected for the PRS-MIX manuscript
([statsleelab/PRSMIX_manuscript](https://github.com/statsleelab/PRSMIX_manuscript)).
Everything here is frozen and committed. The manuscript imports from this
folder only, through its explicit import step. It never reads working outputs
(`inst/scripts/results/`), scratch directories or cluster outputs.

One folder per study:

```text
paper/<study>/
  inputs/          frozen copies of the metric tables the assets are computed from
  *.tex            generated LaTeX table bodies (no caption; the manuscript adds it)
  tables.html      browser preview of the same tables (open locally; not imported)
  tables.tsv       every number in the tables, unrounded, with n per cell
  manifest.json    program, command, commits, settings, input and output sha256
```

Rules:

- Generate assets with a committed program, from committed inputs. Commit the
  program first, run it, then commit the outputs, so `manifest.json` records
  the exact program commit (`asset_commit`) and the commit of the inputs
  (`experiment_commit`).
- Never edit generated files by hand. Change the program or inputs and
  regenerate.
- Keep panels, genotypes, posterior draws and logs out of this folder.

## Studies

| Study | Program | Contents |
| --- | --- | --- |
| `null_sim_eur_asn_250kb` | `inst/scripts/null_sim/12_paper_tables.R` | chr22, EUR+ASN 70:30 GWAS, 250 kb blocks, 100 null and 100 non-null replicates, phi auto: null and non-null tables over the full weight sweep, paired differences against the PRS-CS seed rerun, and settings (LaTeX and HTML) |
| `sim250_sweep_250kb` | `inst/scripts/null_sim/18_sweep_paper_tables.R`; `figures/` by `19_sweep_paper_figures.R` | chr22, 250 kb blocks, phi auto, 4 ancestry pairs x 3 GWAS compositions (90:10, 70:30, 50:50), 100 replicates per setting: accuracy (r^2, MSE), top-decile overlap against the seed rerun, multi-cutoff overlap and settings tables; `figures/` holds the EUR-pair top-10% overlap figure with its own manifest. Replaces the earlier 20-replicate placeholder study. |
