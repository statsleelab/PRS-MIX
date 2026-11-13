# Paired three-layout pilot (prespecified, 2026-11-05)

Written before any pilot output existed. Exploratory: its results decide what
to test next and do not by themselves change the manuscript's main layout.

## Question

In the 250 kb sweep, every method loses about 15% of its r² against PRS-CS's EUR
cut points (`results/null_sim_layout_comparison.tsv`, EUR+ASN 70:30), and in
the EUR+AFR and ASN+AFR settings minority-only LD predicts better than the
generating mixture. There are two competing explanations, and both may hold:

- **Boundary truncation.** Block-diagonal LD discards correlation between
  blocks, and 250 kb tiles cut more of it.
- **Reference mismatch.** Reference-panel LD differs from the covariance the
  GWAS was generated under. Each ancestry's generating covariance is the
  in-sample LD of its discovery individuals: only 205 AFR individuals, against
  5,070 EUR and 2,701 ASN.

## Design

- **Settings.** Two sweep settings, EUR+AFR 50:50 and EUR+ASN 50:50, replicates
  1–20. Their inputs are frozen copies with sha256 manifests: the split,
  genotypes, null and non-null z, sumstats and true effects, the existing 250 kb
  outputs, and the original panels.
- **Layouts.** Fixed 250 kb (rebuilt), fixed 1 Mb (`fixed1000kb`), and PRS-CS
  EUR cut points (`eurcuts`). Every panel is restricted to the frozen analysis
  SNPs, and every arm runs on one shared SNP set, so the layouts differ only in
  their blocks.
- **Panel sources, per layout.**
  - *Reference* panels, as in the sweep. Arms: both single panels, and PRS-MIX
    at w = 0.1–0.9.
  - *Discovery* panels, built from the GWAS's own discovery individuals, which
    differ from the generating covariance only by the panel estimator. Arms:
    both single panels, and PRS-MIX at w₀. This is a diagnostic, not a usable
    method.
- **Seed controls.** PRS-MIX at w₀ and first-ancestry PRS-CS rerun with seed
  2000 + rep.
- **Longer-chain check.** Replicates 1–5, both single panels and PRS-MIX at w₀,
  with 4,000 iterations, 2,000 burn-in and thinning 5. The defaults are 1,000,
  500 and 5.
- **Sampler settings.** Passed explicitly to every arm: a = 1, b = 0.5, φ
  learned, main seed 1000 + rep. Wall time and peak memory are recorded per arm.
- **Sampler-free diagnostics** (`16_layout_pilot_diagnostics.R`):
  - the share of the generating covariance's squared correlation, within 3 Mb,
    that falls between blocks (raw, and with 1/n noise removed);
  - the within-block distance between each panel and the generating covariance;
  - the null-z LD-fit loss for each layout and source, with the block-diagonal
    generating covariance as the oracle.

## Readout and what each explanation predicts

The primary measures are paired r² and MSE differences against the true genetic
value, by target ancestry, with 95% intervals over 20 replicates.

| Observation | Boundary truncation | Reference mismatch |
|---|---|---|
| 1 Mb and EUR cuts vs 250 kb, reference panels | r² rises for all methods, most in EUR/ASN targets, where LD reaches further | Little change |
| Cross-block share of the generating covariance | Large at 250 kb, smaller at 1 Mb and EUR cuts | Not informative |
| Discovery vs reference panels, same layout | Similar gap at every layout | Discovery panels clearly better, most for AFR |
| Best weight with discovery panels | — | Near w₀ (the AFR-only advantage disappears) |
| Panel-vs-generating distance | — | Much larger for AFR reference than for EUR |

Attribution is partial. 1 Mb tiles and EUR cuts differ in both block size and
boundary placement, so a gain at EUR cuts beyond 1 Mb is only suggestive of
placement. Separating the two needs a size-matched, shifted-boundary control.

The longer-chain check passes if the long-minus-default difference in r² is
within the default seed-rerun spread for every layout. If it is not, the
sampler has not converged at the defaults, and that confounds every other
comparison.

## Deviation (2026-11-09, after results)

The longer-chain check failed for EUR+ASN 50:50. In all 5 replicates, every
arm and every layout, 4,000 iterations gave higher r² than the defaults, by
+0.005 to +0.034. The per-replicate seed-rerun SD is 0.009–0.011. The learned φ
also fell, from 0.00043 to 0.00026 at 250 kb. EUR+AFR 50:50 passed: the
differences were within ±0.005, against a seed SD of 0.0025, and φ was
unchanged.

Added, without changing any other arm: replicates 1–5 of EUR+ASN, both single
panels and PRS-MIX at w₀, at 8,000 iterations, 4,000 burn-in, thinning 5 and
seed 1000 + rep (`MODE=xlong`, arms `xlong_*`). If 4,000 and 8,000 agree within
the seed spread, 4,000 is taken as converged. The default-vs-4,000 shift in the
paired differences then decides whether the EUR+ASN comparisons need rerunning
at longer chains.

## Running

Job scripts are in `cluster/jobs/layout_pilot/`, which is not in git:
`bash submit.sh prep <settings>` freezes the inputs and builds the panels, SNPs
and diagnostics. After the prep output is checked, `bash submit.sh prscs
<settings>` runs the arrays and metrics.

Outputs go to `/shared/leed13_shared/grace/make_ld/outputs/layout_pilot/<setting>/`.
