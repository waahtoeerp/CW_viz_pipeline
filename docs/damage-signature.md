# Authentic Ancient-DNA Damage Signature

This is the project's headline authentication result: direct molecular evidence that the recovered coffee DNA is genuinely of historical origin, independent of how little of it survived. Source: `the_coffee_wrecks/poster/NOTES.md`.

## The figure

![Ancient-DNA deamination signature, all 4 samples, post-dedup](figures/damage_curves_all_samples.png)

Built from `mapDamage-out/*_mapDamage_dedup/{5pCtoT,3pGtoA}_freq.txt` — **post-dedup**, so PCR copies of one damaged original molecule can't inflate the apparent signal. Red is the 5′ C→T substitution frequency (deamination at the start of each read); blue is the 3′ G→A frequency (deamination at the end). Position 1 is the very end of the read.

## Why this is the authentication result

Ancient DNA degrades in a specific, well-characterized way: cytosine bases near the ends of fragments spontaneously deaminate to uracil over time, which read as thymine (C→T) during sequencing — concentrated right at the fragment terminus and decaying back to a flat baseline within the first ~8–10bp. Modern contaminant DNA doesn't show this pattern. All 4 samples show it clearly:

| Sample | Position-1 C→T | Baseline (pos 11–15 avg) | Ratio |
|---|---:|---:|---:|
| R0002 | 3.16% | 0.84% | 3.8x |
| R0003 | 4.24% | 0.61% | 7.0x |
| R0004 | 4.92% | 1.06% | 4.6x |
| R0005 | 5.44% | 1.10% | 5.0x |

That decaying spike — elevated right at the read end, dropping to baseline within ~10bp — is the signature itself. It's present in every sample, which is the point: however little endogenous DNA survived (well under 1% of raw reads in every sample — see each sample's own Sample Quality Control card for its exact read funnel), what did survive is authentically ancient.

## Why the 5′ and 3′ curves don't mirror each other

The red (5′ C→T) curve is a clear decaying spike in every sample; the blue (3′ G→A) curve is visibly present but much weaker and noisier — worth stating plainly rather than glossing over, since asymmetric 5′/3′ damage strength is itself a real, reportable pattern, not a data-quality problem.

Confirmed directly from the raw frequency values, not just the plot: e.g. R0002's position-1 G→A is 0.0078 vs. C→T's 0.0316 — about 4x lower — and stays flat/noisy across all 25 positions instead of showing a clear decaying peak the way C→T does. Same pattern in all 4 samples.

In a *double-stranded* aDNA library, both curves are normally expected to mirror each other — both strands carry terminal deamination and get sequenced from both directions (the classic symmetric picture; see Briggs et al. 2007). A strong 5′ C→T with a much weaker, flat 3′ G→A is instead the known signature of **single-stranded library preparation** (or a double-stranded prep with an end-repair/fill-in step), where the overhang chemistry isn't symmetric the way double-stranded ligation is — the most likely explanation here. Partial UDG/USER treatment could also flatten the signal, but would be expected to truncate *both* ends, which isn't what's observed (C→T stays strong) — so it's a less likely explanation.

**Open item:** the library prep kit/protocol isn't documented anywhere in this pipeline's records. Worth confirming with the sequencing facility before stating this more definitively — "consistent with [specific protocol]" rather than left as an inferred best guess.
