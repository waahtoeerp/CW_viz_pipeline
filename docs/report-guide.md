# Understanding the Sample Report

This page explains the concepts and terms used on each sample's report page — what the numbers mean and how to read the charts. It doesn't cover how the underlying data was produced; see the Pipeline page for that.

## The basics

**SNP** (single-nucleotide polymorphism) — a single DNA position where the sample's sequence differs from the reference genome. Each SNP is a *candidate*: something the variant caller flagged as a possible real difference, not a confirmed one.

**VCF** (Variant Call Format) — the standard file format variant callers write candidate SNPs into. Every SNP in the report's table and charts comes from the sample's filtered VCF.

**Reference genome** — the sample's DNA is compared against a known, assembled genome for the species (here, *Coffea arabica*). A SNP is a difference relative to this reference, not necessarily an error or a mutation — some are just normal genetic variation.

## Reading a SNP's quality

Every candidate SNP carries a few quality signals, all shown in the SNP table:

- **QUAL** — the variant caller's confidence score that this position is a real variant, not sequencing noise. Higher is better.
- **DP (depth)** — how many sequencing reads actually covered this position. Low depth means the call rests on very few reads, which is much less trustworthy.
- **Genotype (GT)** — the called genotype at this position, e.g. `1/1` (homozygous for the alternate allele) or `0/1` (heterozygous — reference and alternate present). Ancient, degraded samples are usually shallow enough that heterozygous calls are hard to trust.
- **FILTER / Status** — whether the SNP passed every quality check:
  - **PASS** — passed depth, mapping quality, and quality-by-depth thresholds.
  - **flagged: low mapping quality** — the reads at this position didn't align confidently (repetitive or ambiguous region of the genome).
  - **flagged: low depth** — too few reads covered this position to trust the call.

## Why the PASS fraction is often small

These are ancient, degraded coffee-bean DNA samples, not fresh high-coverage sequencing. Depth is often well under 1× across most of the genome — meaning most positions are covered by less than one read on average. At that depth, most candidate SNPs simply don't have enough supporting reads to clear every filter. A low "passing all filters" fraction is expected for this kind of sample, not necessarily a sign something went wrong — it's a reason to weight PASS calls much more heavily than flagged ones when interpreting results.

## Chromosomes, organelles, and scaffolds

The reference genome isn't just one sequence — it's split into different kinds of contigs, and the report treats them differently:

- **Chromosomes** — the main nuclear chromosomes, large sequences (≥5 Mb here). Most analysis focuses on these.
- **Organelles** — the chloroplast and mitochondrial genomes. These are naturally present in many more copies per cell than nuclear DNA, so much higher read depth there is expected and not a quality signal.
- **Scaffolds** — smaller, often incompletely assembled pieces of the genome that didn't get placed onto a full chromosome. Generally less reliable for analysis than chromosomes.

## The charts, section by section

**"Where the SNPs sit across the genome"** — one row per chromosome, drawn to scale by length, with a tick for every candidate SNP. Tick color shows PASS/flagged status. This gives a quick sense of whether SNPs cluster in particular regions or spread evenly.

**"Coverage & depth per contig"** — mean read depth per chromosome/contig, on a log scale (since depth commonly spans orders of magnitude between nuclear and organelle contigs). Faded bars mark contigs with low mean mapping quality (MAPQ < 10) — a sign of repetitive or ambiguous alignment rather than genuinely low coverage.

**Zoom sections** — the report automatically picks two interesting spots to zoom into: the single highest-confidence SNP, and the densest cluster of nearby SNPs on one contig. Each zoom shows local read depth around the SNP(s) and, where available, any NCBI RefSeq gene annotations overlapping that window — so you can see whether an interesting SNP falls inside a known gene.

**"All candidate SNPs"** — the full table behind the charts, filterable by PASS/flagged status.
