# Understanding Metagenomic Screening

This page explains what metagenomic screening is and why it's part of this project. It doesn't cover the SNP/coverage side of the pipeline — see the Sample Report guide for that.

## Why this exists

Well under 1% of each sample's raw reads actually align to the *Coffea arabica* reference before any variant is called (see each sample's own report page for its exact read funnel). That's expected for degraded, shipwreck-recovered DNA sequenced without target enrichment — but it raises the question of what the other 99%+ actually is. Metagenomic screening answers that directly, by classifying a subsample of those unmapped reads against a database of known organisms so a sample can be checked for DNA beyond the target reference — contaminants, environmental microbes, and the like.

## How it works, roughly

- It breaks every read into overlapping k-mers (short, fixed-length chunks — 35 base pairs by default).
- Each k-mer is looked up in a precomputed database: essentially a giant hash table mapping k-mers to the lowest common ancestor (LCA) in the NCBI taxonomy tree of every reference genome that contains that k-mer. A k-mer found only in *E. coli* maps to *E. coli*; one shared across all Proteobacteria maps up to that broader group instead.
- A read gets classified to whichever taxon accumulates the most k-mer "votes" along the tree.

## Reading a sample's screening report

Each sample's metagenomic report page shows two passes against different databases, answering two different questions:

- **Environmental / contamination profile** — classified against a broad, prebuilt database (bacteria, archaea, viruses, protozoa, fungi, human). Answers "what is this off-target DNA?"
- **Direct coffee-match check** — classified against a dedicated database built from just *Coffea arabica* and its two diploid parent species. Answers "does any of the off-target pool actually look like coffee?" directly — the broad database alone can't, since it has no plant genomes in it.

A read classifying as *Coffea* isn't automatically real coffee signal, either — Kraken2 can call a match from a single 35bp k-mer, so any such hits get re-aligned to the reference and checked directly (see the "Are those coffee-classified reads real?" section on a sample's report page, where present).
