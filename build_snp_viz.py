#!/usr/bin/env python3
"""
Build the SNP / coverage / depth visualization HTML page from the small
extract files produced by extract_snp_viz_data.sh:

  viz-data/<SAMPLE>_coverage_by_contig.txt
  viz-data/<SAMPLE>_snps.tsv
  viz-data/<SAMPLE>_snp_local_depth.tsv

Re-derives the "interesting" loci fresh each run instead of hardcoding them:
  - zoom 1: the highest-confidence SNP (PASS if one exists, else best by depth/qual)
  - zoom 2: the densest cluster of nearby SNPs on any single contig

Optionally annotates both zoom windows with NCBI RefSeq gene features for the
reference accession (network access to eutils.ncbi.nlm.nih.gov; skipped with
a warning if unreachable).

Usage:
  ./build_snp_viz.py [--viz-data-dir viz-data] [--out viz-data/SAMPLE_snp_viz.html] [--no-network]
"""
import argparse
import glob
import json
import os
import re
import sys
import urllib.request
import urllib.error

KNOWN_SAMPLES = ["CT600-007R0002", "CT600-007R0003", "CT600-007R0004", "CT600-007R0005"]

CHROMOSOME_MIN_LEN = 5_000_000  # NC_ contigs at/above this are treated as nuclear chromosomes
CLUSTER_WINDOW_BP = 3000        # SNPs within this distance of a neighbor join the same cluster
ZOOM_PADDING_BP = 1500          # extra context on each side of a zoom window
GENE_FETCH_PADDING_BP = 2000
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def find_samples(viz_dir):
    snp_files = sorted(glob.glob(os.path.join(viz_dir, "*_snps.tsv")))
    if not snp_files:
        sys.exit(f"No *_snps.tsv found in {viz_dir} — run extract_snp_viz_data.sh first "
                  f"and copy its output here.")
    return [os.path.basename(f)[: -len("_snps.tsv")] for f in snp_files]


def parse_contigs(path):
    contigs = []
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq = parts
            contigs.append({
                "name": rname,
                "length": int(endpos),
                "numreads": int(numreads),
                "covbases": int(covbases),
                "coverage_pct": float(coverage),
                "meandepth": float(meandepth),
                "meanmapq": float(meanmapq),
            })
    return contigs


def parse_snps(path):
    snps = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 8:
                continue
            chrom, pos, ref, alt, qual, filt, dp, gt = parts
            snps.append({
                "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                "qual": float(qual), "filter": filt, "dp": int(dp), "gt": gt,
            })
    return snps


def parse_local_depth(path):
    """Returns dict: (chrom) -> list of [pos, depth], deduplicated by position."""
    by_chrom = {}
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            chrom, snp_pos, bin_start, mean_depth = parts
            by_chrom.setdefault(chrom, {})[int(bin_start)] = float(mean_depth)
    return {chrom: sorted([p, d] for p, d in bins.items()) for chrom, bins in by_chrom.items()}


# C. arabica reference GC content, computed once directly from the reference
# fasta (GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna) on Roihu:
# 224,740,116 G + 224,548,132 C over 1,198,235,291 ACGT bases = 37.496%.
# Fixed per reference, not per sample -- recompute only if the reference changes.
REFERENCE_GC = 0.37496

BWA_ALN_PARAMS = "bwa aln -l 16500 -n 0.01"


def parse_signal_funnel(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def parse_gc_content(path, sample):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        all_samples = json.load(f)
    return all_samples.get(sample)


def classify_contigs(contigs):
    chromosomes, organelles, scaffolds = [], [], []
    for c in contigs:
        if c["name"].startswith("NC_") and c["length"] >= CHROMOSOME_MIN_LEN:
            chromosomes.append(c)
        elif c["name"].startswith("NC_"):
            organelles.append(c)
        else:
            scaffolds.append(c)
    chromosomes.sort(key=lambda c: c["name"])
    for i, c in enumerate(chromosomes):
        c["label"] = f"chr{i+1:02d}"
    for c in organelles:
        c["label"] = "organelle (" + c["name"] + ")"
    return chromosomes, organelles, scaffolds


def pick_best_snp(snps):
    passing = [s for s in snps if s["filter"] == "PASS"]
    if passing:
        return max(passing, key=lambda s: s["qual"]), True
    return max(snps, key=lambda s: (s["dp"], s["qual"])), False


def find_densest_cluster(snps):
    """Chain-cluster SNPs within CLUSTER_WINDOW_BP of a neighbor on the same contig;
    return the members of the largest cluster (ties broken by first found)."""
    by_chrom = {}
    for s in snps:
        by_chrom.setdefault(s["chrom"], []).append(s)
    best_cluster = []
    for chrom, group in by_chrom.items():
        group.sort(key=lambda s: s["pos"])
        current = [group[0]]
        for s in group[1:]:
            if s["pos"] - current[-1]["pos"] <= CLUSTER_WINDOW_BP:
                current.append(s)
            else:
                if len(current) > len(best_cluster):
                    best_cluster = current
                current = [s]
        if len(current) > len(best_cluster):
            best_cluster = current
    return best_cluster


def build_zoom(chrom, center_pos, curve_by_chrom, contig_len, snps, extra_span=0):
    curve = curve_by_chrom.get(chrom, [])
    lo = max(1, center_pos - ZOOM_PADDING_BP - extra_span)
    hi = min(contig_len, center_pos + ZOOM_PADDING_BP + extra_span)
    window_curve = [p for p in curve if lo <= p[0] <= hi]
    return {"chrom": chrom, "window": [lo, hi], "curve": window_curve}


def fetch_gene_annotations(chrom, window):
    lo = max(1, window[0] - GENE_FETCH_PADDING_BP)
    hi = window[1] + GENE_FETCH_PADDING_BP
    url = (f"{NCBI_EUTILS}?db=nuccore&id={chrom}&rettype=gb&retmode=text"
           f"&seq_start={lo}&seq_stop={hi}")
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        print(f"  (gene annotation fetch skipped for {chrom}: {e})", file=sys.stderr)
        return []

    m = re.search(r"^FEATURES.*?\n(.*?)^ORIGIN", text, re.S | re.M)
    if not m:
        return []
    body = m.group(1)

    feature_re = re.compile(r"^ {5}(\S+)\s+(.*)$")
    genes, current = [], None
    KEEP = {"CDS", "tRNA", "rRNA", "mRNA", "gene"}
    for line in body.splitlines():
        fm = feature_re.match(line)
        if fm:
            if current:
                genes.append(current)
            key, loc = fm.group(1), fm.group(2)
            current = {"type": key, "loc": loc, "quals": ""} if key in KEEP else None
        elif current is not None:
            current["quals"] += " " + line.strip()
    if current:
        genes.append(current)

    out = []
    seen_spans = set()
    type_priority = {"CDS": 0, "tRNA": 1, "rRNA": 1, "mRNA": 2, "gene": 3}
    for g in genes:
        nums = [int(n) for n in re.findall(r"\d+", g["loc"])]
        if len(nums) < 2:
            continue
        start = lo + min(nums) - 1
        end = lo + max(nums) - 1
        span = (start, end)
        gene_m = re.search(r'/gene="([^"]*)"', g["quals"])
        product_m = re.search(r'/product="([^"]*)"', g["quals"])
        name = gene_m.group(1) if gene_m else (product_m.group(1) if product_m else g["type"])
        out.append({
            "name": name, "start": start, "end": end, "type": g["type"],
            "product": product_m.group(1) if product_m else "",
            "priority": type_priority.get(g["type"], 9),
        })
    # collapse duplicate spans (gene/mRNA/CDS often share the same range), keep most specific
    out.sort(key=lambda g: (g["start"], g["end"], g["priority"]))
    deduped = []
    for g in out:
        if deduped and deduped[-1]["start"] == g["start"] and deduped[-1]["end"] == g["end"]:
            continue
        deduped.append(g)
    return deduped


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--viz-data-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "viz-data"))
    ap.add_argument("--sample", default=None, help="build only this sample instead of every sample found")
    ap.add_argument("--out", default=None, help="output path (only valid with --sample)")
    ap.add_argument("--no-network", action="store_true", help="skip fetching NCBI gene annotations")
    args = ap.parse_args()

    viz_dir = args.viz_data_dir
    if args.out and not args.sample:
        sys.exit("--out only makes sense together with --sample")

    samples = [args.sample] if args.sample else find_samples(viz_dir)
    print(f"Samples: {', '.join(samples)}")
    for sample in samples:
        build_sample(sample, viz_dir, args.no_network, args.out)


def build_sample(sample, viz_dir, no_network, out_override=None):
    print(f"--- Sample: {sample} ---")

    contigs = parse_contigs(os.path.join(viz_dir, f"{sample}_coverage_by_contig.txt"))
    snps = parse_snps(os.path.join(viz_dir, f"{sample}_snps.tsv"))
    curve_by_chrom = parse_local_depth(os.path.join(viz_dir, f"{sample}_snp_local_depth.tsv"))
    length_by_chrom = {c["name"]: c["length"] for c in contigs}

    chromosomes, organelles, scaffolds = classify_contigs(contigs)
    print(f"{len(chromosomes)} chromosome(s), {len(organelles)} organelle contig(s), "
          f"{len(scaffolds)} scaffold(s); {len(snps)} candidate SNPs")

    best_snp, best_is_pass = pick_best_snp(snps)
    cluster = find_densest_cluster(snps)

    zoom1 = build_zoom(best_snp["chrom"], best_snp["pos"], curve_by_chrom,
                        length_by_chrom.get(best_snp["chrom"], best_snp["pos"] + ZOOM_PADDING_BP), snps)
    zoom1["title"] = ("Zoom — the highest-confidence SNP"
                       if best_is_pass else
                       "Zoom — the best-supported SNP (none passed every filter)")
    zoom1["subject_pos"] = best_snp["pos"]

    zoom2 = None
    if len(cluster) >= 2:
        span = cluster[-1]["pos"] - cluster[0]["pos"]
        center = (cluster[0]["pos"] + cluster[-1]["pos"]) // 2
        zoom2 = build_zoom(cluster[0]["chrom"], center, curve_by_chrom,
                            length_by_chrom.get(cluster[0]["chrom"], center + ZOOM_PADDING_BP),
                            snps, extra_span=span // 2)
        zoom2["title"] = f"Zoom — densest SNP cluster ({len(cluster)} SNPs on {cluster[0]['chrom']})"

    if not no_network:
        print("Fetching NCBI gene annotations for zoom windows...")
        zoom1["genes"] = fetch_gene_annotations(zoom1["chrom"], zoom1["window"])
        if zoom2:
            zoom2["genes"] = fetch_gene_annotations(zoom2["chrom"], zoom2["window"])
    else:
        zoom1["genes"] = []
        if zoom2:
            zoom2["genes"] = []

    def gene_summary(genes):
        if not genes:
            return "No annotated features found in this window."
        names = ", ".join(sorted({g["name"] for g in genes})[:6])
        return f"Overlapping annotated feature(s): {names}."

    zoom1["description"] = (
        f"{best_snp['chrom']}:{best_snp['pos']:,} ({best_snp['ref']}→{best_snp['alt']}, "
        f"depth {best_snp['dp']}×, {best_snp['filter']}). " + gene_summary(zoom1["genes"])
    )
    if zoom2:
        zoom2["description"] = (
            f"{cluster[0]['chrom']}:{cluster[0]['pos']:,}–{cluster[-1]['pos']:,}, "
            f"{len(cluster)} SNPs. " + gene_summary(zoom2["genes"])
        )

    n_pass = sum(1 for s in snps if s["filter"] == "PASS")
    main_depth = (sum(c["meandepth"] for c in chromosomes) / len(chromosomes)) if chromosomes else 0.0

    signal_funnel = parse_signal_funnel(os.path.join(viz_dir, f"{sample}_signal_funnel.json"))
    gc = parse_gc_content(os.path.join(viz_dir, "read_gc_content.json"), sample)

    metagenomic_dir = os.path.join(os.path.dirname(os.path.abspath(viz_dir)), "metagenomic-data")

    nav_samples = []
    for s in KNOWN_SAMPLES:
        nav_samples.append({
            "id": s,
            "current": s == sample,
            "types": {
                "snp": {
                    "file": f"{s}_snp_viz.html",
                    "ready": os.path.exists(os.path.join(viz_dir, f"{s}_snps.tsv")),
                },
                "metagenomic": {
                    "file": f"../metagenomic-data/{s}_metagenomic_viz.html",
                    "ready": os.path.exists(os.path.join(metagenomic_dir, f"{s}_kraken2_report.txt")),
                },
            },
        })

    nav_docs = [
        {"label": "Report Guide", "file": "../report-guide.html"},
        {"label": "Pipeline", "file": "../pipeline.html"},
        {"label": "File Chart", "file": "../file-chart.html"},
    ]

    data = {
        "sample": sample,
        "view_type": "snp",
        "nav_samples": nav_samples,
        "nav_docs": nav_docs,
        "main_chr": chromosomes,
        "organelles": organelles,
        "scaffolds": scaffolds,
        "all_contigs": contigs,
        "snps": snps,
        "zoom1": zoom1,
        "zoom2": zoom2,
        "stats": {
            "n_chromosomes": len(chromosomes),
            "n_organelles": len(organelles),
            "n_snps": len(snps),
            "n_pass": n_pass,
            "mean_depth_main": main_depth,
        },
        "qc": {
            "signal_funnel": signal_funnel,
            "gc_sample": gc["gc_after"] if gc else None,
            "reference_gc": REFERENCE_GC,
            "bwa_aln_params": BWA_ALN_PARAMS,
            "n_filtered_variants": len(snps),
        },
    }

    out_path = out_override or os.path.join(viz_dir, f"{sample}_snp_viz.html")
    html = render_html(data)
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Wrote {out_path} ({len(html):,} bytes) — open it directly in a browser.")


def render_html(data):
    data_json = json.dumps(data)
    return TEMPLATE.replace("__DATA_JSON__", data_json)


TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>SNP &amp; coverage overview</title>
<style>
.viz-root {
  color-scheme: light;
  --font-header:    'Times New Roman', Times, serif;
  --font-body:      Arial, Helvetica, sans-serif;
  --accent:         #fed95e;
  --surface-1:      #ffffff;
  --page:           #f2e9dc;
  --text-primary:   #211d16;
  --text-secondary: #6b6459;
  --text-muted:     #8a8073;
  --gridline:       #e4dbcb;
  --baseline:       #cdc3b0;
  --border:         rgba(33,29,22,0.12);
  --series-1:       #2a78d6;
  --series-1-wash:  rgba(42,120,214,0.10);
  --good:           #0ca30c;
  --warning:        #fab219;
  --critical:       #d03b3b;
}
.viz-root { font-family: var(--font-body); background: var(--page);
  color: var(--text-primary); padding: 32px 20px 80px; min-height: 100vh; box-sizing: border-box; }
.wrap { max-width: 1080px; margin: 0 auto; }
h1 { font-family: var(--font-header); font-weight: 400; font-size: 26px; margin: 0 0 4px; }
.subtitle { color: var(--text-secondary); font-size: 14px; margin: 0 0 24px; }
.card { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px;
  padding: 20px 22px; margin-bottom: 20px; overflow-x: auto; }
.card h2 { font-family: var(--font-header); font-weight: 400; font-size: 18px; margin: 0 0 4px; }
.card .desc { color: var(--text-secondary); font-size: 13px; margin: 0 0 16px; max-width: 680px; }
.stat-row { display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px; margin-bottom: 20px; }
.stat-tile { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px; padding: 14px 16px; }
.stat-tile .label { color: var(--text-secondary); font-size: 12px; margin-bottom: 6px; }
.stat-tile .value { font-size: 24px; font-weight: 600; line-height: 1.1; }
.stat-tile .value.warn { color: var(--warning); }
.stat-tile .value.good { color: var(--good); }
.banner { background: var(--surface-1); border: 1px solid var(--border); border-left: 4px solid var(--warning);
  border-radius: 10px; padding: 14px 18px; margin-bottom: 20px; font-size: 13.5px; line-height: 1.55; }
.legend { display: flex; gap: 18px; flex-wrap: wrap; font-size: 12px; color: var(--text-secondary); margin-bottom: 10px; }
.legend .item { display: flex; align-items: center; gap: 6px; }
.legend .swatch { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
svg { display: block; overflow: visible; font-family: inherit; }
.axis-label { fill: var(--text-muted); font-size: 10px; }
.chrom-label { fill: var(--text-secondary); font-size: 11px; }
.gene-label { fill: var(--text-secondary); font-size: 10.5px; }
.tooltip { position: fixed; pointer-events: none; background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 8px; padding: 8px 10px; font-size: 12px; line-height: 1.5; color: var(--text-primary);
  box-shadow: 0 4px 16px rgba(0,0,0,0.18); z-index: 100; opacity: 0; transition: opacity .1s; max-width: 260px; }
.tooltip .tt-value { font-weight: 600; }
.tooltip .tt-sub { color: var(--text-secondary); }
table.snp-table { width: 100%; border-collapse: collapse; font-size: 12.5px; }
table.snp-table thead th { position: sticky; top: 0; background: var(--surface-1); text-align: left;
  padding: 8px 10px; color: var(--text-secondary); font-weight: 500; border-bottom: 1px solid var(--gridline); }
table.snp-table td { padding: 7px 10px; border-bottom: 1px solid var(--gridline); font-variant-numeric: tabular-nums; }
table.snp-table tbody tr:hover { background: var(--series-1-wash); }
.table-scroll { max-height: 420px; overflow-y: auto; border: 1px solid var(--gridline); border-radius: 8px; }
.pill { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px; font-weight: 500; }
.pill-good { color: var(--good); background: color-mix(in srgb, var(--good) 14%, transparent); }
.pill-warning { color: var(--warning); background: color-mix(in srgb, var(--warning) 18%, transparent); }
.pill-muted { color: var(--text-muted); background: color-mix(in srgb, var(--text-muted) 14%, transparent); }
.filters { display: flex; gap: 8px; margin-bottom: 12px; flex-wrap: wrap; }
.filter-btn { font: inherit; font-size: 12.5px; cursor: pointer; border: 1px solid var(--gridline);
  background: var(--surface-1); color: var(--text-secondary); border-radius: 999px; padding: 6px 14px; }
.filter-btn.active { border-color: var(--accent); color: var(--text-primary); background: color-mix(in srgb, var(--accent) 30%, transparent); }
.accent-rule { height: 4px; width: 60px; background: var(--accent); margin: 0 0 24px; border-radius: 2px; }
.callout-note { font-size: 12px; color: var(--text-secondary); margin-top: 10px; line-height: 1.5; }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11.5px; }
footer.credits { color: var(--text-muted); font-size: 11.5px; margin-top: 8px; line-height: 1.6; }
.site-ribbon { display: flex; align-items: center; gap: 20px; background: #ffffff;
  border-bottom: 1px solid rgba(33,29,22,0.12); padding: 8px 24px; position: sticky; top: 0; z-index: 50; }
.site-ribbon .ribbon-home { display: flex; align-items: center; }
.site-ribbon .ribbon-home img { height: 32px; width: auto; display: block; }
.site-ribbon .ribbon-links { display: flex; align-items: center; gap: 6px; flex-wrap: wrap; font-family: var(--font-body); font-size: 15px; }
.site-ribbon .ribbon-links > a { text-decoration: none; padding: 7px 15px; border-radius: 999px; color: #6b6459; }
.site-ribbon .ribbon-links > a:hover { background: rgba(33,29,22,0.06); }
.site-ribbon .ribbon-divider { width: 1px; align-self: stretch; background: rgba(33,29,22,0.12); margin: 4px 2px; }
.ribbon-sample { position: relative; }
.ribbon-sample-btn { font: inherit; font-size: 15px; cursor: pointer; border: none; background: none;
  padding: 7px 15px; border-radius: 999px; color: #6b6459; }
.ribbon-sample-btn:hover { background: rgba(33,29,22,0.06); }
.ribbon-sample-btn.current { background: #fed95e; color: #4a3900; font-weight: bold; }
.ribbon-menu { position: absolute; top: 100%; left: 0; margin-top: 4px; background: #ffffff;
  border: 1px solid rgba(33,29,22,0.12); border-radius: 10px; box-shadow: 0 8px 24px rgba(33,29,22,0.14);
  padding: 6px; min-width: 200px; z-index: 60; display: none; }
.ribbon-menu.open { display: block; }
.ribbon-menu a { display: block; text-decoration: none; padding: 8px 12px; border-radius: 8px;
  font-size: 14px; color: #211d16; }
.ribbon-menu a:hover { background: rgba(33,29,22,0.06); }
.ribbon-menu a.current { font-weight: bold; color: #6b5300; background: color-mix(in srgb, #fed95e 30%, transparent); }
.ribbon-menu a.pending { pointer-events: none; opacity: 0.5; }
footer.site-footer { text-align: center; margin-top: 32px; padding-top: 24px; border-top: 1px solid var(--border); }
footer.site-footer img { width: 273px; height: auto; opacity: 0.85; }
</style>
</head>
<body>
<nav class="site-ribbon">
  <a class="ribbon-home" href="../index.html" title="Back to front page">
    <img src="../theme/logo-compass.png" alt="The Coffee Wrecks — back to front page">
  </a>
  <div class="ribbon-links" id="ribbon-links"></div>
</nav>
<div class="viz-root">
<div class="wrap">
  <h1>SNP &amp; coverage overview</h1>
  <p class="subtitle mono" id="subtitle"></p>
  <div class="accent-rule"></div>
  <div class="stat-row" id="stat-row"></div>
  <div class="banner" id="banner"></div>

  <div class="card" id="qc-card" style="display:none">
    <h2>Sample QC: where did the reads go?</h2>
    <p class="desc" id="qc-desc"></p>
    <div class="stat-row" id="qc-funnel-row"></div>
    <div class="stat-row" id="qc-gc-row"></div>
    <p class="callout-note" id="qc-params-note"></p>
  </div>

  <div class="card">
    <h2>Where the SNPs sit across the genome</h2>
    <p class="desc">Each row is a chromosome, drawn to scale by length. Every tick is one candidate SNP call, colored by whether it passed quality filtering.</p>
    <div class="legend">
      <span class="item"><span class="swatch" style="background:var(--good)"></span>PASS</span>
      <span class="item"><span class="swatch" style="background:var(--warning)"></span>flagged: low mapping quality</span>
      <span class="item"><span class="swatch" style="background:var(--text-muted)"></span>flagged: low depth</span>
    </div>
    <div id="ideogram"></div>
    <p class="callout-note">Hover a tick for details.</p>
  </div>

  <div class="card">
    <h2>Coverage &amp; depth per contig</h2>
    <p class="desc">Mean read depth per chromosome/contig (log scale — depth commonly spans several orders of magnitude). Bar length and the label both show the actual value.</p>
    <div id="coverage-chart"></div>
    <p class="callout-note">Faded bars = contigs where reads mapped but mean mapping quality was low (MAPQ &lt; 10) — a sign of repetitive/ambiguous alignment rather than confident coverage. Organelle contigs are naturally high-copy, so higher depth there is expected.</p>
  </div>

  <div class="card">
    <h2 id="zoom1-title"></h2>
    <p class="desc" id="zoom1-desc"></p>
    <div id="zoom1"></div>
  </div>

  <div class="card" id="zoom2-card" style="display:none">
    <h2 id="zoom2-title"></h2>
    <p class="desc" id="zoom2-desc"></p>
    <div id="zoom2"></div>
  </div>

  <div class="card">
    <h2>All candidate SNPs</h2>
    <p class="desc">The full call set behind the charts above.</p>
    <div class="filters" id="table-filters"></div>
    <div class="table-scroll">
      <table class="snp-table">
        <thead><tr><th>Chromosome</th><th>Position</th><th>Change</th><th>Qual</th><th>Depth</th><th>Genotype</th><th>Status</th></tr></thead>
        <tbody id="snp-tbody"></tbody>
      </table>
    </div>
    <footer class="credits" id="footer-credits"></footer>
  </div>
  <footer class="site-footer">
    <img src="../theme/logo-compass.png" alt="The Coffee Wrecks Project compass mark">
  </footer>
</div>
<div class="tooltip" id="tooltip"></div>
</div>

<script type="application/json" id="viz-data">__DATA_JSON__</script>
<script>
(function(){
  const data = JSON.parse(document.getElementById('viz-data').textContent);
  const tooltip = document.getElementById('tooltip');
  const SVGNS = 'http://www.w3.org/2000/svg';

  function elx(tag, attrs, parent) {
    const e = document.createElementNS(SVGNS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function fmt(n, d) { return Number(n).toLocaleString(undefined, {maximumFractionDigits: d ?? 2}); }
  function showTip(evt, lines) {
    tooltip.innerHTML = '';
    lines.forEach(function(l, i){
      const row = document.createElement('div');
      row.className = i === 0 ? 'tt-value' : 'tt-sub';
      row.textContent = l;
      tooltip.appendChild(row);
    });
    tooltip.style.opacity = 1;
    positionTip(evt);
  }
  function positionTip(evt) {
    const x = evt.clientX, y = evt.clientY, pad = 14;
    let left = x + pad, top = y + pad;
    if (left + 260 > window.innerWidth) left = x - 260 - pad;
    if (top + 120 > window.innerHeight) top = y - 120 - pad;
    tooltip.style.left = left + 'px';
    tooltip.style.top = top + 'px';
  }
  function hideTip() { tooltip.style.opacity = 0; }
  function bucket(filter) {
    if (filter === 'PASS') return 'good';
    if (filter === 'LowMQ') return 'warning';
    return 'muted';
  }
  function bucketLabel(filter) {
    if (filter === 'PASS') return 'PASS';
    if (filter === 'LowMQ') return 'flagged: low mapping quality';
    if (filter === 'LowDP') return 'flagged: low depth';
    return 'flagged: low depth & mapping quality';
  }

  document.getElementById('subtitle').textContent = 'Sample ' + data.sample;

  (function ribbon(){
    const wrap = document.getElementById('ribbon-links');
    const home = document.createElement('a');
    home.href = '../index.html';
    home.textContent = 'Front page';
    wrap.appendChild(home);

    const REPORT_TYPES = [['snp', 'SNP & Coverage'], ['metagenomic', 'Metagenomic Screening']];

    function closeAllMenus() {
      document.querySelectorAll('.ribbon-menu.open').forEach(function(m){ m.classList.remove('open'); });
    }

    data.nav_samples.forEach(function(s){
      const box = document.createElement('div');
      box.className = 'ribbon-sample';

      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'ribbon-sample-btn' + (s.current ? ' current' : '');
      btn.textContent = s.id;

      const menu = document.createElement('div');
      menu.className = 'ribbon-menu';
      REPORT_TYPES.forEach(function(pair){
        const key = pair[0], label = pair[1];
        const info = s.types[key];
        const a = document.createElement('a');
        a.textContent = label;
        if (s.current && key === data.view_type) {
          a.className = 'current';
          a.href = '#';
        } else if (info.ready) {
          a.href = info.file;
        } else {
          a.className = 'pending';
          a.href = '#';
        }
        menu.appendChild(a);
      });

      btn.addEventListener('click', function(evt){
        evt.stopPropagation();
        const wasOpen = menu.classList.contains('open');
        closeAllMenus();
        if (!wasOpen) menu.classList.add('open');
      });

      box.appendChild(btn);
      box.appendChild(menu);
      wrap.appendChild(box);
    });
    document.addEventListener('click', closeAllMenus);

    const divider = document.createElement('div');
    divider.className = 'ribbon-divider';
    wrap.appendChild(divider);
    data.nav_docs.forEach(function(d){
      const a = document.createElement('a');
      a.textContent = d.label;
      a.href = d.file;
      wrap.appendChild(a);
    });
  })();

  (function qc(){
    const q = data.qc;
    const sf = q && q.signal_funnel;
    if (!sf && q.gc_sample == null) return;
    document.getElementById('qc-card').style.display = '';

    if (sf) {
      const funnelRow = document.getElementById('qc-funnel-row');
      const funnelTiles = [
        {label: 'Raw reads', value: fmt(sf.raw_reads, 0)},
        {label: 'Mapped (endogenous)', value: fmt(sf.mapped_reads, 0) + '  (' + (sf.mapped_pct_of_raw*100).toFixed(3) + '%)', cls: sf.mapped_pct_of_raw < 0.01 ? 'warn' : ''},
        {label: 'Unique after dedup', value: fmt(sf.unique_reads, 0), cls: sf.unique_reads < 10000 ? 'warn' : ''},
        {label: 'Est. library size', value: fmt(sf.estimated_library_size, 0), cls: sf.estimated_library_size < 10000 ? 'warn' : ''},
        {label: 'Filtered SNPs', value: fmt(q.n_filtered_variants, 0)},
      ];
      funnelTiles.forEach(function(t){
        const tile = document.createElement('div');
        tile.className = 'stat-tile';
        const l = document.createElement('div'); l.className = 'label'; l.textContent = t.label;
        const v = document.createElement('div'); v.className = 'value ' + (t.cls||''); v.textContent = t.value;
        tile.appendChild(l); tile.appendChild(v);
        funnelRow.appendChild(tile);
      });
    }

    if (q.gc_sample != null) {
      const gcRow = document.getElementById('qc-gc-row');
      const diff = (q.gc_sample - q.reference_gc) * 100;
      const gcTiles = [
        {label: 'GC content, this sample\\'s reads', value: (q.gc_sample*100).toFixed(1) + '%'},
        {label: 'GC content, C. arabica reference', value: (q.reference_gc*100).toFixed(1) + '%'},
        {label: 'Difference', value: (diff >= 0 ? '+' : '') + diff.toFixed(1) + ' pts', cls: Math.abs(diff) > 10 ? 'warn' : ''},
      ];
      gcTiles.forEach(function(t){
        const tile = document.createElement('div');
        tile.className = 'stat-tile';
        const l = document.createElement('div'); l.className = 'label'; l.textContent = t.label;
        const v = document.createElement('div'); v.className = 'value ' + (t.cls||''); v.textContent = t.value;
        tile.appendChild(l); tile.appendChild(v);
        gcRow.appendChild(tile);
      });
    }

    let desc = 'How many reads made it through each stage, from raw sequencing output to the variants shown below.';
    if (q.gc_sample != null) {
      desc += ' A large GC-content gap between the raw reads and the reference is itself a contamination signal — coffee DNA reads should land close to the reference\\'s own GC content, so a big gap points to non-target DNA making up most of the read pool (confirmed directly for this project by metagenomic screening — see the Pipeline doc).';
    }
    document.getElementById('qc-desc').textContent = desc;
    document.getElementById('qc-params-note').textContent =
      'Alignment: ' + q.bwa_aln_params + ' (bwa aln, not bwa mem — the standard choice for short/damaged aDNA reads; see the Pipeline doc\\'s "bwa aln parameters" section for why).';
  })();

  (function stats(){
    const row = document.getElementById('stat-row');
    const s = data.stats;
    const tiles = [
      {label: 'Chromosomes analyzed', value: s.n_chromosomes + ' + ' + s.n_organelles + ' organelle(s)'},
      {label: 'Mean depth, main chromosomes', value: s.mean_depth_main.toFixed(4) + '×', cls: s.mean_depth_main < 1 ? 'warn' : ''},
      {label: 'Candidate SNPs', value: fmt(s.n_snps, 0)},
      {label: 'Passing all filters', value: s.n_pass + ' / ' + s.n_snps, cls: s.n_pass < s.n_snps * 0.1 ? 'warn' : 'good'},
    ];
    tiles.forEach(function(t){
      const tile = document.createElement('div');
      tile.className = 'stat-tile';
      const l = document.createElement('div'); l.className = 'label'; l.textContent = t.label;
      const v = document.createElement('div'); v.className = 'value ' + (t.cls||''); v.textContent = t.value;
      tile.appendChild(l); tile.appendChild(v);
      row.appendChild(tile);
    });
  })();

  (function banner(){
    const s = data.stats;
    const b = document.getElementById('banner');
    const passFrac = s.n_snps ? s.n_pass / s.n_snps : 0;
    let msg = 'Mean depth across the ' + s.n_chromosomes + ' main chromosome(s) is ' + s.mean_depth_main.toFixed(4) +
      '×. Of ' + s.n_snps + ' candidate SNPs, ' + s.n_pass + ' pass every quality filter (depth, mapping quality, quality-by-depth).';
    if (passFrac < 0.1) {
      msg += ' That is a small fraction — at this depth, most single-read SNP calls cannot be reliably distinguished from sequencing or DNA-damage noise.';
    }
    b.textContent = msg;
  })();

  (function ideogram(){
    const container = document.getElementById('ideogram');
    const chrs = data.main_chr;
    if (!chrs.length) { container.textContent = 'No chromosome-scale contigs found.'; return; }
    const maxLen = Math.max.apply(null, chrs.map(function(c){ return c.length; }));
    const labelW = 90, rightPad = 60, rowH = 20, topPad = 10;
    const width = 900, plotW = width - labelW - rightPad;
    const height = topPad + rowH * chrs.length + 24;
    const svg = elx('svg', {width: '100%', viewBox: '0 0 ' + width + ' ' + height, height: height});
    container.appendChild(svg);
    const snpsByChrom = {};
    data.snps.forEach(function(s){ (snpsByChrom[s.chrom] = snpsByChrom[s.chrom] || []).push(s); });
    chrs.forEach(function(c, i){
      const y = topPad + i * rowH + rowH/2;
      const x2 = labelW + (c.length / maxLen) * plotW;
      elx('text', {x: 0, y: y+4, class: 'chrom-label'}, svg).textContent = c.label;
      const countText = elx('text', {x: labelW - 8, y: y+4, class: 'chrom-label', 'text-anchor': 'end'}, svg);
      countText.textContent = 'n=' + (snpsByChrom[c.name] || []).length;
      countText.setAttribute('fill', 'var(--text-muted)');
      elx('line', {x1: labelW, x2: x2, y1: y, y2: y, stroke: 'var(--baseline)', 'stroke-width': 2, 'stroke-linecap': 'round'}, svg);
      (snpsByChrom[c.name] || []).forEach(function(s){
        const sx = labelW + (s.pos / c.length) * (x2 - labelW);
        const b = bucket(s.filter);
        const line = elx('line', {x1: sx, x2: sx, y1: y-7, y2: y+7, stroke: 'var(--' + b + ')', 'stroke-width': b === 'good' ? 3 : 2}, svg);
        const hit = elx('rect', {x: sx-6, y: y-10, width: 12, height: 20, fill: 'transparent'}, svg);
        function onEnter(evt){
          showTip(evt, [s.chrom + ' : ' + fmt(s.pos, 0), s.ref + ' → ' + s.alt + '   qual ' + s.qual,
            'depth ' + s.dp + '×   ' + s.gt, bucketLabel(s.filter)]);
        }
        [line, hit].forEach(function(n){
          n.addEventListener('pointerenter', onEnter);
          n.addEventListener('pointermove', positionTip);
          n.addEventListener('pointerleave', hideTip);
        });
      });
    });
    const axisY = topPad + rowH * chrs.length + 14;
    [0, 0.25, 0.5, 0.75, 1].forEach(function(f){
      const x = labelW + f * plotW;
      elx('text', {x: x, y: axisY, class: 'axis-label', 'text-anchor': f===1 ? 'end' : (f===0?'start':'middle')}, svg)
        .textContent = fmt(f * maxLen / 1e6, 0) + ' Mb';
    });
  })();

  (function coverage(){
    const container = document.getElementById('coverage-chart');
    const list = data.main_chr.concat(data.organelles).slice().sort(function(a,b){ return b.meandepth - a.meandepth; });
    const labelW = 190, rightPad = 90, rowH = 22, topPad = 10;
    const width = 900, plotW = width - labelW - rightPad;
    const height = topPad + rowH * list.length + 26;
    const svg = elx('svg', {width: '100%', viewBox: '0 0 ' + width + ' ' + height, height: height});
    container.appendChild(svg);
    const values = list.map(function(c){ return Math.log10(Math.max(c.meandepth, 1e-5)); });
    const minV = Math.min.apply(null, values), maxV = Math.max.apply(null, values);
    list.forEach(function(c, i){
      const y = topPad + i * rowH;
      const v = Math.log10(Math.max(c.meandepth, 1e-5));
      const barW = ((v - minV) / (maxV - minV || 1)) * plotW;
      const isAnomaly = c.meanmapq < 10 && c.numreads > 10;
      const isOrganelle = data.organelles.some(function(o){ return o.name === c.name; });
      elx('text', {x: labelW - 8, y: y + rowH/2 + 4, class: 'chrom-label', 'text-anchor': 'end'}, svg).textContent = c.label || c.name;
      elx('rect', {x: labelW, y: y + 3, width: Math.max(barW, 2), height: rowH - 8, rx: 4,
        fill: isOrganelle ? 'var(--warning)' : 'var(--series-1)', opacity: isAnomaly && !isOrganelle ? 0.55 : 1}, svg);
      elx('text', {x: labelW + barW + 8, y: y + rowH/2 + 4, class: 'axis-label'}, svg).textContent = c.meandepth.toFixed(4) + '×';
      const hit = elx('rect', {x: labelW, y: y, width: Math.max(barW,2)+70, height: rowH, fill: 'transparent'}, svg);
      hit.addEventListener('pointerenter', function(evt){
        showTip(evt, [(c.label || c.name) + '  (' + c.name + ')',
          'mean depth ' + c.meandepth.toFixed(5) + '×,  ' + c.coverage_pct.toFixed(3) + '% of bases covered',
          c.numreads + ' reads,  mean MAPQ ' + c.meanmapq.toFixed(1),
          isAnomaly ? 'Low mean MAPQ for its read count — likely multi-mapping/repetitive region' : ''].filter(Boolean));
      });
      hit.addEventListener('pointermove', positionTip);
      hit.addEventListener('pointerleave', hideTip);
    });
  })();

  function renderZoom(containerId, zoom) {
    const container = document.getElementById(containerId);
    if (!zoom || !zoom.curve || !zoom.curve.length) { container.textContent = 'No depth data in this window.'; return; }
    const width = 900;
    const margin = {left: 50, right: 20};
    const plotW = width - margin.left - margin.right;
    const depthH = 70, geneH = 16, geneRowGap = 4, tickH = 22, gap = 18;
    const genes = zoom.genes || [];
    const window_ = zoom.window;
    const x0 = window_[0], x1 = window_[1];
    function X(pos) { return margin.left + ((pos - x0) / (x1 - x0)) * plotW; }
    const geneRows = genes.length;
    const height = depthH + gap + geneRows * (geneH + geneRowGap) + gap + tickH + 30;
    const svg = elx('svg', {width: '100%', viewBox: '0 0 ' + width + ' ' + height, height: height});
    container.appendChild(svg);

    const curve = zoom.curve;
    const maxDepth = Math.max.apply(null, curve.map(function(p){ return p[1]; }), 1);
    function Y(depth) { return depthH - (depth / maxDepth) * (depthH - 6) - 2; }
    let dLine = 'M ' + X(curve[0][0]) + ' ' + Y(curve[0][1]);
    let dArea = 'M ' + X(curve[0][0]) + ' ' + depthH;
    curve.forEach(function(p){ dLine += ' L ' + X(p[0]) + ' ' + Y(p[1]); dArea += ' L ' + X(p[0]) + ' ' + Y(p[1]); });
    dArea += ' L ' + X(curve[curve.length-1][0]) + ' ' + depthH + ' Z';
    elx('path', {d: dArea, fill: 'var(--series-1-wash)', stroke: 'none'}, svg);
    elx('path', {d: dLine, fill: 'none', stroke: 'var(--series-1)', 'stroke-width': 2, 'stroke-linejoin':'round'}, svg);
    elx('text', {x: margin.left, y: 10, class: 'axis-label'}, svg).textContent = 'depth (max ' + maxDepth.toFixed(1) + '×)';
    elx('line', {x1: margin.left, x2: width - margin.right, y1: depthH, y2: depthH, stroke: 'var(--gridline)', 'stroke-width': 1}, svg);

    let gy = depthH + gap;
    genes.forEach(function(g){
      const gx1 = X(Math.max(g.start, x0)), gx2 = X(Math.min(g.end, x1));
      const thick = g.type === 'CDS';
      elx('rect', {x: gx1, y: gy + (thick?0:geneH/2-2), width: Math.max(gx2-gx1,2), height: thick?geneH:4, rx: 3,
        fill: 'var(--series-1)', opacity: thick?0.8:0.45}, svg);
      const labelX = Math.max(gx1, margin.left);
      const label = g.name.length > 46 ? g.name.slice(0,44)+'…' : g.name;
      elx('text', {x: labelX, y: gy + geneH + 9, class: 'gene-label'}, svg).textContent = label + ' (' + g.type + ')';
      gy += geneH + geneRowGap + 12;
    });

    const tickY = gy + 10;
    elx('line', {x1: margin.left, x2: width - margin.right, y1: tickY, y2: tickY, stroke: 'var(--baseline)', 'stroke-width': 1}, svg);
    const snpsHere = data.snps.filter(function(s){ return s.chrom === zoom.chrom && s.pos >= x0 && s.pos <= x1; });
    snpsHere.forEach(function(s){
      const sx = X(s.pos);
      const b = bucket(s.filter);
      const r = b === 'good' ? 6 : 4.5;
      const dot = elx('circle', {cx: sx, cy: tickY, r: r, fill: 'var(--' + b + ')', stroke: 'var(--surface-1)', 'stroke-width': 2}, svg);
      const hit = elx('circle', {cx: sx, cy: tickY, r: 14, fill: 'transparent'}, svg);
      function onEnter(evt){
        showTip(evt, [s.chrom + ' : ' + fmt(s.pos, 0), s.ref + ' → ' + s.alt + '   qual ' + s.qual,
          'depth ' + s.dp + '×   ' + s.gt, bucketLabel(s.filter)]);
      }
      [dot, hit].forEach(function(n){
        n.addEventListener('pointerenter', onEnter);
        n.addEventListener('pointermove', positionTip);
        n.addEventListener('pointerleave', hideTip);
      });
      if (b === 'good') elx('text', {x: sx, y: tickY - 14, class: 'axis-label', 'text-anchor': 'middle'}, svg).textContent = 'PASS';
    });

    const axisY = height - 6;
    [0, 0.5, 1].forEach(function(f){
      const pos = x0 + f * (x1 - x0);
      elx('text', {x: X(pos), y: axisY, class: 'axis-label', 'text-anchor': f===1?'end':(f===0?'start':'middle')}, svg)
        .textContent = fmt(pos, 0);
    });
  }

  document.getElementById('zoom1-title').textContent = data.zoom1.title;
  document.getElementById('zoom1-desc').textContent = data.zoom1.description;
  renderZoom('zoom1', data.zoom1);

  if (data.zoom2) {
    document.getElementById('zoom2-card').style.display = '';
    document.getElementById('zoom2-title').textContent = data.zoom2.title;
    document.getElementById('zoom2-desc').textContent = data.zoom2.description;
    renderZoom('zoom2', data.zoom2);
  }

  document.getElementById('footer-credits').textContent =
    'Gene annotations (if shown) from NCBI RefSeq, fetched at build time. Regenerate this page with build_snp_viz.py after refreshing viz-data/ from Puhti.';

  (function table(){
    const tbody = document.getElementById('snp-tbody');
    const filterBar = document.getElementById('table-filters');
    const filters = [
      {key: 'all', label: 'All (' + data.snps.length + ')'},
      {key: 'PASS', label: 'PASS (' + data.snps.filter(function(s){return s.filter==='PASS';}).length + ')'},
      {key: 'LowMQ', label: 'Low mapping quality only (' + data.snps.filter(function(s){return s.filter==='LowMQ';}).length + ')'},
      {key: 'LowDP', label: 'Low depth (' + data.snps.filter(function(s){return s.filter.indexOf('LowDP')>=0;}).length + ')'},
    ];
    let active = 'all';
    function matches(s) {
      if (active === 'all') return true;
      if (active === 'LowDP') return s.filter.indexOf('LowDP') >= 0;
      return s.filter === active;
    }
    function renderRows() {
      tbody.innerHTML = '';
      data.snps.filter(matches).forEach(function(s){
        const tr = document.createElement('tr');
        const b = bucket(s.filter);
        const pillClass = b === 'good' ? 'pill-good' : (b === 'warning' ? 'pill-warning' : 'pill-muted');
        [s.chrom, fmt(s.pos,0), s.ref + '→' + s.alt, s.qual, s.dp + '×', s.gt].forEach(function(c){
          const td = document.createElement('td'); td.textContent = c; tr.appendChild(td);
        });
        const statusTd = document.createElement('td');
        const pill = document.createElement('span');
        pill.className = 'pill ' + pillClass; pill.textContent = s.filter;
        statusTd.appendChild(pill); tr.appendChild(statusTd);
        tbody.appendChild(tr);
      });
    }
    filters.forEach(function(f){
      const btn = document.createElement('button');
      btn.className = 'filter-btn' + (f.key === active ? ' active' : '');
      btn.textContent = f.label;
      btn.addEventListener('click', function(){
        active = f.key;
        Array.prototype.forEach.call(filterBar.children, function(b2){ b2.classList.remove('active'); });
        btn.classList.add('active');
        renderRows();
      });
      filterBar.appendChild(btn);
    });
    renderRows();
  })();
})();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
