#!/usr/bin/env python3
"""
Build the metagenomic screening visualization HTML page from the Kraken2
report files fetched from CW_metagenomic_screening's Roihu /scratch results:

  metagenomic-data/<SAMPLE>_kraken2_report.txt          (PlusPF-8, broad/contamination)
  metagenomic-data/<SAMPLE>_kraken2_coffee_report.txt    (dedicated 3-genome Coffea DB)
  metagenomic-data/coffee_hit_repeat_regions.json        (shared, see below)

Background: the_coffee_wrecks aligns <1% of reads to the C. arabica reference.
CW_metagenomic_screening classifies a 2,000,000-read subsample of the
unmapped/off-target reads with Kraken2 against two databases to check whether
that >99% is contamination (expected) or missed real coffee signal. See that
repo's LABDIARY.md for the full narrative this page summarizes.

coffee_hit_repeat_regions.json backs the "where are those repeated areas"
table: every Coffea-classified read was re-aligned with bwa mem (soft-clip
tolerant, unlike the pipeline's stricter bwa aln) against the C. arabica
reference, chain-clustered into genomic windows (<=500bp between neighbors),
and the reference sequence at each cluster was pulled directly (samtools
faidx) and checked by eye for recognizable motifs. Ten of the twelve
largest clusters turned out to be unambiguous bacterial 16S rRNA gene
fragments -- e.g. the NC_092313.1:77046430 cluster opens with
AGAGTTTGATCCTGGCTCAG, the universal 27F primer sequence -- duplicated
near-verbatim across several coffee nuclear chromosomes and the
mitochondrial contig (NC_008535.1). That's a more specific, verified
explanation than CW_metagenomic_screening/LABDIARY.md's first guess
(microsatellite repeats): reads that are themselves ordinary
environmental-bacteria 16S fragments pick up a stray "Coffea" k-mer hit
because the coffee reference assembly carries its own near-identical
copies of this highly conserved gene. This table isn't produced by any
script yet -- it was built manually on Roihu (bwa mem + samtools faidx)
and would need re-running by hand for new samples.

Usage:
  ./build_metagenomic_viz.py [--data-dir metagenomic-data] [--sample X]
"""
import argparse
import glob
import json
import os
import sys

KNOWN_SAMPLES = ["CT600-007R0002", "CT600-007R0003", "CT600-007R0004", "CT600-007R0005"]
SUBSAMPLE_N = 2_000_000

# Read-alignment follow-up isn't in any report file -- it's a manual check
# recorded in CW_metagenomic_screening/LABDIARY.md (2026-09-19 entry): every
# Coffea-classified hit was re-aligned to the C. arabica reference and none
# of them actually aligned. Fill in per-sample once that check has been run;
# a sample missing from this dict just doesn't get the follow-up section.
ALIGNMENT_FOLLOWUP = {
    "CT600-007R0003": {"aligned": 0},
    "CT600-007R0004": {"aligned": 0},
}


def load_repeat_regions(data_dir):
    path = os.path.join(data_dir, "coffee_hit_repeat_regions.json")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


def parse_kraken_report(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                continue
            pct, clade, direct, rank, taxid, name = parts
            rows.append({
                "pct": float(pct), "clade": int(clade), "direct": int(direct),
                "rank": rank, "taxid": int(taxid), "name": name.strip(),
            })
    return rows


def find_samples(data_dir):
    files = sorted(glob.glob(os.path.join(data_dir, "*_kraken2_report.txt")))
    return [os.path.basename(f)[: -len("_kraken2_report.txt")] for f in files]


def build_broad_summary(rows):
    unclassified = next(r for r in rows if r["rank"] == "U")
    total = unclassified["clade"] + next(r for r in rows if r["rank"] == "R")["clade"]
    human = next((r for r in rows if r["taxid"] == 9606), None)
    # Species/genus-level calls only -- these are the mutually-exclusive
    # "leaf" taxa (the higher ranks like Bacteria/Pseudomonadota are
    # cumulative parents of these and would double-count if mixed in).
    leaf = sorted(
        (r for r in rows if r["rank"] in ("S", "G") and r["direct"] > 0),
        key=lambda r: -r["direct"],
    )[:12]
    return {
        "total_reads": total,
        "unclassified_reads": unclassified["clade"],
        "unclassified_pct": unclassified["pct"],
        "human_pct": human["pct"] if human else 0.0,
        "human_reads": human["clade"] if human else 0,
        "top_taxa": [
            {"name": r["name"], "rank": r["rank"], "pct": r["pct"], "reads": r["direct"],
             "is_human": r["taxid"] == 9606}
            for r in leaf
        ],
    }


def build_coffee_summary(rows, sample, repeat_regions):
    total = next(r for r in rows if r["rank"] == "U")["clade"] + next(r for r in rows if r["rank"] == "R")["clade"]
    genus = next((r for r in rows if r["rank"] == "G"), None)
    species = sorted(
        (r for r in rows if r["rank"] == "S"),
        key=lambda r: -r["clade"],
    )
    classified = genus["clade"] if genus else 0
    unassigned_to_species = classified - sum(s["clade"] for s in species)
    followup = ALIGNMENT_FOLLOWUP.get(sample)
    if followup is not None:
        followup = dict(followup, repeat_regions=repeat_regions)
    return {
        "total_reads": total,
        "classified_reads": classified,
        "classified_pct": (classified / total * 100) if total else 0.0,
        "species": [{"name": s["name"], "reads": s["clade"], "pct": s["pct"]} for s in species],
        "unassigned_to_species": unassigned_to_species,
        "followup": followup,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "metagenomic-data"))
    ap.add_argument("--sample", default=None)
    args = ap.parse_args()

    samples = [args.sample] if args.sample else find_samples(args.data_dir)
    if not samples:
        sys.exit(f"No *_kraken2_report.txt found in {args.data_dir}")
    print(f"Samples: {', '.join(samples)}")
    for sample in samples:
        build_sample(sample, args.data_dir)


def build_sample(sample, data_dir):
    print(f"--- Sample: {sample} ---")
    broad_rows = parse_kraken_report(os.path.join(data_dir, f"{sample}_kraken2_report.txt"))
    coffee_rows = parse_kraken_report(os.path.join(data_dir, f"{sample}_kraken2_coffee_report.txt"))

    broad = build_broad_summary(broad_rows)
    coffee = build_coffee_summary(coffee_rows, sample, load_repeat_regions(data_dir))

    viz_dir = os.path.join(os.path.dirname(os.path.abspath(data_dir)), "viz-data")

    nav_samples = []
    for s in KNOWN_SAMPLES:
        nav_samples.append({
            "id": s,
            "current": s == sample,
            "types": {
                "snp": {
                    "file": f"../viz-data/{s}_snp_viz.html",
                    "ready": os.path.exists(os.path.join(viz_dir, f"{s}_snps.tsv")),
                },
                "metagenomic": {
                    "file": f"{s}_metagenomic_viz.html",
                    "ready": os.path.exists(os.path.join(data_dir, f"{s}_kraken2_report.txt")),
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
        "view_type": "metagenomic",
        "subsample_n": SUBSAMPLE_N,
        "nav_samples": nav_samples,
        "nav_docs": nav_docs,
        "broad": broad,
        "coffee": coffee,
    }

    out_path = os.path.join(data_dir, f"{sample}_metagenomic_viz.html")
    html = render_html(data)
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Wrote {out_path} ({len(html):,} bytes) — open it directly in a browser.")


def render_html(data):
    return TEMPLATE.replace("__DATA_JSON__", json.dumps(data))


TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Metagenomic screening overview</title>
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
  --series-2:       #0ca36e;
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
.banner { background: var(--surface-1); border: 1px solid var(--border); border-left: 4px solid var(--good);
  border-radius: 10px; padding: 14px 18px; margin-bottom: 20px; font-size: 13.5px; line-height: 1.55; }
.legend { display: flex; gap: 18px; flex-wrap: wrap; font-size: 12px; color: var(--text-secondary); margin-bottom: 10px; }
.legend .item { display: flex; align-items: center; gap: 6px; }
.legend .swatch { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
svg { display: block; overflow: visible; font-family: inherit; }
.axis-label { fill: var(--text-muted); font-size: 10px; }
.chrom-label { fill: var(--text-secondary); font-size: 11px; }
.tooltip { position: fixed; pointer-events: none; background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 8px; padding: 8px 10px; font-size: 12px; line-height: 1.5; color: var(--text-primary);
  box-shadow: 0 4px 16px rgba(0,0,0,0.18); z-index: 100; opacity: 0; transition: opacity .1s; max-width: 280px; }
.tooltip .tt-value { font-weight: 600; }
.tooltip .tt-sub { color: var(--text-secondary); }
.accent-rule { height: 4px; width: 60px; background: var(--accent); margin: 0 0 24px; border-radius: 2px; }
.callout-note { font-size: 12px; color: var(--text-secondary); margin-top: 10px; line-height: 1.5; }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11.5px; }
.followup { background: color-mix(in srgb, var(--good) 8%, var(--surface-1)); border: 1px solid color-mix(in srgb, var(--good) 25%, var(--border));
  border-radius: 10px; padding: 14px 18px; font-size: 13.5px; line-height: 1.6; margin-bottom: 16px; }
.followup .headline { font-weight: bold; margin-bottom: 6px; }
table.region-table { width: 100%; border-collapse: collapse; font-size: 12.5px; }
table.region-table thead th { text-align: left; padding: 8px 10px; color: var(--text-secondary);
  font-weight: 500; border-bottom: 1px solid var(--gridline); }
table.region-table td { padding: 7px 10px; border-bottom: 1px solid var(--gridline); vertical-align: top; }
table.region-table tbody tr:hover { background: var(--series-1-wash); }
table.region-table .region-cell { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11.5px; white-space: nowrap; }
table.region-table .num { font-variant-numeric: tabular-nums; text-align: right; }
.pill { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11px; font-weight: 500; }
.pill-good { color: var(--good); background: color-mix(in srgb, var(--good) 14%, transparent); }
.pill-muted { color: var(--text-muted); background: color-mix(in srgb, var(--text-muted) 14%, transparent); }
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
  <h1>Metagenomic screening overview</h1>
  <p class="subtitle mono" id="subtitle"></p>
  <div class="accent-rule"></div>
  <div class="stat-row" id="stat-row"></div>
  <div class="banner" id="banner"></div>

  <div class="card">
    <h2>Environmental / contamination profile</h2>
    <p class="desc">A 2,000,000-read subsample of the reads that didn't align to the <em>C. arabica</em> reference, classified against Kraken2's PlusPF-8 database (bacteria, archaea, viruses, protozoa, fungi, human). Bars show the species/genus-level taxa Kraken2 could actually pin down — the higher ranks above them in the tree (Bacteria, Pseudomonadota, ...) are cumulative parents of these and are left out here to avoid double-counting.</p>
    <div class="legend">
      <span class="item"><span class="swatch" style="background:var(--series-1)"></span>bacteria</span>
      <span class="item"><span class="swatch" style="background:var(--warning)"></span>human (handling contamination)</span>
    </div>
    <div id="broad-chart"></div>
    <p class="callout-note" id="broad-note"></p>
  </div>

  <div class="card">
    <h2>Direct coffee-match check</h2>
    <p class="desc">The same subsample, classified instead against a dedicated database built from just three genomes: <em>C. arabica</em> (the target) and its two diploid parent species <em>C. canephora</em> and <em>C. eugenioides</em> (arabica is a natural hybrid of the two). This asks directly: does any of the off-target read pool actually look like coffee?</p>
    <div id="coffee-chart"></div>
    <p class="callout-note" id="coffee-note"></p>
  </div>

  <div class="card" id="followup-card" style="display:none">
    <h2>Are those coffee-classified reads real?</h2>
    <p class="desc">A k-mer hit isn't the same as a confirmed match — Kraken2 can classify a read from a single 35bp k-mer match. Every Coffea-classified hit was re-aligned to the <em>C. arabica</em> reference with the pipeline's own bwa parameters to check directly, then re-aligned again with a soft-clip-tolerant aligner (bwa mem) and clustered by genomic position to see exactly where on the reference they're landing.</p>
    <div class="followup" id="followup-body"></div>
    <p class="desc" id="region-table-desc" style="margin-top:18px"></p>
    <div class="table-scroll" style="max-height:420px; overflow-y:auto; border:1px solid var(--gridline); border-radius:8px;">
      <table class="region-table">
        <thead><tr><th>Genomic region</th><th>Reads here (this sample)</th><th>MAPQ 0</th><th>Mean mismatches</th><th>What's actually there</th></tr></thead>
        <tbody id="region-tbody"></tbody>
      </table>
    </div>
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
    if (left + 280 > window.innerWidth) left = x - 280 - pad;
    if (top + 120 > window.innerHeight) top = y - 120 - pad;
    tooltip.style.left = left + 'px';
    tooltip.style.top = top + 'px';
  }
  function hideTip() { tooltip.style.opacity = 0; }

  document.getElementById('subtitle').textContent = 'Sample ' + data.sample + '  ·  ' + fmt(data.subsample_n, 0) + '-read subsample of off-target reads';

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

  (function stats(){
    const row = document.getElementById('stat-row');
    const b = data.broad, c = data.coffee;
    const tiles = [
      {label: 'Unclassified (PlusPF-8)', value: b.unclassified_pct.toFixed(1) + '%'},
      {label: 'Human (handling contamination)', value: b.human_pct.toFixed(2) + '%', cls: b.human_pct < 1 ? 'good' : 'warn'},
      {label: 'Classified as Coffea', value: c.classified_pct.toFixed(3) + '%  (' + fmt(c.classified_reads,0) + ' reads)'},
      {label: 'Of those, aligned to reference', value: c.followup ? (c.followup.aligned + ' / ' + fmt(c.classified_reads,0)) : 'not yet checked', cls: c.followup && c.followup.aligned === 0 ? 'warn' : ''},
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
    const b = data.broad, c = data.coffee;
    const bnr = document.getElementById('banner');
    let msg = fmt(b.unclassified_reads,0) + ' of ' + fmt(data.subsample_n,0) + ' subsampled reads (' + b.unclassified_pct.toFixed(1) +
      '%) don\\'t match anything in PlusPF-8 at all. Of what does classify, it\\'s dominated by soil/sediment/groundwater bacteria — consistent with an unenriched, shipwreck-recovered sample — plus a trace of human DNA (' + b.human_pct.toFixed(2) + '%, normal handling contamination).';
    if (c.classified_pct < 1) {
      msg += ' Only ' + c.classified_pct.toFixed(3) + '% of reads classify as Coffea against the dedicated coffee database';
      msg += c.followup && c.followup.aligned === 0 ? ', and none of those actually align to the reference on closer inspection (see below) — nothing here suggests real coffee signal is being thrown away.' : '.';
    }
    bnr.textContent = msg;
  })();

  function hbarChart(containerId, items, opts) {
    const container = document.getElementById(containerId);
    if (!items.length) { container.textContent = 'No taxa to show.'; return; }
    const labelW = opts.labelW || 210, rightPad = 90, rowH = 22, topPad = 10;
    const width = 900, plotW = width - labelW - rightPad;
    const height = topPad + rowH * items.length + 10;
    const svg = elx('svg', {width: '100%', viewBox: '0 0 ' + width + ' ' + height, height: height});
    container.appendChild(svg);
    const maxV = Math.max.apply(null, items.map(function(it){ return it.value; }));
    items.forEach(function(it, i){
      const y = topPad + i * rowH;
      const barW = maxV > 0 ? (it.value / maxV) * plotW : 0;
      elx('text', {x: labelW - 8, y: y + rowH/2 + 4, class: 'chrom-label', 'text-anchor': 'end'}, svg).textContent = it.label;
      elx('rect', {x: labelW, y: y + 3, width: Math.max(barW, 2), height: rowH - 8, rx: 4, fill: it.color || 'var(--series-1)'}, svg);
      elx('text', {x: labelW + barW + 8, y: y + rowH/2 + 4, class: 'axis-label'}, svg).textContent = it.valueLabel;
      const hit = elx('rect', {x: labelW, y: y, width: Math.max(barW,2)+80, height: rowH, fill: 'transparent'}, svg);
      hit.addEventListener('pointerenter', function(evt){ showTip(evt, it.tooltip || [it.label, it.valueLabel]); });
      hit.addEventListener('pointermove', positionTip);
      hit.addEventListener('pointerleave', hideTip);
    });
  }

  (function broadChart(){
    const items = data.broad.top_taxa.map(function(t){
      return {
        label: t.name, value: t.reads,
        valueLabel: fmt(t.reads,0) + ' reads (' + t.pct.toFixed(3) + '%)',
        color: t.is_human ? 'var(--warning)' : 'var(--series-1)',
        tooltip: [t.name, fmt(t.reads,0) + ' reads directly assigned here', t.pct.toFixed(3) + '% of the ' + fmt(data.subsample_n,0) + '-read subsample', 'rank: ' + (t.rank === 'S' ? 'species' : 'genus')],
      };
    });
    hbarChart('broad-chart', items, {});
    document.getElementById('broad-note').textContent =
      'Top ' + items.length + ' species/genus-level calls by read count. Nothing plant-like appears in this list — PlusPF-8 has no plant genomes in it at all, so it couldn\\'t surface coffee DNA here even if present; that\\'s what the coffee-specific database to the right checks for.';
  })();

  (function coffeeChart(){
    const c = data.coffee;
    const items = c.species.map(function(s){
      return {
        label: s.name, value: s.reads,
        valueLabel: fmt(s.reads,0) + ' reads (' + s.pct.toFixed(4) + '%)',
        color: 'var(--series-2)',
        tooltip: [s.name, fmt(s.reads,0) + ' reads', s.pct.toFixed(4) + '% of the subsample'],
      };
    });
    if (c.unassigned_to_species > 0) {
      items.push({
        label: 'Coffea (genus only)', value: c.unassigned_to_species,
        valueLabel: fmt(c.unassigned_to_species,0) + ' reads',
        color: 'var(--baseline)',
        tooltip: ['Coffea, genus-level only', fmt(c.unassigned_to_species,0) + ' reads', 'classified to the genus but not confidently to one of the three species'],
      });
    }
    hbarChart('coffee-chart', items, {labelW: 190});
    document.getElementById('coffee-note').textContent =
      fmt(c.classified_reads,0) + ' of ' + fmt(c.total_reads,0) + ' reads (' + c.classified_pct.toFixed(3) + '%) classify somewhere under Coffea, every one of them inside Rubiaceae/Coffea specifically — nothing scattered elsewhere in the tree. Interestingly, more hits land on C. eugenioides than on the target C. arabica itself.';
  })();

  (function followup(){
    const c = data.coffee;
    if (!c.followup) return;
    document.getElementById('followup-card').style.display = '';
    const body = document.getElementById('followup-body');
    const f = c.followup;
    const regions = f.repeat_regions || [];
    const totalInRegions = regions.reduce(function(s,r){ return s + r.reads_total; }, 0);
    const characterized = regions.filter(function(r){ return r.content.indexOf('not characterized') < 0 && r.content.indexOf('not independently') < 0; });
    let html = '<div class="headline">' + f.aligned + ' of ' + fmt(c.classified_reads,0) + ' Coffea-classified reads aligned to the C. arabica reference with the pipeline\\'s own strict parameters (bwa aln -l 16500 -n 0.01) \\u2014 none did.</div>';
    html += 'Re-run with a soft-clip-tolerant aligner (bwa mem) instead, ' + (regions.length ? 'most reads do land somewhere on the reference, but ' : '') +
      'they pile up into a small number of narrow genomic windows rather than spreading out \\u2014 the table below is exactly those windows, found independently in both samples, with the reference sequence at each one pulled directly and checked for what it actually is (not guessed).' +
      (characterized.length ? ' ' + characterized.length + ' of the ' + regions.length + ' largest windows (covering ' + fmt(totalInRegions,0) + ' reads combined across both samples) turned out to be bacterial 16S ribosomal RNA gene fragments \\u2014 duplicated near-identically across several coffee nuclear chromosomes and the mitochondrial contig \\u2014 not coffee-specific sequence at all. A read that is itself an ordinary environmental-bacteria 16S fragment can pick up a stray \\u201cCoffea\\u201d k-mer hit simply because the coffee reference assembly carries its own near-identical copy of this highly conserved gene.' : '') +
      ' <strong>Conclusion: no meaningful pool of recoverable coffee signal in the unmapped fraction</strong> \\u2014 the base pipeline\\u2019s alignment step excluded the off-target reads correctly.';
    body.innerHTML = html;

    const isR0003 = data.sample.slice(-4) === 'R0003';
    document.getElementById('region-table-desc').textContent = regions.length ?
      'Where those Coffea-classified reads actually land on the reference, chain-clustered by position (\\u2264500bp between neighbors), across both R0003 and R0004 combined \\u2014 sorted by total read count.' : '';
    const tbody = document.getElementById('region-tbody');
    regions.forEach(function(r){
      const tr = document.createElement('tr');
      const sampleReads = isR0003 ? r.reads_r0003 : r.reads_r0004;
      const cells = [
        {text: r.region, cls: 'region-cell'},
        {text: fmt(sampleReads,0) + ' (' + fmt(r.reads_total,0) + ' total)', cls: 'num'},
        {text: (r.mapq0_frac*100).toFixed(0) + '%', cls: 'num'},
        {text: r.mean_nm.toFixed(1), cls: 'num'},
      ];
      cells.forEach(function(c2){
        const td = document.createElement('td');
        if (c2.cls) td.className = c2.cls;
        td.textContent = c2.text;
        tr.appendChild(td);
      });
      const contentTd = document.createElement('td');
      const notChar = r.content.indexOf('not characterized') >= 0 || r.content.indexOf('not independently') >= 0;
      const pill = document.createElement('span');
      pill.className = 'pill ' + (notChar ? 'pill-muted' : 'pill-good');
      pill.textContent = notChar ? 'not characterized' : '16S rRNA gene';
      contentTd.appendChild(pill);
      const detail = document.createElement('div');
      detail.style.marginTop = '4px';
      detail.style.color = 'var(--text-secondary)';
      detail.textContent = r.content.replace(/^bacterial /, '').replace(/^short window, low read count . not independently characterized/, 'too few reads / too short a window to check by hand');
      contentTd.appendChild(detail);
      tr.appendChild(contentTd);
      tbody.appendChild(tr);
    });
  })();

  document.querySelectorAll('.credits').forEach(function(){});
})();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
