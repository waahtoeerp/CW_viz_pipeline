#!/usr/bin/env python3
"""
Build themed HTML pages from the markdown docs in docs/ (PIPELINE.md,
simple_file_chart.md). Markdown -> HTML and mermaid diagram rendering both
happen client-side (marked.js / mermaid.js via CDN) so this script just
embeds the raw markdown source and the shared page chrome.

Usage:
  ./build_docs_pages.py
"""
import json
import os

BASE = os.path.dirname(os.path.abspath(__file__))

KNOWN_SAMPLES = ["CT600-007R0002", "CT600-007R0003", "CT600-007R0004", "CT600-007R0005"]

DOCS = [
    {
        "title": "Pipeline Run Order",
        "md_file": "docs/PIPELINE.md",
        "out_file": "pipeline.html",
        "nav_label": "Pipeline",
    },
    {
        "title": "Simple File Chart",
        "md_file": "docs/simple_file_chart.md",
        "out_file": "file-chart.html",
        "nav_label": "File Chart",
    },
]


def build_nav(current_out_file):
    viz_dir = os.path.join(BASE, "viz-data")
    nav_samples = []
    for s in KNOWN_SAMPLES:
        html_file = f"{s}_snp_viz.html"
        nav_samples.append({
            "id": s,
            "file": f"viz-data/{html_file}",
            "ready": os.path.exists(os.path.join(viz_dir, html_file)),
        })
    nav_docs = [
        {"label": d["nav_label"], "file": d["out_file"], "current": d["out_file"] == current_out_file}
        for d in DOCS
    ]
    return nav_samples, nav_docs


def render_html(title, markdown_text, nav_samples, nav_docs):
    payload = json.dumps({"markdown": markdown_text, "nav_samples": nav_samples, "nav_docs": nav_docs})
    return TEMPLATE.replace("__TITLE__", title).replace("__DATA_JSON__", payload)


TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ — The Coffee Wrecks</title>
<style>
:root {
  --font-header:    'Times New Roman', Times, serif;
  --font-body:      Arial, Helvetica, sans-serif;
  --accent:         #fed95e;
  --surface-1:      #ffffff;
  --page:           #f2e9dc;
  --text-primary:   #211d16;
  --text-secondary: #6b6459;
  --text-muted:     #8a8073;
  --gridline:       #e4dbcb;
  --border:         rgba(33,29,22,0.12);
  --link:           #7a5f00;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--text-primary); font-family: var(--font-body); }
.site-ribbon { display: flex; align-items: center; gap: 20px; background: #ffffff;
  border-bottom: 1px solid var(--border); padding: 8px 24px; position: sticky; top: 0; z-index: 50; }
.site-ribbon .ribbon-home { display: flex; align-items: center; }
.site-ribbon .ribbon-home img { height: 32px; width: auto; display: block; }
.site-ribbon .ribbon-links { display: flex; gap: 6px; flex-wrap: wrap; font-size: 13px; }
.site-ribbon .ribbon-links a { text-decoration: none; padding: 6px 14px; border-radius: 999px; color: var(--text-secondary); }
.site-ribbon .ribbon-links a:hover { background: rgba(33,29,22,0.06); }
.site-ribbon .ribbon-links a.current { background: var(--accent); color: #4a3900; font-weight: bold; }
.site-ribbon .ribbon-links a.pending { pointer-events: none; opacity: 0.5; }
.site-ribbon .ribbon-divider { width: 1px; align-self: stretch; background: var(--border); margin: 4px 2px; }
.page { max-width: 900px; margin: 0 auto; padding: 36px 24px 80px; }
h1.doc-title { font-family: var(--font-header); font-weight: 400; font-size: 28px; margin: 0 0 4px; }
.accent-rule { height: 4px; width: 60px; background: var(--accent); margin: 0 0 24px; border-radius: 2px; }
.card { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px; padding: 24px 28px; overflow-x: auto; }
.markdown-body h1 { font-family: var(--font-header); font-weight: 400; font-size: 22px; margin: 28px 0 12px; }
.markdown-body h1:first-child { margin-top: 0; }
.markdown-body h2 { font-family: var(--font-header); font-weight: 400; font-size: 19px; margin: 26px 0 10px; }
.markdown-body h3 { font-family: var(--font-header); font-weight: 400; font-size: 16px; margin: 22px 0 8px; }
.markdown-body p, .markdown-body li { font-size: 14px; line-height: 1.7; color: var(--text-primary); }
.markdown-body ul, .markdown-body ol { padding-left: 22px; }
.markdown-body li { margin-bottom: 6px; }
.markdown-body a { color: var(--link); }
.markdown-body code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12.5px;
  background: color-mix(in srgb, var(--accent) 18%, transparent); padding: 1px 5px; border-radius: 4px; }
.markdown-body pre { background: #211d16; color: #f2e9dc; padding: 14px 16px; border-radius: 8px; overflow-x: auto; }
.markdown-body pre code { background: none; color: inherit; padding: 0; }
.markdown-body table { border-collapse: collapse; width: 100%; margin: 12px 0 20px; font-size: 13px; }
.markdown-body th, .markdown-body td { border: 1px solid var(--gridline); padding: 8px 10px; text-align: left; vertical-align: top; }
.markdown-body th { background: color-mix(in srgb, var(--accent) 25%, transparent); font-weight: bold; }
.markdown-body hr { border: none; border-top: 1px solid var(--gridline); margin: 24px 0; }
.markdown-body pre.mermaid { background: var(--surface-1); border: 1px solid var(--gridline); padding: 20px; text-align: center; }
</style>
</head>
<body>
<nav class="site-ribbon">
  <a class="ribbon-home" href="index.html" title="Back to front page">
    <img src="theme/logo-compass.png" alt="The Coffee Wrecks — back to front page">
  </a>
  <div class="ribbon-links" id="ribbon-links"></div>
</nav>
<div class="page">
  <h1 class="doc-title">__TITLE__</h1>
  <div class="accent-rule"></div>
  <div class="card">
    <div class="markdown-body" id="markdown-body"></div>
  </div>
</div>

<script type="application/json" id="page-data">__DATA_JSON__</script>
<script src="https://cdn.jsdelivr.net/npm/marked@9.1.6/marked.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/mermaid@10.9.1/dist/mermaid.min.js"></script>
<script>
(function(){
  const data = JSON.parse(document.getElementById('page-data').textContent);

  (function ribbon(){
    const wrap = document.getElementById('ribbon-links');
    const home = document.createElement('a');
    home.href = 'index.html';
    home.textContent = 'Front page';
    wrap.appendChild(home);
    data.nav_samples.forEach(function(s){
      const a = document.createElement('a');
      a.textContent = s.id;
      if (s.ready) { a.href = s.file; } else { a.className = 'pending'; a.href = '#'; }
      wrap.appendChild(a);
    });
    const divider = document.createElement('div');
    divider.className = 'ribbon-divider';
    wrap.appendChild(divider);
    data.nav_docs.forEach(function(d){
      const a = document.createElement('a');
      a.textContent = d.label;
      a.href = d.file;
      if (d.current) a.className = 'current';
      wrap.appendChild(a);
    });
  })();

  function escapeHtml(s) {
    return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  const renderer = new marked.Renderer();
  const originalCode = renderer.code.bind(renderer);
  renderer.code = function(code, infostring) {
    if ((infostring || '').trim() === 'mermaid') {
      return '<pre class="mermaid">' + escapeHtml(code) + '</pre>';
    }
    return originalCode(code, infostring);
  };
  marked.use({ renderer: renderer });

  document.getElementById('markdown-body').innerHTML = marked.parse(data.markdown);

  mermaid.initialize({ startOnLoad: false, theme: 'neutral', securityLevel: 'loose' });
  mermaid.run({ querySelector: '.mermaid' });
})();
</script>
</body>
</html>
"""


def main():
    for doc in DOCS:
        with open(os.path.join(BASE, doc["md_file"])) as f:
            markdown_text = f.read()
        nav_samples, nav_docs = build_nav(doc["out_file"])
        html = render_html(doc["title"], markdown_text, nav_samples, nav_docs)
        out_path = os.path.join(BASE, doc["out_file"])
        with open(out_path, "w") as f:
            f.write(html)
        print(f"Wrote {out_path} ({len(html):,} bytes)")


if __name__ == "__main__":
    main()
