"""Download candidate datasets (CSV) and their documentation pages from Rdatasets (agridat package)."""
import sys, re, html, urllib.request, pathlib
BASE = "https://vincentarelbundock.github.io/Rdatasets"
root = pathlib.Path(__file__).parent
def fetch(item, pkg="agridat"):
    for sub in ("raw", "docs"):
        (root / sub).mkdir(exist_ok=True)
    for kind, ext, sub in (("csv", "csv", "raw"), ("doc", "html", "docs")):
        out = root / sub / f"{item}.{ext}"
        if not out.exists():
            out.write_bytes(urllib.request.urlopen(f"{BASE}/{kind}/{pkg}/{item}.{ext}", timeout=60).read())
def doc_text(item):
    t = (root / "docs" / f"{item}.html").read_text(encoding="utf-8", errors="replace")
    t = re.sub(r"(?s)<(script|style).*?</\1>", "", t)
    t = re.sub(r"(?s)<h3>Examples</h3>.*", "", t)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return re.sub(r"[ \t]+", " ", re.sub(r"\n\s*\n+", "\n", t)).strip()
if __name__ == "__main__":
    import pandas as pd
    for item in sys.argv[1:]:
        fetch(item)
        d = pd.read_csv(root / "raw" / f"{item}.csv").drop(columns="rownames", errors="ignore")
        print("=" * 100); print(item, d.shape)
        print(doc_text(item)[:1800])
        for c in d.columns:
            u = d[c].nunique()
            if u <= 12: print(f"  {c}: {dict(d[c].value_counts().sort_index())}")
            else: print(f"  {c}: {u} unique, min={d[c].min()}, max={d[c].max()}" + (f", NA={d[c].isna().sum()}" if d[c].isna().any() else ""))
