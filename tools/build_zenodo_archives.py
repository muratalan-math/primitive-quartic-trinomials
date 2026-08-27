#!/usr/bin/env python3
"""
Build logically partitioned Zenodo upload archives from a completed run_all
result directory.

Run this script from inside the completed result directory, or pass the result
directory as the first argument.

Required precondition:
  RUN_ALL_INTEGRITY_REPORT.txt exists and contains "OVERALL: PASS".

The script never modifies raw result files. It creates ./zenodo_upload/.
"""
from pathlib import Path
import sys, zipfile, hashlib, shutil, re

RESULT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd().resolve()
REPORT = RESULT / "RUN_ALL_INTEGRITY_REPORT.txt"
if not REPORT.exists():
    raise SystemExit("Missing RUN_ALL_INTEGRITY_REPORT.txt")
if "OVERALL: PASS" not in REPORT.read_text(encoding="utf-8-sig", errors="replace"):
    raise SystemExit("Integrity report does not contain OVERALL: PASS")

OUT = RESULT / "zenodo_upload"
OUT.mkdir(exist_ok=True)

groups = [
    ("certificates_100m_5b.zip", 100_000_000, 5_000_000_000),
    ("certificates_5b_10b.zip", 5_000_000_000, 10_000_000_000),
    ("certificates_10b_15b.zip", 10_000_000_000, 15_000_000_000),
    ("certificates_15b_20b.zip", 15_000_000_000, 20_000_000_000),
    ("certificates_20b_26b.zip", 20_000_000_000, 26_000_000_000),
]

parts = {}
rx = re.compile(r"^part_(\d+)\.csv\.gz$")
for p in RESULT.glob("part_*.csv.gz"):
    m = rx.match(p.name)
    if m:
        parts[int(m.group(1))] = p

if len(parts) != 518:
    raise SystemExit(f"Expected 518 part CSV.GZ files, found {len(parts)}")

def add_stored(zf, path):
    zf.write(path, arcname=path.name, compress_type=zipfile.ZIP_STORED)

for zname, lo, hi in groups:
    selected = sorted(q for q in parts if lo <= q < hi)
    expected = (hi - lo) // 50_000_000
    if len(selected) != expected:
        raise SystemExit(f"{zname}: expected {expected} parts, found {len(selected)}")
    zp = OUT / zname
    with zipfile.ZipFile(zp, "w", allowZip64=True) as zf:
        for q in selected:
            gz = RESULT / f"part_{q}.csv.gz"
            st = RESULT / f"part_{q}.stats"
            if not st.exists():
                raise SystemExit(f"Missing {st.name}")
            add_stored(zf, gz)
            add_stored(zf, st)

ppzip = OUT / "prime_powers_100m_26b.zip"
with zipfile.ZipFile(ppzip, "w", allowZip64=True) as zf:
    for name in ("pp.csv.gz", "pp.stats"):
        p = RESULT / name
        if not p.exists():
            raise SystemExit(f"Missing {name}")
        add_stored(zf, p)

for name in (
    "summary.txt",
    "failures.csv",
    "RUN_ALL_INTEGRITY_REPORT.txt",
    "MANIFEST_SHA256_ALL.txt",
    "MANIFEST_SHA256_ZENODO.txt",
    "ZENODO_FILE_LIST.txt",
):
    p = RESULT / name
    if p.exists():
        shutil.copy2(p, OUT / name)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

entries = []
for p in sorted(OUT.iterdir()):
    if p.is_file() and p.name != "ZENODO_ARCHIVE_SHA256.txt":
        entries.append((sha256(p), p.stat().st_size, p.name))

manifest = OUT / "ZENODO_ARCHIVE_SHA256.txt"
with manifest.open("w", encoding="utf-8") as f:
    for h, n, name in entries:
        f.write(f"{h}  {n}  {name}\n")

print(f"Built Zenodo upload staging directory: {OUT}")
for _, n, name in entries:
    print(f"{name:42s} {n:,} bytes")
print("ZENODO_ARCHIVE_SHA256.txt created.")
