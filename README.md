# Primitive quartic trinomials over finite fields — computational artifacts

Computational code, exact verification tools, integrity checks, and archival metadata for certified computations concerning primitive quartic trinomials over finite fields.

Authors: **Murat Alan** and **Kadriye Simsek Alan**.

## Main computational result

The certified finite sweep covers prime powers `q = p^m` of characteristic `p > 3` through

`q = 2.6 × 10^10`.

The exceptional field sizes are

`q = 5, 7, 13`.

For the interval `(10^8, 2.6×10^10]`, the completed production run processed:

- 1,127,956,365 primes;
- 13,992 proper prime powers;
- 1,127,970,357 fields in total.

Every field in this interval passed the partial-factorisation criterion (`C4F`), with no `C4`, `DS`, or unresolved case.

## Repository contents

```text
code/
  verify_exact.py
  sweep_quartic.py
  run_all.py

data/
  expected_output_verify.txt
  DATA.md
  LICENSE

tools/
  verify_result_integrity.ps1
  VERIFY_RESULT.bat
  build_zenodo_archives.py
  BUILD_ZENODO_ARCHIVES.bat
  README_TOOLS.md

CITATION.cff
LICENSE
.zenodo.json
.gitignore
```

## Requirements

Python **3.8+**. No third-party Python dependencies are required.

## Exact verification

```bash
cd code
python3 verify_exact.py > got.txt
```

Compare `got.txt` byte-for-byte with:

```text
data/expected_output_verify.txt
```

## Small sweep sanity test

```bash
cd code
python3 sweep_quartic.py sanity
```

## Full complementary sweep

```bash
cd code
python3 run_all.py --lo 100000000 --hi 26000000000 --procs 8 --outdir results
```

This covers the interval `(10^8, 2.6×10^10]`.

A completed run contains 1559 raw output files:

- 518 `part_*.csv.gz`;
- 518 matching `.stats`;
- 518 matching `.done`;
- `pp.csv.gz`, `pp.stats`, `pp.done`;
- `summary.txt`;
- an empty `failures.csv`.

## Integrity verification

The tools in `tools/` verify:

- the expected 1559-file structure;
- the summary totals;
- complete decompression of all 519 `.csv.gz` files;
- SHA-256 manifests.

The completed production run returned:

```text
OVERALL: PASS
```

## Certificate data

The certificate datasets are archived separately on Zenodo. The Git repository contains the code, verification tools, and documentation needed to reproduce and validate them.

Zenodo version DOI (v1.0.0): **10.5281/zenodo.22135054**  
Zenodo concept DOI (all versions): **10.5281/zenodo.22135053**

## Licensing

- Code: MIT (`LICENSE`).
- Certificate data: CC-BY-4.0 (`data/LICENSE`).

## Citation

See `CITATION.cff`.
