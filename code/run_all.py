#!/usr/bin/env python3
"""
run_all.py -- parallel, resumable driver for the remaining range of the
primitive-quartic-trinomial sweep:  1e8 < q <= 2.6e10.

Usage (from the directory containing sweep_quartic.py):

    python3 run_all.py --lo 100000000 --hi 26000000000 --procs 8 --outdir results

Re-running the same command resumes: finished parts (marked .done) are
skipped, an interrupted part is redone from its start.  Output per part is
gzipped CSV (multi-member gzip is valid), plus a .stats sidecar; prime
powers p^m (m>=2) are handled separately and written to pp.csv.gz with a
pp.stats sidecar.  At the end a summary.txt (primes, prime powers, and
total) and failures.csv are written in --outdir.

Archival outputs: keep summary.txt, failures.csv, the per-part certificate
streams, pp.csv.gz, and their statistics/checksum metadata. The large
certificate data are intended for a persistent data archive such as Zenodo.
"""
import argparse, gzip, json, math, os, sys, time
from multiprocessing import Pool
from sweep_quartic import handle_prime_v2, handle_pp_v2

PART_WIDTH = 50_000_000
SEG = 5_000_000

def base_primes(limit):
    bs = bytearray([1]) * (limit + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if bs[i]:
            bs[i * i:limit + 1:i] = b"\x00" * len(range(i * i, limit + 1, i))
    return [i for i in range(2, limit + 1) if bs[i]]

def iter_primes_segmented(lo, hi, bps):
    """Primes in (lo, hi], segmented sieve with base primes bps (>= sqrt(hi))."""
    lo = max(lo, 1)
    seg_lo = lo + 1
    while seg_lo <= hi:
        seg_hi = min(seg_lo + SEG - 1, hi)
        n = seg_hi - seg_lo + 1
        bs = bytearray([1]) * n
        for p in bps:
            if p * p > seg_hi:
                break
            start = max(p * p, ((seg_lo + p - 1) // p) * p)
            bs[start - seg_lo::p] = b"\x00" * len(range(start, seg_hi + 1, p))
        for i in range(n):
            if bs[i]:
                yield seg_lo + i
        seg_lo = seg_hi + 1

def do_part(args):
    plo, phi, outdir = args
    tag = f"part_{plo}"
    done = os.path.join(outdir, tag + ".done")
    if os.path.exists(done):
        return tag, None
    bps = base_primes(int(math.isqrt(phi)) + 1)
    t0 = time.time()
    stats = {"lo": plo, "hi": phi, "n": 0, "C4F": 0, "C4": 0, "DS": 0,
             "other": 0, "fails": []}
    with gzip.open(os.path.join(outdir, tag + ".csv.gz"), "wt",
                   compresslevel=1) as f:
        for q in iter_primes_segmented(plo, phi, bps):
            line = handle_prime_v2(q)
            f.write(line + "\n")
            m = line.split(",")[4]
            stats["n"] += 1
            if m in ("C4F", "C4", "DS"):
                stats[m] += 1
                if m != "C4F":
                    stats["fails"].append(line)
            else:
                stats["other"] += 1
                stats["fails"].append(line)
    stats["secs"] = round(time.time() - t0, 1)
    with open(os.path.join(outdir, tag + ".stats"), "w") as f:
        json.dump(stats, f)
    open(done, "w").close()
    return tag, stats

def do_pp(lo, hi, outdir):
    """Certify all prime powers p^m (p>3, m>=2) in (lo, hi]."""
    done = os.path.join(outdir, "pp.done")
    stats_path = os.path.join(outdir, "pp.stats")
    if os.path.exists(done) and os.path.exists(stats_path):
        return json.load(open(stats_path))
    bps = base_primes(int(math.isqrt(hi)) + 1)
    t0 = time.time()
    stats = {"n": 0, "C4F": 0, "C4": 0, "DS": 0, "other": 0, "fails": []}
    with gzip.open(os.path.join(outdir, "pp.csv.gz"), "wt",
                   compresslevel=1) as f:
        for p in bps:
            if p <= 3:
                continue
            v, m = p * p, 2
            while v <= hi:
                if v > lo:
                    line = handle_pp_v2(p, m)
                    f.write(line + "\n")
                    meth = line.split(",")[4]
                    stats["n"] += 1
                    if meth in ("C4F", "C4", "DS"):
                        stats[meth] += 1
                        if meth != "C4F":
                            stats["fails"].append(line)
                    else:
                        stats["other"] += 1
                        stats["fails"].append(line)
                v *= p
                m += 1
    stats["secs"] = round(time.time() - t0, 1)
    with open(stats_path, "w") as fh:
        json.dump(stats, fh)
    open(done, "w").close()
    return stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lo", type=int, default=100_000_000)
    ap.add_argument("--hi", type=int, default=26_000_000_000)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    parts = []
    lo = a.lo
    while lo < a.hi:
        hi = min(lo + PART_WIDTH, a.hi)
        parts.append((lo, hi, a.outdir))
        lo = hi
    todo = [p for p in parts
            if not os.path.exists(os.path.join(a.outdir, f"part_{p[0]}.done"))]
    print(f"{len(parts)} parts total, {len(todo)} to do, {a.procs} processes",
          flush=True)

    t0 = time.time()
    ndone = len(parts) - len(todo)
    with Pool(a.procs) as pool:
        for tag, st in pool.imap_unordered(do_part, todo):
            ndone += 1
            el = time.time() - t0
            eta = el / max(1, ndone - (len(parts) - len(todo))) * \
                (len(parts) - ndone)
            print(f"[{ndone}/{len(parts)}] {tag} done "
                  f"({st['n']} primes, {st['secs']}s)  ETA {eta/3600:.1f}h",
                  flush=True)

    print("prime parts complete; running prime powers...", flush=True)
    pp = do_pp(a.lo, a.hi, a.outdir)
    print(f"prime powers done ({pp['n']} fields, {pp.get('secs', 0)}s)",
          flush=True)

    prime_tot = {"n": 0, "C4F": 0, "C4": 0, "DS": 0, "other": 0}
    fails = []
    for p in parts:
        st = json.load(open(os.path.join(a.outdir, f"part_{p[0]}.stats")))
        for k in prime_tot:
            prime_tot[k] += st[k]
        fails += st["fails"]
    fails += pp["fails"]
    tot = {k: prime_tot[k] + pp[k] for k in prime_tot}

    with open(os.path.join(a.outdir, "failures.csv"), "w") as f:
        f.write("\n".join(fails) + ("\n" if fails else ""))
    with open(os.path.join(a.outdir, "summary.txt"), "w") as f:
        f.write(
            f"range ({a.lo}, {a.hi}]\n"
            f"\n"
            f"primes processed        : {prime_tot['n']}\n"
            f"prime powers (m>=2)     : {pp['n']}\n"
            f"total fields processed  : {tot['n']}\n"
            f"\n"
            f"C4F (partial-factorisation criterion) : {tot['C4F']}\n"
            f"C4  (full-factorisation criterion)    : {tot['C4']}\n"
            f"DS  (explicit witness)                : {tot['DS']}\n"
            f"other/UNRESOLVED                      : {tot['other']}\n"
            f"lines in failures.csv (non-C4F)       : {len(fails)}\n"
        )
    print(open(os.path.join(a.outdir, "summary.txt")).read())

if __name__ == "__main__":
    main()
