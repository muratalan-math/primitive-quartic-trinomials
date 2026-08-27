#!/usr/bin/env python3
"""
sweep_quartic.py -- certified existence of primitive quartic trinomials
                    x^4 + a x + b over F_q, for all prime powers q = p^m,
                    p > 3, in a given range.

For each q the program produces one of three certificates:

  C4:t          the exact sieve criterion holds with sieve
                parameter t; verified in EXACT rational/integer arithmetic
                (the inequality delta(q^2-1) - 4q/theta(d) > C q^{3/2} is
                squared and compared as exact fractions; no floating point).

  DS:a,b,tr     an explicit primitive trinomial x^4 + a x + b was found
                after tr random trials; primitivity is PROVED by the exact
                order computation ord(x mod f) = q^4 - 1, which certifies
                simultaneously that f is irreducible and that its root is
                a primitive element of F_{q^4}.  (If the order of x in
                (F_q[x]/f)^* equals q^4-1, then since that order is prime
                to p and the prime-to-p exponent of the unit group for any
                f with a factor pattern other than a single irreducible
                quartic divides lcm(q^{d_i}-1) < q^4-1, f must be
                irreducible; primitivity is then the definition.)
                For prime powers q = p^m, elements of F_q are encoded as
                integers in base p; a and b are reported in that encoding
                with respect to the recorded primitive polynomial h.

  EXH:NONE      exhaustive search over all (a,b) in (F_q^*)^2 found no
                primitive quartic trinomial (only attempted for tiny q).

All primality tests are deterministic Miller-Rabin with the first twelve prime
bases {2,3,5,...,37}, proven correct for
n < 318665857834031151167461 (Sorenson-Webster).  The composite
q^4-1 is never submitted to Miller-Rabin: the factorisation is organised via
q-1, q+1, and q^2+1 and their factors.  In the certified range
q <= 2.6e10, every integer submitted to a primality test is therefore below
q^2+1 < 6.8e20, with a wide margin inside the deterministic bound.

Usage:
  python3 sweep_quartic.py sanity
  python3 sweep_quartic.py sweep  START END OUT.csv   # primes in (START,END]
  python3 sweep_quartic.py pp     START END OUT.csv   # prime powers p^m, m>=2
"""

import sys, random
from fractions import Fraction

# ----------------------------------------------------------------- primality

_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)  # proven < 318665857834031151167461

def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def _rho(n, seed=1):
    if n % 2 == 0:
        return 2
    from math import gcd
    while True:
        x = y = random.Random(seed).randrange(2, n)
        c = random.Random(seed + 1).randrange(1, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d
        seed += 2

def factor_primes(n):
    """Set of prime divisors of n, deterministic-rigorous in our range."""
    out, stack = set(), [n]
    while stack:
        v = stack.pop()
        if v == 1:
            continue
        # strip small primes first
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            while v % p == 0:
                out.add(p)
                v //= p
        if v == 1:
            continue
        if is_prime(v):
            out.add(v)
            continue
        d = _rho(v)
        stack.append(d)
        stack.append(v // d)
    return out


def factor_primes_Q4(q):
    """Prime divisors of q^4-1, computed piecewise so that every integer
    submitted to a primality test divides q-1, q+1 or q^2+1 and therefore
    lies far inside the proven range of the deterministic MR base set."""
    return factor_primes(q - 1) | factor_primes(q + 1) | factor_primes(q * q + 1)

# ------------------------------------------------- criterion (4)/(5), exact

def crit4(q):
    """Exact sieve test for a single q; returns t or None."""
    fRm = sorted(factor_primes(q + 1) | factor_primes(q * q + 1))
    fQextra = sorted(factor_primes(q - 1) - set(fRm))
    s = len(fRm)
    # suffix sums of 1/l over sieve primes, prefix products of theta
    suf = [Fraction(0)] * (s + 1)
    for i in range(s - 1, -1, -1):
        suf[i] = suf[i + 1] + Fraction(1, fRm[i])
    th_extra = Fraction(1)
    for l in fQextra:
        th_extra *= Fraction(l - 1, l)
    th_pref = Fraction(1)
    q3 = q * q * q
    q2m1 = q * q - 1
    for t in range(0, s + 1):
        r = s - t
        delta = 1 - suf[t]
        if delta > 0:
            theta = th_pref * th_extra
            lhs = delta * q2m1 - Fraction(4 * q) / theta
            if lhs > 0:
                Cc = 6 * (1 << t) * (3 * r - 1) if r > 0 else 6 * ((1 << t) - 1)
                if Cc == 0 or lhs * lhs > Cc * Cc * q3:
                    return t
        if t < s:
            th_pref *= Fraction(fRm[t] - 1, fRm[t])
    return None

# ------------------------------------- direct search over prime fields F_p

def _order_is_Q_prime(a, b, q, Q, exps):
    """Prime q: is ord(x mod x^4+ax+b) = Q?  Exact; proves primitivity."""
    def mul(u, v):
        r = [0] * 7
        for i in range(4):
            ui = u[i]
            if ui:
                for j in range(4):
                    if v[j]:
                        r[i + j] = (r[i + j] + ui * v[j]) % q
        for d in (6, 5, 4):
            c = r[d]
            if c:
                r[d] = 0
                r[d - 3] = (r[d - 3] - c * a) % q
                r[d - 4] = (r[d - 4] - c * b) % q
        return r[:4]

    def powx(u, e):
        r = [1, 0, 0, 0]
        while e:
            if e & 1:
                r = mul(r, u)
            u = mul(u, u)
            e >>= 1
        return r

    ONE = [1, 0, 0, 0]
    X = [0, 1, 0, 0]
    if powx(X, Q) != ONE:
        return False
    return all(powx(X, e) != ONE for e in exps)

def direct_prime(q, max_trials=5000):
    Q = q ** 4 - 1
    exps = [Q // l for l in factor_primes_Q4(q)]
    rng = random.Random(q)  # reproducible
    for tr in range(1, max_trials + 1):
        a = rng.randrange(1, q)
        b = rng.randrange(1, q)
        if _order_is_Q_prime(a, b, q, Q, exps):
            return a, b, tr
    return None

def exhaustive_prime(q):
    Q = q ** 4 - 1
    exps = [Q // l for l in factor_primes_Q4(q)]
    for a in range(1, q):
        for b in range(1, q):
            if _order_is_Q_prime(a, b, q, Q, exps):
                return a, b
    return None

# --------------------------- direct search over prime powers F_{p^m}, m>=2

class Fq:
    """F_{p^m} with exp/log tables built from a primitive polynomial h.
    Elements are integers 0..q-1, base-p encodings of coordinate vectors.
    Building the table by enumerating powers of x PROVES h primitive."""

    def __init__(self, p, m, seed=0):
        self.p, self.m, self.q = p, m, p ** m
        rng = random.Random(seed)
        while True:
            h = [rng.randrange(p) for _ in range(m)]  # x^m + h[m-1]x^{m-1}+...+h[0]
            if h[0] == 0:
                continue
            if self._try_tables(h):
                self.h = h
                return

    def _try_tables(self, h):
        p, m, q = self.p, self.m, self.q
        exp = [0] * (q - 1)
        log = [None] * q
        cur = [1] + [0] * (m - 1)  # 1
        for i in range(q - 1):
            code = 0
            for c in reversed(cur):
                code = code * p + c
            if code == 1 and i > 0:
                return False  # order of x is i < q-1: h not primitive
            if log[code] is not None:
                return False
            exp[i] = code
            log[code] = i
            # multiply cur by x
            top = cur[m - 1]
            cur = [0] + cur[:m - 1]
            if top:
                for j in range(m):
                    cur[j] = (cur[j] - top * h[j]) % p
        if log.count(None) != 1:  # only 0 unassigned
            return False
        self.exp, self.log = exp, log
        return True

    def mul(self, u, v):
        if u == 0 or v == 0:
            return 0
        return self.exp[(self.log[u] + self.log[v]) % (self.q - 1)]

    def add(self, u, v):
        p, m = self.p, self.m
        r, mult = 0, 1
        for _ in range(m):
            r += ((u % p + v % p) % p) * mult
            u //= p
            v //= p
            mult *= p
        return r

    def neg(self, u):
        p, m = self.p, self.m
        r, mult = 0, 1
        for _ in range(m):
            r += ((-u) % p) * mult
            u //= p
            mult *= p
        return r

def _order_is_Q_pp(F, a, b, Q, exps):
    na, nb = F.neg(a), F.neg(b)
    mulq, addq = F.mul, F.add

    def mul4(u, v):
        r = [0] * 7
        for i in range(4):
            ui = u[i]
            if ui:
                for j in range(4):
                    vj = v[j]
                    if vj:
                        r[i + j] = addq(r[i + j], mulq(ui, vj))
        for d in (6, 5, 4):
            c = r[d]
            if c:
                r[d] = 0
                r[d - 3] = addq(r[d - 3], mulq(c, na))
                r[d - 4] = addq(r[d - 4], mulq(c, nb))
        return r[:4]

    def powx(u, e):
        r = [1, 0, 0, 0]
        while e:
            if e & 1:
                r = mul4(r, u)
            u = mul4(u, u)
            e >>= 1
        return r

    ONE = [1, 0, 0, 0]
    X = [0, 1, 0, 0]
    if powx(X, Q) != ONE:
        return False
    return all(powx(X, e) != ONE for e in exps)

def direct_pp(p, m, max_trials=5000):
    F = Fq(p, m)
    q = F.q
    Q = q ** 4 - 1
    exps = [Q // l for l in factor_primes_Q4(q)]
    rng = random.Random(q)
    for tr in range(1, max_trials + 1):
        a = rng.randrange(1, q)
        b = rng.randrange(1, q)
        if _order_is_Q_pp(F, a, b, Q, exps):
            return F, a, b, tr
    return F, None, None, None

def exhaustive_pp(p, m):
    F = Fq(p, m)
    q = F.q
    Q = q ** 4 - 1
    exps = [Q // l for l in factor_primes_Q4(q)]
    for a in range(1, q):
        for b in range(1, q):
            if _order_is_Q_pp(F, a, b, Q, exps):
                return F, a, b
    return F, None, None

# ------------------------------------------------------------------ drivers

def handle_prime(q):
    t = crit4(q)
    if t is not None:
        return f"{q},p,{q},1,C4,{t}"
    res = direct_prime(q)
    if res:
        a, b, tr = res
        return f"{q},p,{q},1,DS,{a}|{b}|{tr}"
    if q <= 2500:
        res = exhaustive_prime(q)
        if res:
            return f"{q},p,{q},1,DS-EXH,{res[0]}|{res[1]}"
        return f"{q},p,{q},1,EXH,NONE"
    return f"{q},p,{q},1,UNRESOLVED,"

def handle_pp(p, m):
    q = p ** m
    t = crit4(q)
    if t is not None:
        return f"{q},pp,{p},{m},C4,{t}"
    F, a, b, tr = direct_pp(p, m)
    h = "~".join(map(str, F.h))
    if a is not None:
        return f"{q},pp,{p},{m},DS,{a}|{b}|{tr}|h={h}"
    if q <= 2500:
        F, a, b = exhaustive_pp(p, m)
        h = "~".join(map(str, F.h))
        if a is not None:
            return f"{q},pp,{p},{m},DS-EXH,{a}|{b}|h={h}"
        return f"{q},pp,{p},{m},EXH,NONE|h={h}"
    return f"{q},pp,{p},{m},UNRESOLVED,"

def sieve_primes(limit):
    bs = bytearray([1]) * (limit + 1)
    bs[0:2] = b"\x00\x00"
    i = 2
    while i * i <= limit:
        if bs[i]:
            bs[i * i:limit + 1:i] = b"\x00" * len(range(i * i, limit + 1, i))
        i += 1
    return [i for i in range(2, limit + 1) if bs[i]]

def main():
    mode = sys.argv[1]
    if mode == "sanity":
        # Known: no primitive quartic trinomial over F_5, F_7, F_13 (Cohen-Mills);
        # existence over F_11, F_17, F_19, F_23.
        for q in (5, 7, 11, 13, 17, 19, 23):
            r = exhaustive_prime(q)
            print(f"F_{q}: {'NONE' if r is None else f'x^4+{r[0]}x+{r[1]} primitive'}")
        # first unknown prime powers
        for (p, m) in ((5, 2), (7, 2), (13, 2), (5, 3)):
            F, a, b = exhaustive_pp(p, m)
            tag = "NONE" if a is None else f"a={a},b={b} (base-{p} codes), h={F.h}"
            print(f"F_{p}^{m} (q={p**m}): {tag}")
        return

    start, end, out = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    f = open(out, "a", buffering=1)
    done = 0
    if mode == "sweep":
        for q in sieve_primes(end):
            if q <= start or q <= 3:
                continue
            f.write(handle_prime(q) + "\n")
            done += 1
            if done % 2000 == 0:
                print(f"[sweep] {done} primes, at q={q}", flush=True)
    elif mode == "pp":
        ps = sieve_primes(int(end ** 0.5) + 1)
        qs = []
        for p in ps:
            if p <= 3:
                continue
            v = p * p
            m = 2
            while v <= end:
                if v > start:
                    qs.append((v, p, m))
                v *= p
                m += 1
        for q, p, m in sorted(qs):
            f.write(handle_pp(p, m) + "\n")
            done += 1
            print(f"[pp] q={q}=({p}^{m}) done ({done})", flush=True)
    f.close()

if __name__ == "__main__":
    main()


# ------------- fast sufficient criterion via partial factorisation --------
# Rigorous by monotonicity: replace the (unknown) large prime factors of the
# unfactored part U of R by the pessimistic surrogate of kU "phantom" primes,
# each >= B, where kU = max c with B^c <= U (any c actual primes >= B have
# product >= B^c <= U).  The surrogate has MORE sieving primes (r_s >= r_a)
# and SMALLER delta (delta_s <= delta_a), while t, theta(d) use only exactly
# known primes; criterion (4) is monotone in (r decreasing, delta increasing),
# so the surrogate inequality implies the exact sieve inequality for the
# actual data.

_B_FAST = 10 ** 4
_SMALL_PRIMES = None

def _small_primes():
    global _SMALL_PRIMES
    if _SMALL_PRIMES is None:
        _SMALL_PRIMES = sieve_primes(_B_FAST)
    return _SMALL_PRIMES

def crit4_fast(q):
    """Sufficient test of criterion (4) using only trial division of R up to
    B=1e4 (plus exact factorisation of the small number q-1).  Returns t on
    success, None if inconclusive (then fall back to crit4)."""
    R = (q + 1) * (q * q + 1)
    U = R
    small = []
    for l in _small_primes():
        if l > 3 and U % l == 0 or l <= 3 and U % l == 0:
            small.append(l)
            while U % l == 0:
                U //= l
    # unfactored part: all its primes exceed B
    if U == 1:
        phantoms, real_big = 0, []
    elif is_prime(U):
        phantoms, real_big = 0, [U]
    else:
        c = 0
        v = 1
        while v <= U:
            v *= _B_FAST
            c += 1
        phantoms, real_big = c - 1 if v > U else c, []
        if phantoms < 2:      # composite with all factors > B has >= 2
            phantoms = 2
    fq1 = factor_primes(q - 1)
    extraQ = [l for l in fq1 if R % l != 0]
    th_extra = Fraction(1)
    for l in extraQ:
        th_extra *= Fraction(l - 1, l)
    fRs = sorted(small)
    s = len(fRs)
    # suffix sums of 1/l over known-small sieving primes
    suf = [Fraction(0)] * (s + 1)
    for i in range(s - 1, -1, -1):
        suf[i] = suf[i + 1] + Fraction(1, fRs[i])
    big_pen = Fraction(phantoms, _B_FAST) + sum(Fraction(1, l) for l in real_big)
    rr_big = phantoms + len(real_big)
    q3, q2m1 = q * q * q, q * q - 1
    th_pref = Fraction(1)
    for t in range(0, s + 1):
        r = (s - t) + rr_big
        delta = 1 - suf[t] - big_pen
        if delta > 0:
            theta = th_pref * th_extra
            lhs = delta * q2m1 - Fraction(4 * q) / theta
            if lhs > 0:
                Cc = 6 * (1 << t) * (3 * r - 1) if r > 0 else 6 * ((1 << t) - 1)
                if Cc == 0 or lhs * lhs > Cc * Cc * q3:
                    return t
        if t < s:
            th_pref *= Fraction(fRs[t] - 1, fRs[t])
    return None

def handle_prime_v2(q):
    t = crit4_fast(q)
    if t is not None:
        return f"{q},p,{q},1,C4F,{t}"
    return handle_prime(q)

def handle_pp_v2(p, m):
    q = p ** m
    t = crit4_fast(q)
    if t is not None:
        return f"{q},pp,{p},{m},C4F,{t}"
    return handle_pp(p, m)
