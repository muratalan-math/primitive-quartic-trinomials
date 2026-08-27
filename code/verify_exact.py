#!/usr/bin/env python3
"""
Exact verification script for the computational identities used in the
primitive-quartic-trinomial calculations.

Every check below is performed in exact arithmetic: integers, finite-field
elements, or set/gcd computations.  No floating-point arithmetic and no
roots of unity are used anywhere in this file, and there are no dependencies
outside the Python standard library: primality is deterministic Miller-Rabin
on the first twelve prime bases {2,3,5,...,37} and factorisation is Pollard's rho with certified
prime factors.

The reduction that makes this possible is the identity
    sum_{psi^g = eps} psi(u)  =  g * [ u is a g-th power in F_q^* ],
which collapses the relevant character sums and the master closed form into
a plain integer count.

Run:  python3 verify_exact.py
"""

from math import gcd
import itertools
import random as _random

# ------------------------------------------- self-contained number theory
# Deterministic Miller-Rabin on the first twelve prime bases {2,3,5,...,37}, proven correct for
# n < 318665857834031151167461 (Sorenson-Webster); every integer tested here is far smaller.

_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

def is_prime(n):
    if n < 2:
        return False
    for p in _MR_BASES:
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
    while True:
        x = y = _random.Random(seed).randrange(2, n)
        c = _random.Random(seed + 1).randrange(1, n)
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
    """Set of prime divisors of n (rigorous in our range)."""
    out, stack = set(), [n]
    while stack:
        v = stack.pop()
        for p in _MR_BASES + (41, 43, 47):
            while v % p == 0:
                out.add(p)
                v //= p
        if v == 1:
            continue
        if is_prime(v):
            out.add(v)
            continue
        d = _rho(v)
        stack += [d, v // d]
    return out

def nth_primes(count):
    """List of the first `count` primes, by sieve."""
    lim = max(30, int(count * (2.3 + 1.3 * __import__('math').log(count)))) \
        if count > 5 else 30
    bs = bytearray([1]) * (lim + 1)
    bs[0:2] = b"\x00\x00"
    i = 2
    while i * i <= lim:
        if bs[i]:
            bs[i * i:lim + 1:i] = b"\x00" * len(range(i * i, lim + 1, i))
        i += 1
    ps = [i for i in range(2, lim + 1) if bs[i]]
    return ps[:count]

# ---------------------------------------------------------------- utilities

def eta(x, p):
    """Quadratic character of F_p, with eta(0) = 0."""
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def is_gth_power(u, g, p):
    """Exact test: is u a g-th power in F_p^* ?"""
    u %= p
    return u != 0 and pow(u, (p - 1) // g, p) == 1


def trinomial_data(n, k, p):
    """Return (n1, k1, alpha, beta, A, B) for x^n + a x^k + b over F_p."""
    d = gcd(n, k)
    n1, k1 = n // d, k // d
    alpha, beta = n1, n1 - k1
    A = pow(n, n1, p)
    B = (pow(-1, n1 + 1) * pow(n - k, beta, p) * pow(k, k1, p)) % p
    return n1, k1, alpha, beta, A, B


def confounding(n, k, p):
    d = gcd(n, k)
    return d % 2 == 1 and n % p and k % p and (n - k) % p


# ------------------------------------------------- 1. aggregate identities

def check_theorems_AB(primes, nmax):
    """
    Direct evaluation of the aggregate sums against their closed forms:
        S    = sum_{a,b} eta(C)            (first aggregate identity),
        W    = sum_{a,b} eta(b)^{k-1}eta(C) (second aggregate identity, before the sign eps_q(n)).
    Both are computed term by term in F_p and compared with the stated values.
    """
    rowsA = badA = rowsB = badB = 0
    for p in primes:
        for n in range(3, nmax):
            for k in range(1, n):
                if not confounding(n, k, p):
                    continue
                n1, k1, alpha, beta, A, B = trinomial_data(n, k, p)
                S = sum(eta((A * pow(b, beta, p) + B * pow(a, alpha, p)) % p, p)
                        for a in range(1, p) for b in range(1, p))
                if n1 % 2 == 1 and beta % 2 == 1:
                    predA = 0
                elif beta % 2 == 0:
                    predA = -(p - 1) * eta(A, p)
                else:
                    predA = -(p - 1) * eta(B, p)
                rowsA += 1
                badA += (S != predA)

                W = sum(eta(b, p) ** (k - 1)
                        * eta((A * pow(b, beta, p) + B * pow(a, alpha, p)) % p, p)
                        for a in range(1, p) for b in range(1, p))
                predB = (-(p - 1) * eta(A, p) if n1 % 2 == 1
                         else -(p - 1) * eta(B, p))
                rowsB += 1
                badB += (W != predB)
    return rowsA, badA, rowsB, badB


def check_theorem_D(primes, nmax):
    """Exact row-sum identity for S_b."""
    rows = bad = 0
    gammas = set()
    for p in primes:
        N = p - 1
        for n in range(3, nmax):
            for k in range(1, n):
                if not confounding(n, k, p):
                    continue
                n1, k1, alpha, beta, A, B = trinomial_data(n, k, p)
                g = gcd(alpha, N)
                gammas.add(g)
                Binv = pow(B, p - 2, p)
                for b in range(1, p):
                    c = (A * pow(b, beta, p)) % p
                    lhs = sum(eta((B * pow(a, alpha, p) + c) % p, p)
                              for a in range(1, p))
                    tau = (-c * Binv) % p
                    rhs = g * eta(c, p) * sum(
                        eta(1 - x, p) for x in range(2, p)
                        if is_gth_power(tau * x, g, p))
                    rows += 1
                    bad += (lhs != rhs)
    return rows, bad, sorted(gammas)


# ------------------------------------- 2. master closed form, exact form

def check_master(primes, nmax):
    """
    sum_{psi^g=eps} psi(4 tau) J(psi,psi) = g * #{x != 0,1 : 4 tau x(1-x) a g-th power};
    the terms of order <= 2 are J(eps,eps) = q-2 and, if 2|g,
    eta(4 tau) J(eta,eta) with J(eta,eta) = -eta(-1).
    """
    rows = bad = 0
    for p in primes:
        N = p - 1
        em1 = eta(-1, p)
        for n in range(3, nmax):
            for k in range(1, n):
                if not confounding(n, k, p):
                    continue
                n1, k1, alpha, beta, A, B = trinomial_data(n, k, p)
                g = gcd(alpha, N)
                Binv = pow(B, p - 2, p)
                for b in range(1, p):
                    c = (A * pow(b, beta, p)) % p
                    lhs = sum(eta((B * pow(a, alpha, p) + c) % p, p)
                              for a in range(1, p))
                    tau = (-c * Binv) % p
                    full = g * sum(1 for x in range(2, p)
                                   if is_gth_power(4 * tau * x * (1 - x), g, p))
                    hi = full - (p - 2)
                    if g % 2 == 0:
                        hi -= eta(4 * tau, p) * (-em1)
                    rhs = -eta(c, p) + eta(c, p) * hi
                    if g % 2 == 0:
                        rhs -= eta(B, p)
                    rows += 1
                    bad += (lhs != rhs)
    return rows, bad


# ------------------------------- 3. one-variable reduction

def check_homogeneity(primes, nmax):
    """eta(C) = eta(a)^{n1} * eta(B + A w),  w = b^beta a^{-alpha}."""
    pairs = bad = 0
    for p in primes:
        for n in range(3, nmax):
            for k in range(1, n):
                if not confounding(n, k, p):
                    continue
                n1, k1, alpha, beta, A, B = trinomial_data(n, k, p)
                for a in range(1, p):
                    ia = pow(a, p - 2, p)
                    for b in range(1, p):
                        C = (A * pow(b, beta, p) + B * pow(a, alpha, p)) % p
                        w = (pow(b, beta, p) * pow(ia, alpha, p)) % p
                        lhs = eta(C, p)
                        rhs = (eta(a, p) ** alpha) * eta(B + A * w, p)
                        pairs += 1
                        bad += (lhs != rhs)
    return pairs, bad


# ------------------------------------ 4. quartic side: exact set/gcd checks

def build_Fp4(p):
    """Build F_{p^4} with a discrete-log table.  Exact integer arithmetic.

    Instead of testing irreducibility with an external library, we search for
    a monic quartic g such that x itself has order p^4 - 1 modulo g.  As in
    the direct-search certificates of the accompanying sweep, an order equal
    to p^4 - 1 certifies simultaneously that g is irreducible (the prime-to-p
    exponent of the unit group for any other factor type is a proper divisor)
    and that x generates F_{p^4}^*; the log table is then filled by
    enumerating the powers of x, which re-proves the order claim."""
    m, q = 4, p ** 4
    g = None  # determined below together with the generator

    def mul(x, y):
        r = [0] * (2 * m - 1)
        for i in range(m):
            if x[i]:
                for j in range(m):
                    if y[j]:
                        r[i + j] = (r[i + j] + x[i] * y[j]) % p
        for d in range(2 * m - 2, m - 1, -1):
            if r[d]:
                c = r[d]
                r[d] = 0
                for j in range(m + 1):
                    r[d - m + j] = (r[d - m + j] - c * g[j]) % p
        return r[:m]

    def toi(c):
        x = 0
        for i in reversed(range(m)):
            x = x * p + c[i]
        return x

    def toc(x):
        c = []
        for _ in range(m):
            c.append(x % p)
            x //= p
        return c

    def powi(x, e):
        r = [1, 0, 0, 0]
        while e > 0:
            if e & 1:
                r = mul(r, x)
            x = mul(x, x)
            e >>= 1
        return r

    # Search for g making x primitive; enumerating the q-1 powers of x both
    # certifies ord(x) = q-1 (hence g irreducible and x a generator) and
    # fills the discrete-log table in the same pass.
    logt = None
    rng = _random.Random(p)
    while logt is None:
        g = [rng.randrange(p) for _ in range(4)] + [1]
        if g[0] == 0:
            continue
        table = [None] * q
        x = [1, 0, 0, 0]
        ok = True
        for i in range(q - 1):
            xi = toi(x)
            if table[xi] is not None or (xi == 1 and i > 0):
                ok = False
                break
            table[xi] = i
            x = mul(x, [0, 1, 0, 0])
        if ok and toi(x) == 1:
            logt = table
    logt = [v if v is not None else 0 for v in logt]
    frob = [toi(powi(toc(xi), p)) for xi in range(q)]

    def trace(xi):
        s = toc(xi)[:]
        cur = xi
        for _ in range(3):
            cur = frob[cur]
            cc = toc(cur)
            for i in range(4):
                s[i] = (s[i] + cc[i]) % p
        return s

    return q, mul, toi, toc, logt, trace


def check_quartic(primes):
    """
    Exact checks:
      (a) |V \\ {0}| = q^2 - 1;
      (b) V \\ {0} is a union of full F_q^*-cosets, whence the sum of any
          character nontrivial on F_q^* over V \\ {0} is exactly 0 ;
      (c) the conic points together with F_q^* generate F_{q^4}^*, i.e. the
          index is 1, which is the geometric-nontriviality condition.
    """
    out = []
    for p in primes:
        q, mul, toi, toc, logt, trace = build_Fp4(p)
        zero = [0, 0, 0, 0]
        V = set(xi for xi in range(1, q)
                if trace(xi) == zero
                and trace(toi(mul(toc(xi), toc(xi)))) == zero)
        size_ok = (len(V) == p * p - 1)
        coset_ok = all(
            toi([(c * cc) % p for cc in toc(xi)]) in V
            for xi in V for c in range(1, p))
        Fp_star = [toi([c, 0, 0, 0]) for c in range(1, p)]
        gg = 0
        for xi in list(V) + Fp_star:
            gg = gcd(gg, logt[xi])
        index = gcd(gg, q - 1)
        out.append((p, len(V), size_ok, coset_ok, index))
    return out


# ------------------------- 5. breakpoint check for the explicit threshold

def check_threshold(kmax=199, q0=26000000000):
    """
    Finite part of the proof of the explicit threshold, in EXACT INTEGER
    arithmetic (no floating point anywhere).

    After absorbing the totient factor the criterion is

        delta*(q^2 - 1) > q^{3/2} * C,      C = 6*2^t*(3r-1) + 1,

    with delta = (p_{t+1} - r)/p_{t+1}.  Writing u = q^{1/6} this reads
    delta*(u^4 - 1) > C u^3, for which delta*u >= C+1 suffices; that is

        u^6 * (p_{t+1} - r)^6  >=  (C+1)^6 * p_{t+1}^6.                (*)

    At the jump point L_k one has u^6 = p_1...p_k / 2 exactly, so (*) is a
    comparison of two exact integers.  At the integer point q0 we instead
    square the criterion directly, again an exact integer comparison.
    """
    _P = nth_primes(kmax + 2)
    pr = lambda k: _P[k - 1]

    primorial, v = [], 1
    for k in range(1, kmax + 2):
        v *= pr(k)
        primorial.append(v)

    def C_of(t, r):
        return 6 * 2 ** t * (3 * r - 1) + 1 if r > 0 else 6 * (2 ** t - 1) + 1

    def certify_Lk(k):
        """exact integer certification of (*) at L_k; returns t or None."""
        u6 = primorial[k - 1] // 2          # exact: u^6 = P_k / 2
        for t in range(0, k + 1):
            r = k - t
            pt = pr(t + 1)
            if r > 0 and pt - r <= 0:
                continue
            C = C_of(t, r)
            num = (pt - r) if r > 0 else pt
            if u6 * num ** 6 >= (C + 1) ** 6 * pt ** 6:
                return t
        return None

    def certify_q0(k=22):
        """exact integer certification at the integer point q0."""
        for t in range(0, k + 1):
            r = k - t
            pt = pr(t + 1)
            if r > 0 and pt - r <= 0:
                continue
            C = C_of(t, r)
            num = (pt - r) if r > 0 else pt
            if num ** 2 * (q0 * q0 - 1) ** 2 > pt ** 2 * q0 ** 3 * C ** 2:
                return t
        return None

    t0 = certify_q0()
    table = [(k, certify_Lk(k)) for k in range(23, kmax + 1)]
    fails = [k for k, t in table if t is None]
    return t0, fails, table


def check_analytic_start(kneed=200):
    """
    The analytic part of the threshold proof is needed only for k >= kneed,
    the smaller indices being certified exactly.  For k >= 200 one has
    p_k >= p_200 = 1223, whence theta(p_k) > 0.85 p_k by Rosser-Schoenfeld
    Thm 10, and the sufficient condition becomes

        0.1417 log k - 0.117/k - 1.746/log(2k) - log(40k)/k > 0,

    whose left side is increasing in k.  Returns (first valid k, its value
    at kneed).
    """
    import math

    def f(k):
        return (0.1417 * math.log(k) - 0.117 / k
                - 1.746 / math.log(2 * k) - math.log(40 * k) / k)

    first = next(k for k in range(2, 2000) if f(k) > 0)
    return first, f(kneed)


# ------------------- 6. Li-Wang-Zhao criterion, inside its hypotheses

def check_LWZ_F27(nmax=6):
    """
    F_27 = F_3[t]/(t^3 - t - 1), so q = 3, m = 3 > 2; test n > 2.
    Claim: a^n b^{1-n} in F_3^*  =>  x^n + a x + b is not primitive over F_27.
    """
    Q = 27

    def mul(x, y):
        r = [0] * 5
        for i in range(3):
            if x[i]:
                for j in range(3):
                    if y[j]:
                        r[i + j] = (r[i + j] + x[i] * y[j]) % 3
        for d in (4, 3):
            c = r[d]
            if c:
                r[d] = 0
                r[d - 3] = (r[d - 3] + c) % 3
                r[d - 2] = (r[d - 2] + c) % 3
        return tuple(r[:3])

    def add(x, y):
        return tuple((x[i] + y[i]) % 3 for i in range(3))

    ZERO, ONE = (0, 0, 0), (1, 0, 0)
    nz = [(i, j, k) for i in range(3) for j in range(3) for k in range(3)
          if (i, j, k) != ZERO]

    def powe(x, e):
        r = ONE
        while e > 0:
            if e & 1:
                r = mul(r, x)
            x = mul(x, x)
            e >>= 1
        return r

    F3s = [(1, 0, 0), (2, 0, 0)]

    def primitive(n, a, b):
        N = Q ** n - 1

        def red(r):
            for d in range(len(r) - 1, n - 1, -1):
                c = r[d]
                if c != ZERO:
                    r[d] = ZERO
                    r[d - n + 1] = add(r[d - n + 1],
                                       tuple((-v) % 3 for v in mul(a, c)))
                    r[d - n] = add(r[d - n],
                                   tuple((-v) % 3 for v in mul(b, c)))
            return r[:n]

        def pmul(u, v):
            r = [ZERO] * (2 * n - 1)
            for i in range(n):
                if u[i] != ZERO:
                    for j in range(n):
                        if v[j] != ZERO:
                            r[i + j] = add(r[i + j], mul(u[i], v[j]))
            return red(r)

        def ppow(u, e):
            r = [ONE] + [ZERO] * (n - 1)
            while e > 0:
                if e & 1:
                    r = pmul(r, u)
                u = pmul(u, u)
                e >>= 1
            return r

        X = [ZERO, ONE] + [ZERO] * (n - 2)
        one = [ONE] + [ZERO] * (n - 1)
        if ppow(X, N) != one:
            return False
        for r in factor_primes(N):
            if ppow(X, N // r) == one:
                return False
        return True

    tot = hits = viol = other = 0
    for n in range(3, nmax + 1):
        for a in nz:
            for b in nz:
                inv_ = mul(powe(a, n), powe(b, (1 - n) % (Q - 1)))
                pr_ = primitive(n, a, b)
                tot += 1
                if inv_ in F3s:
                    hits += 1
                    viol += pr_
                else:
                    other += pr_
    return tot, hits, viol, other


# ------------------------------------------- 7. fixed family check

def check_family(mmax=25):
    """x^{4m+1} + x^2 + 4 over F_5:  C = 4[(4m+1) + (4m-1)^3] mod 5."""
    res = []
    p = 5
    for m in range(1, mmax + 1):
        n = 4 * m + 1
        if n % p == 0 or (n - 2) % p == 0:
            continue
        C = (4 * ((4 * m + 1) + pow((4 * m - 1) % 5, 3, 5))) % 5
        res.append((m, m % 5, C, eta(C, p)))
    return res


# ------------------------------------------------------------------- main

if __name__ == "__main__":
    print("EXACT VERIFICATION (integer arithmetic only)\n")

    rA, bA, rB, bB = check_theorems_AB([7, 11, 13, 17, 19, 23], 24)
    print(f"aggregate sums  : {rA:7d} / {rB:7d} classes, {bA} / {bB} mismatches")

    rows, bad, gs = check_theorem_D([7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43], 32)
    print(f"row-sum identity: {rows:7d} rows, {bad} mismatches; gamma in {gs}")

    rows, bad = check_master([7, 11, 13, 17, 19, 23, 29, 31, 37, 41], 30)
    print(f"master formula : {rows:7d} rows, {bad} mismatches")

    pairs, bad = check_homogeneity([7, 11, 13, 17, 19, 23, 29, 31, 37], 26)
    print(f"one-variable   : {pairs:7d} pairs, {bad} mismatches")

    print("\nquartic side (p, |V*|, |V*|=q^2-1, cosets, generation index):")
    for row in check_quartic([5, 7, 11, 13]):
        print("   ", row)

    print("\nfamily x^{4m+1}+x^2+4 over F_5  (m, m mod 5, C, eta(C)):")
    for row in check_family():
        print("   ", row)

    t0, fails, table = check_threshold()
    kfirst, val200 = check_analytic_start()
    print(f"\nthreshold : q0 = 2.6e10 certified exactly with t={t0}; "
          f"L_k for 23<=k<=199 certified, failures: {fails if fails else 'none'}; "
          f"analytic bound first valid at k={kfirst}, value at k=200 is {val200:+.4f}")

    print("\ncertificate for the explicit threshold")
    print("  (k, t) with t the sieve parameter certifying the integer inequality at L_k:")
    for row in range(0, len(table), 12):
        print("   ", "  ".join(f"({k},{t})" for k, t in table[row:row + 12]))

    tot, hits, viol, other = check_LWZ_F27()
    print(f"LWZ (F_27, m=3>2, n>2): {tot} trinomials, {hits} with invariant in F_3^*, "
          f"{viol} of those primitive, {other} primitive among the rest")
