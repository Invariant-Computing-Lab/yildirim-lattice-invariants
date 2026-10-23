#!/usr/bin/env python3
"""Reproducibility driver for the invariant-geometry synthetic experiments.

The script is intentionally self-contained. Exact integer/Fraction arithmetic is used
for discrete claims. SciPy is used only as an independent numerical convex-hull check.

Outputs are written under ../results relative to this file unless --out is supplied.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import random
import sys
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

SEED_MASTER = 20260812

# ----------------------------- exact linear algebra -----------------------------

def det3cols(cols: Sequence[Sequence[int | Fraction]]) -> int | Fraction:
    a, b, c = cols
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    )


def gcd_nonzero(vals: Iterable[int]) -> int:
    g = 0
    for v in vals:
        g = gcd(g, abs(int(v)))
    return g


def primitive(v: Sequence[int]) -> Tuple[int, int, int]:
    g = gcd_nonzero(v)
    if g == 0:
        return tuple(v)  # type: ignore
    return tuple(int(x // g) for x in v)  # type: ignore


def generic_rank3(G: Sequence[Sequence[int | Fraction]]) -> bool:
    return all(det3cols([G[i], G[j], G[k]]) != 0 for i, j, k in itertools.combinations(range(len(G)), 3))


def rand_generic_primitive_G(m: int, bound: int, rng: random.Random) -> List[Tuple[int, int, int]]:
    while True:
        G: List[Tuple[int, int, int]] = []
        seen = set()
        while len(G) < m:
            v = tuple(rng.randint(-bound, bound) for _ in range(3))
            if v == (0, 0, 0):
                continue
            v = primitive(v)
            w = v
            for x in w:
                if x != 0:
                    if x < 0:
                        w = tuple(-z for z in w)
                    break
            if w in seen:
                continue
            seen.add(w)
            G.append(v)
        if generic_rank3(G):
            return G


# ----------------------------- chirotope quotient ------------------------------

def parity(vals: Sequence[int]) -> int:
    inv = 0
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            inv += int(vals[i] > vals[j])
    return inv & 1


def chir_bits(G: Sequence[Sequence[int | Fraction]], rank: int = 3) -> int:
    bits = 0
    for bi, B in enumerate(itertools.combinations(range(len(G)), rank)):
        d = det3cols([G[i] for i in B])
        if d == 0:
            raise ValueError("Chirotope bit encoding requires a uniform/generic realization")
        if d > 0:
            bits |= 1 << bi
    return bits


def precompute_chiro_quotient(m: int, rank: int = 3):
    bases = list(itertools.combinations(range(m), rank))
    bidx = {B: i for i, B in enumerate(bases)}
    full = (1 << len(bases)) - 1
    reor_masks = []
    for e_bits in range(1 << m):
        mask = 0
        for bi, B in enumerate(bases):
            if sum((e_bits >> i) & 1 for i in B) & 1:
                mask |= 1 << bi
        reor_masks.append(mask)
    perms = []
    for p in itertools.permutations(range(m)):
        src = []
        flipmask = 0
        for bi, B in enumerate(bases):
            vals = [p[i] for i in B]
            src.append(bidx[tuple(sorted(vals))])
            if parity(vals):
                flipmask |= 1 << bi
        perms.append((src, flipmask))
    return bases, full, tuple(reor_masks), tuple(perms)


def permute_bits(bits: int, src: Sequence[int], flipmask: int) -> int:
    out = flipmask
    for dest, si in enumerate(src):
        if (bits >> si) & 1:
            out ^= 1 << dest
    return out


_PRECOMP6 = precompute_chiro_quotient(6, 3)
_PRECOMP7 = precompute_chiro_quotient(7, 3)

# Fast m=6 quotient table: reorientation/global-sign then relabeling.
_BASES6, _FULL6, _REOR6, _PERMS6 = _PRECOMP6
_arr = np.arange(1 << len(_BASES6), dtype=np.uint32)
_reor_best6 = np.full_like(_arr, np.uint32(_FULL6))
for _rm in _REOR6:
    _x = _arr ^ np.uint32(_rm)
    _reor_best6 = np.minimum(_reor_best6, np.minimum(_x, _x ^ np.uint32(_FULL6)))

_perm_low6 = np.zeros((len(_PERMS6), 1024), dtype=np.uint32)
_perm_high6 = np.zeros((len(_PERMS6), 1024), dtype=np.uint32)
_flip6 = np.zeros(len(_PERMS6), dtype=np.uint32)
for pi, (src, flipmask) in enumerate(_PERMS6):
    _flip6[pi] = flipmask
    contrib = [0] * len(_BASES6)
    for dest, si in enumerate(src):
        contrib[si] |= 1 << dest
    for x in range(1024):
        out = 0
        xx = x
        while xx:
            lsb = xx & -xx
            idx = lsb.bit_length() - 1
            out ^= contrib[idx]
            xx -= lsb
        _perm_low6[pi, x] = out
        out = 0
        xx = x
        while xx:
            lsb = xx & -xx
            idx = lsb.bit_length() - 1
            out ^= contrib[10 + idx]
            xx -= lsb
        _perm_high6[pi, x] = out


def canon_chiro(bits: int, m: int) -> int:
    if m == 6:
        lo = bits & 1023
        hi = (bits >> 10) & 1023
        vals = _perm_low6[:, lo] ^ _perm_high6[:, hi] ^ _flip6
        return int(_reor_best6[vals].min())
    bases, full, reors, perms = _PRECOMP7 if m == 7 else precompute_chiro_quotient(m, 3)
    best = None
    for src, flipmask in perms:
        pb = permute_bits(bits, src, flipmask)
        for rm in reors:
            x = pb ^ rm
            x = min(x, x ^ full)
            if best is None or x < best:
                best = x
    assert best is not None
    return best


# --------------------------- arithmetic multiplicities --------------------------

def subset_rank3_multiplicity(G: Sequence[Sequence[int]], S: Sequence[int], symmetric: bool = True) -> int:
    # Symmetric P=sum[-g_i,g_i] is an integer translate of Z(2g_i), so use v_i=2g_i.
    cols = [tuple((2 * int(x) if symmetric else int(x)) for x in G[i]) for i in S]
    j = len(S)
    if j == 0:
        return 1
    if j == 1:
        return gcd_nonzero(cols[0])
    if j == 2:
        a, b = cols
        return gcd_nonzero([
            a[0] * b[1] - a[1] * b[0],
            a[0] * b[2] - a[2] * b[0],
            a[1] * b[2] - a[2] * b[1],
        ])
    # For rank 3, multiplicity is gcd of all 3x3 minors of the sublist.
    ds = [det3cols([cols[a], cols[b], cols[c]]) for a, b, c in itertools.combinations(range(j), 3)]
    return gcd_nonzero(ds)


def arithmetic_profile(G: Sequence[Sequence[int]], include_dependent: bool = True) -> Tuple[int, ...]:
    m = len(G)
    vals = []
    max_size = m if include_dependent else 3
    for j in range(1, max_size + 1):
        for S in itertools.combinations(range(m), j):
            vals.append(subset_rank3_multiplicity(G, S, symmetric=True))
    return tuple(vals)


def profile_subsets(m: int, include_dependent: bool = True):
    subs = []
    max_size = m if include_dependent else 3
    for j in range(1, max_size + 1):
        subs.extend(itertools.combinations(range(m), j))
    return tuple(subs)


def precompute_profile_perm_maps(m: int, include_dependent: bool = True):
    subs = profile_subsets(m, include_dependent)
    idx = {S: i for i, S in enumerate(subs)}
    maps = []
    for p in itertools.permutations(range(m)):
        maps.append(tuple(idx[tuple(sorted(p[i] for i in S))] for S in subs))
    return subs, tuple(maps)


_SUBS6_ALL, _PROFILE_PERMS6_ALL = precompute_profile_perm_maps(6, True)


def canon_arithmetic_profile(G: Sequence[Sequence[int]]) -> Tuple[int, ...]:
    if len(G) != 6:
        subs, maps = precompute_profile_perm_maps(len(G), True)
    else:
        maps = _PROFILE_PERMS6_ALL
    prof = arithmetic_profile(G, include_dependent=True)
    best = None
    for mp in maps:
        v = tuple(prof[i] for i in mp)
        if best is None or v < best:
            best = v
    assert best is not None
    return best


def ehrhart_coeffs_symmetric_zonotope(G: Sequence[Sequence[int]]) -> Tuple[int, int, int, int]:
    # Stanley's lattice-zonotope formula, for translated generators 2g_i.
    c = [1, 0, 0, 0]
    m = len(G)
    for j in (1, 2, 3):
        for S in itertools.combinations(range(m), j):
            c[j] += subset_rank3_multiplicity(G, S, symmetric=True)
    return tuple(c)  # type: ignore


def hash_profile(x: Sequence[int]) -> str:
    return hashlib.sha256(json.dumps(list(x), separators=(",", ":")).encode()).hexdigest()[:16]


# ------------------------------- group actions ---------------------------------

def det3mat(A: Sequence[Sequence[int]]) -> int:
    cols = [tuple(A[i][j] for i in range(3)) for j in range(3)]
    return int(det3cols(cols))


def matvec(A: Sequence[Sequence[int]], v: Sequence[int]) -> Tuple[int, int, int]:
    return tuple(sum(int(A[i][j]) * int(v[j]) for j in range(3)) for i in range(3))  # type: ignore


def rand_GL3(rng: random.Random, bound: int = 4):
    while True:
        A = [[rng.randint(-bound, bound) for _ in range(3)] for __ in range(3)]
        if det3mat(A) != 0:
            return A


def rand_unimodular3(rng: random.Random, steps: int = 8):
    A = np.eye(3, dtype=int)
    for _ in range(steps):
        op = rng.choice(["swap", "sign", "add"])
        if op == "swap":
            i, j = rng.sample(range(3), 2)
            A[[i, j]] = A[[j, i]]
        elif op == "sign":
            i = rng.randrange(3)
            A[i] *= -1
        else:
            i, j = rng.sample(range(3), 2)
            k = rng.choice([-3, -2, -1, 1, 2, 3])
            A[i] += k * A[j]
    return A.tolist()


def transform_structural(G: Sequence[Sequence[int]], rng: random.Random):
    A = rand_GL3(rng)
    p = list(range(len(G)))
    rng.shuffle(p)
    eps = [rng.choice([-1, 1]) for _ in G]
    return [matvec(A, tuple(eps[new] * int(x) for x in G[old])) for new, old in enumerate(p)]


def transform_lattice(G: Sequence[Sequence[int]], rng: random.Random):
    U = rand_unimodular3(rng)
    p = list(range(len(G)))
    rng.shuffle(p)
    eps = [rng.choice([-1, 1]) for _ in G]
    return [matvec(U, tuple(eps[new] * int(x) for x in G[old])) for new, old in enumerate(p)]


# --------------------------------- statistics ----------------------------------

def entropy(labels: Sequence[object]) -> float:
    N = len(labels)
    cnt = collections.Counter(labels)
    return -sum((n / N) * math.log2(n / N) for n in cnt.values())


def conditional_entropy(C: Sequence[object], I: Sequence[object]) -> float:
    N = len(C)
    groups: Dict[object, List[object]] = collections.defaultdict(list)
    for c, i in zip(C, I):
        groups[i].append(c)
    H = 0.0
    for vals in groups.values():
        n = len(vals)
        cnt = collections.Counter(vals)
        h = -sum((k / n) * math.log2(k / n) for k in cnt.values())
        H += (n / N) * h
    return H


def collision_rate_given_different(C: Sequence[object], I: Sequence[object]):
    N = len(C)
    total = N * (N - 1) // 2 - sum(n * (n - 1) // 2 for n in collections.Counter(C).values())
    bad = 0
    groups: Dict[object, List[object]] = collections.defaultdict(list)
    for c, i in zip(C, I):
        groups[i].append(c)
    for vals in groups.values():
        n = len(vals)
        allpairs = n * (n - 1) // 2
        same_class = sum(k * (k - 1) // 2 for k in collections.Counter(vals).values())
        bad += allpairs - same_class
    return bad / total if total else 0.0, bad, total


def wilson(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / den
    return (ctr - half, ctr + half)


# -------------------------------- experiments ----------------------------------

def experiment_0_hull(out: Path):
    from scipy.spatial import ConvexHull
    rows = []
    rng = np.random.default_rng(SEED_MASTER)
    for n in range(2, 7):
        pts = np.array(list(itertools.product([-1.0, 0.0, 1.0], repeat=n)), dtype=float)
        hull = ConvexHull(pts)
        # scipy returns extreme point indices
        nv = len(hull.vertices)
        max_support_error = 0.0
        for _ in range(500):
            u = rng.normal(size=n)
            h_grid = float(np.max(pts @ u))
            h_cube = float(np.sum(np.abs(u)))
            max_support_error = max(max_support_error, abs(h_grid - h_cube))
        rows.append({
            "n": n,
            "precursor_points": 3 ** n,
            "hull_vertices_numeric": nv,
            "correct_vertices_2^n": 2 ** n,
            "legacy_vertices_claim": 2 ** (n + 1) if n >= 2 else "NA",
            "max_support_abs_error_500_dirs": max_support_error,
        })
    write_csv(out / "tables/exp0_convex_hull.csv", rows)
    return rows


def experiment_1_exact_identities(out: Path):
    rows = []
    for n in range(1, 51):
        f = [2 ** (n - k) * math.comb(n, k) for k in range(n + 1)]
        total = sum(f)
        boundary_euler = sum(((-1) ** k) * f[k] for k in range(n))
        expected_euler = 1 + (-1) ** (n - 1)
        recurrence_ok = True
        if n >= 2:
            prev = [2 ** (n - 1 - k) * math.comb(n - 1, k) for k in range(n)]
            for k in range(n):
                left = f[k]
                right = 2 * (prev[k] if k < len(prev) else 0) + (prev[k - 1] if k - 1 >= 0 else 0)
                if left != right:
                    recurrence_ok = False
                    break
        rows.append({
            "n": n,
            "sum_nonempty_faces_including_polytope": total,
            "expected_3^n": 3 ** n,
            "ternary_identity_ok": total == 3 ** n,
            "boundary_euler": boundary_euler,
            "expected_boundary_euler": expected_euler,
            "euler_ok": boundary_euler == expected_euler,
            "recurrence_ok": recurrence_ok,
        })
    write_csv(out / "tables/exp1_exact_identities.csv", rows)
    return rows


def experiment_2_quotient(out: Path):
    rows = []
    for m, N, seed in [(6, 1000, SEED_MASTER + 3), (7, 100, SEED_MASTER + 4)]:
        rng = random.Random(seed)
        structural_pass = 0
        raw_pass = 0
        for _ in range(N):
            G = rand_generic_primitive_G(m, 5, rng)
            Gp = transform_structural(G, rng)
            raw_pass += int(chir_bits(G) == chir_bits(Gp))
            structural_pass += int(canon_chiro(chir_bits(G), m) == canon_chiro(chir_bits(Gp), m))
        lo, hi = wilson(structural_pass, N)
        rows.append({
            "regime": "structural_GL(R)_signed-permutation",
            "m": m,
            "N": N,
            "quotient_invariant_pass": structural_pass,
            "quotient_invariant_rate": structural_pass / N,
            "wilson95_low": lo,
            "wilson95_high": hi,
            "raw_chirotope_same": raw_pass,
            "raw_rate": raw_pass / N,
        })
    # lattice arithmetic profile under GL_Z
    rng = random.Random(SEED_MASTER + 5)
    N = 300
    arith_pass = 0
    raw_pass = 0
    for _ in range(N):
        G = rand_generic_primitive_G(6, 4, rng)
        Gp = transform_lattice(G, rng)
        raw_pass += int(arithmetic_profile(G) == arithmetic_profile(Gp))
        arith_pass += int(canon_arithmetic_profile(G) == canon_arithmetic_profile(Gp))
    lo, hi = wilson(arith_pass, N)
    rows.append({
        "regime": "lattice_GL(Z)_signed-permutation",
        "m": 6,
        "N": N,
        "quotient_invariant_pass": arith_pass,
        "quotient_invariant_rate": arith_pass / N,
        "wilson95_low": lo,
        "wilson95_high": hi,
        "raw_chirotope_same": raw_pass,
        "raw_rate": raw_pass / N,
    })
    write_csv(out / "tables/exp2_quotient_preservation.csv", rows)
    return rows


def experiment_3_random_atlas(out: Path, N: int = 2000):
    rng = random.Random(SEED_MASTER + 7)
    rows = []
    target = []
    feats = collections.defaultdict(list)
    om_counts = collections.Counter()
    for idx in range(N):
        G = rand_generic_primitive_G(6, 4, rng)
        om = canon_chiro(chir_bits(G), 6)
        arith = canon_arithmetic_profile(G)
        ehr = ehrhart_coeffs_symmetric_zonotope(G)
        target_class = (om, arith)
        target.append(target_class)
        om_counts[om] += 1
        feat_values = {
            "rank_m": (3, 6),
            "f_vector": (32, 60, 30),
            "oriented_matroid": om,
            "arithmetic_profile": arith,
            "ehrhart": ehr,
            "volume": ehr[3],
            "OM_plus_Ehrhart": (om, ehr),
            "OM_plus_Arithmetic": target_class,
        }
        for name, val in feat_values.items():
            feats[name].append(val)
        rows.append({
            "id": idx,
            "om_class": om,
            "arith_profile_hash": hash_profile(arith),
            "ehrhart_c0": ehr[0],
            "ehrhart_c1": ehr[1],
            "ehrhart_c2": ehr[2],
            "ehrhart_c3_volume": ehr[3],
            "generator_matrix": json.dumps(G, separators=(",", ":")),
        })
    write_csv(out / "tables/exp3_random_atlas.csv", rows)

    Hc = entropy(target)
    metric_rows = []
    for name, values in feats.items():
        Hcond = conditional_entropy(target, values)
        R = 1.0 - Hcond / Hc if Hc else 1.0
        coll, bad, total = collision_rate_given_different(target, values)
        metric_rows.append({
            "invariant": name,
            "unique_signatures": len(set(values)),
            "reconstructibility_R": R,
            "conditional_entropy_bits": Hcond,
            "collision_rate_given_different_target": coll,
            "collision_pairs": bad,
            "different_target_pairs": total,
        })
    write_csv(out / "tables/exp3_reconstructibility_metrics.csv", metric_rows)
    return rows, metric_rows, dict(om_counts)


def box_ehrhart(a: Sequence[int]):
    # prod_i (2 a_i k + 1), coefficients increasing order.
    coeff = [1]
    for ai in a:
        nxt = [0] * (len(coeff) + 1)
        for j, c in enumerate(coeff):
            nxt[j] += c
            nxt[j + 1] += 2 * ai * c
        coeff = nxt
    return tuple(coeff)


def experiment_4_boxes(out: Path):
    all_rows = []
    summaries = []
    for n, maxa in [(3, 10), (4, 7)]:
        tuples = list(itertools.combinations_with_replacement(range(1, maxa + 1), n))
        ehr = [box_ehrhart(a) for a in tuples]
        vol = [e[-1] for e in ehr]
        rows = []
        for a, e, v in zip(tuples, ehr, vol):
            rows.append({"n": n, "a": json.dumps(a), "volume": v, "ehrhart_coefficients": json.dumps(e)})
        all_rows.extend(rows)
        summaries.append({
            "n": n,
            "max_axis_parameter": maxa,
            "N_boxes": len(tuples),
            "unpolarized_toric_signatures": 1,
            "unique_volumes": len(set(vol)),
            "unique_ehrhart": len(set(ehr)),
            "ehrhart_uniqueness_rate": len(set(ehr)) / len(tuples),
        })
    write_csv(out / "tables/exp4_boxes.csv", all_rows)
    write_csv(out / "tables/exp4_boxes_summary.csv", summaries)
    return summaries


def experiment_5_adversarial(out: Path):
    # Exact examples found by deterministic search; retained explicitly for regression tests.
    G_A = [(1, -1, 0), (-4, 4, -3), (0, -1, 2), (-1, -4, 2), (-1, 0, -2), (-2, 3, -3)]
    G_B = [(2, -1, 2), (4, -4, -1), (-3, 4, 0), (-1, 2, -2), (0, 1, -1), (-3, 2, 4)]
    omA, omB = canon_chiro(chir_bits(G_A), 6), canon_chiro(chir_bits(G_B), 6)
    ehrA, ehrB = ehrhart_coeffs_symmetric_zonotope(G_A), ehrhart_coeffs_symmetric_zonotope(G_B)

    G_scale = [tuple(2 * x for x in G_A[0])] + G_A[1:]
    omS = canon_chiro(chir_bits(G_scale), 6)
    ehrS = ehrhart_coeffs_symmetric_zonotope(G_scale)

    # Box volume collision.
    box1, box2 = (1, 1, 6), (1, 2, 3)
    ebox1, ebox2 = box_ehrhart(box1), box_ehrhart(box2)

    # Exact Gram-spectrum counterexample under a unimodular shear.
    # G=I, U=[[1,1,0],[0,1,0],[0,0,1]], Gram eigenvalues are 1, (3±sqrt(5))/2.
    data = {
        "same_ehrhart_different_structural_class": {
            "G_A": G_A,
            "G_B": G_B,
            "OM_A": omA,
            "OM_B": omB,
            "ehrhart_A": ehrA,
            "ehrhart_B": ehrB,
            "verified": bool(ehrA == ehrB and omA != omB),
        },
        "same_structural_class_different_lattice_counting": {
            "G_A": G_A,
            "G_scaled_one_generator": G_scale,
            "OM_A": omA,
            "OM_scaled": omS,
            "ehrhart_A": ehrA,
            "ehrhart_scaled": ehrS,
            "verified": bool(omA == omS and ehrA != ehrS),
        },
        "same_volume_different_ehrhart_boxes": {
            "a": box1,
            "b": box2,
            "ehrhart_a": ebox1,
            "ehrhart_b": ebox2,
            "same_volume": ebox1[-1] == ebox2[-1],
            "different_ehrhart": ebox1 != ebox2,
        },
        "gram_spectrum_not_GLZ_invariant": {
            "G": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "U": [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
            "spectrum_G": [1, 1, 1],
            "spectrum_UG": ["(3-sqrt(5))/2", 1, "(3+sqrt(5))/2"],
            "det_U": 1,
        },
    }
    write_json(out / "tables/exp5_adversarial_cases.json", data)
    return data


def experiment_6_wall(out: Path):
    G5 = [(1, 1, -2), (5, -1, 2), (1, -2, -2), (4, -1, -4), (-5, -2, -4)]
    a = (11, -1, 2)  # g0 + 2 g1
    v = (-3, -3, 4)

    def make(t: Fraction):
        return [tuple(Fraction(x, 1) for x in g) for g in G5] + [
            tuple(Fraction(a[i], 1) + t * Fraction(v[i], 1) for i in range(3))
        ]

    roots = []
    for i, j in itertools.combinations(range(5), 2):
        A = det3cols([G5[i], G5[j], a])
        B = det3cols([G5[i], G5[j], v])
        if B != 0:
            roots.append({"basis": (i, j, 5), "A": int(A), "B": int(B), "root": Fraction(-A, B)})
    roots_sorted = sorted(roots, key=lambda d: float(d["root"]))

    sample_ts = [Fraction(-1, 20), Fraction(-1, 100), Fraction(-1, 1000), Fraction(0, 1), Fraction(1, 1000), Fraction(1, 100), Fraction(1, 20)]
    rows = []
    for t in sample_ts:
        G = make(t)
        zeros = [B for B in itertools.combinations(range(6), 3) if det3cols([G[i] for i in B]) == 0]
        cls = None if zeros else canon_chiro(chir_bits(G), 6)
        target_det = det3cols([G[0], G[1], G[5]])
        rows.append({
            "t": str(t),
            "t_float": float(t),
            "target_minor_det": str(target_det),
            "zero_minor_count": len(zeros),
            "zero_minors": json.dumps(zeros),
            "OM_class": cls if cls is not None else "undefined_on_discriminant",
        })
    write_csv(out / "tables/exp6_wall_crossing.csv", rows)
    write_json(out / "tables/exp6_wall_roots.json", {
        "G5": G5,
        "a": a,
        "v": v,
        "target_minor_formula": "det(g0,g1,g5(t)) = 12 t",
        "roots": [{**r, "root": str(r["root"]), "root_float": float(r["root"])} for r in roots_sorted],
    })
    return rows, roots_sorted


def experiment_7_legacy_table(out: Path):
    # Values as printed in the legacy source; only complete entries are checked.
    legacy = [
        {"n": 2, "f": [8, 8]},
        {"n": 3, "f": [16, 24, 10]},
        {"n": 4, "f": [32, 56, 34, 12]},
    ]
    rows = []
    for item in legacy:
        n = item["n"]
        f = item["f"]
        e = sum(((-1) ** k) * f[k] for k in range(n))
        expected = 1 + (-1) ** (n - 1)
        correct = [2 ** (n - k) * math.comb(n, k) for k in range(n)]
        rows.append({
            "n": n,
            "legacy_f_boundary": json.dumps(f),
            "legacy_euler": e,
            "expected_euler": expected,
            "legacy_euler_pass": e == expected,
            "correct_cube_boundary_f": json.dumps(correct),
        })
    write_csv(out / "tables/exp7_legacy_table_consistency.csv", rows)
    return rows


# ---------------------------------- figures ------------------------------------

def make_figures(out: Path, metric_rows, box_summary, wall_rows):
    import matplotlib.pyplot as plt

    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # Reconstructibility
    names = [r["invariant"] for r in metric_rows]
    vals = [r["reconstructibility_R"] for r in metric_rows]
    plt.figure(figsize=(9, 4.8))
    plt.bar(range(len(names)), vals)
    plt.xticks(range(len(names)), names, rotation=35, ha="right")
    plt.ylim(0, 1.05)
    plt.ylabel("Reconstructibility R")
    plt.title("Invariant reconstructibility in the synthetic rank-3, m=6 atlas")
    plt.tight_layout()
    plt.savefig(figdir / "figure_reconstructibility.png", dpi=220)
    plt.close()

    # Collision rate (log-ish display with floor)
    vals = [max(r["collision_rate_given_different_target"], 1e-6) for r in metric_rows]
    plt.figure(figsize=(9, 4.8))
    plt.bar(range(len(names)), vals)
    plt.yscale("log")
    plt.xticks(range(len(names)), names, rotation=35, ha="right")
    plt.ylabel("Collision rate | different target class")
    plt.title("Adversarial relevance of coarse versus layered invariants")
    plt.tight_layout()
    plt.savefig(figdir / "figure_collision_rates.png", dpi=220)
    plt.close()

    # Box uniqueness
    labels = [f"n={r['n']}" for r in box_summary]
    volume = [r["unique_volumes"] / r["N_boxes"] for r in box_summary]
    ehr = [r["unique_ehrhart"] / r["N_boxes"] for r in box_summary]
    x = np.arange(len(labels))
    w = 0.35
    plt.figure(figsize=(6.5, 4.6))
    plt.bar(x - w / 2, volume, width=w, label="Volume")
    plt.bar(x + w / 2, ehr, width=w, label="Ehrhart")
    plt.xticks(x, labels)
    plt.ylim(0, 1.05)
    plt.ylabel("Unique-signature fraction")
    plt.title("Lattice boxes: volume versus Ehrhart discrimination")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figdir / "figure_box_discrimination.png", dpi=220)
    plt.close()

    # Wall crossing
    xs = [r["t_float"] for r in wall_rows if isinstance(r["OM_class"], int)]
    ys = [r["OM_class"] for r in wall_rows if isinstance(r["OM_class"], int)]
    plt.figure(figsize=(6.5, 4.6))
    plt.step(xs, ys, where="mid")
    plt.axvline(0.0, linestyle="--")
    plt.xlabel("t")
    plt.ylabel("Canonical structural class code")
    plt.title("Exact discriminant-wall transition")
    plt.tight_layout()
    plt.savefig(figdir / "figure_wall_crossing.png", dpi=220)
    plt.close()


# ----------------------------------- I/O ---------------------------------------

def write_csv(path: Path, rows: Sequence[Dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = []
    for row in rows:
        for k in row:
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    ap.add_argument("--atlas-n", type=int, default=2000)
    args = ap.parse_args()
    out = args.out
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    res = {}
    res["exp0_hull"] = experiment_0_hull(out)
    res["exp1_exact"] = experiment_1_exact_identities(out)
    res["exp2_quotient"] = experiment_2_quotient(out)
    _, metric_rows, om_counts = experiment_3_random_atlas(out, N=args.atlas_n)
    res["exp3_metrics"] = metric_rows
    res["exp3_om_counts"] = om_counts
    box_summary = experiment_4_boxes(out)
    res["exp4_boxes"] = box_summary
    res["exp5_adversarial"] = experiment_5_adversarial(out)
    wall_rows, wall_roots = experiment_6_wall(out)
    res["exp6_wall"] = wall_rows
    res["exp7_legacy"] = experiment_7_legacy_table(out)

    make_figures(out, metric_rows, box_summary, wall_rows)

    session = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "master_seed": SEED_MASTER,
        "atlas_n": args.atlas_n,
    }
    try:
        import scipy
        session["scipy"] = scipy.__version__
    except Exception:
        pass
    try:
        import matplotlib
        session["matplotlib"] = matplotlib.__version__
    except Exception:
        pass
    write_json(out / "session_info.json", session)
    write_json(out / "results_summary.json", res)
    print(json.dumps({"status": "ok", "out": str(out), "session": session}, indent=2))


if __name__ == "__main__":
    main()
