#!/usr/bin/env bash
set -euo pipefail

echo "=== Benchmark: Wavefront Alignment vs. Quadratic DP ==="

python3 -c '
import time, random
from verify_logic import gotoh_dp, wfa_solve

bases = ["A", "C", "G", "T"]
for length in [200, 500, 1000]:
    t = "".join(random.choice(bases) for _ in range(length))
    # 2% divergence
    q = list(t)
    for i in range(0, len(q), 50):
        q[i] = random.choice([b for b in bases if b != q[i]])
    q = "".join(q)

    # Benchmark WFA
    t0 = time.perf_counter()
    wfa_s = wfa_solve(q, t)
    t_wfa = time.perf_counter() - t0

    # Benchmark DP
    t0 = time.perf_counter()
    dp_s = gotoh_dp(q, t)
    t_dp = time.perf_counter() - t0

    speedup = t_dp / t_wfa if t_wfa > 0 else 0
    print(f"Length {length:4d} bp: WFA = {t_wfa*1000:6.2f} ms, DP = {t_dp*1000:6.2f} ms | Speedup: {speedup:5.1f}x (Scores: {wfa_s} == {dp_s})")
'

echo "Benchmark completed."
