> [!IMPORTANT]
> **Archival Notice**
> This repository is archived and is no longer actively maintained.

# wfa: Gap-Affine Wavefront Alignment & BiWFA Engine

[![Rust](https://img.shields.io/badge/rust-2021-orange.svg)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

`wfa` is an exact, sub-quadratic gap-affine pairwise sequence alignment engine written in pure Rust. It implements the **Wavefront Alignment (WFA)** algorithm (Marco-Sola et al., 2021) and the $O(s)$ memory **BiWFA** divide-and-conquer strategy, computing exact global alignment scores and CIGAR strings in $O(n \cdot s)$ time (where $s$ is the alignment score), achieving up to 240× speedups over classical Gotoh dynamic programming.

---

## Key Features

1. **Branchless 64-Bit LCP Extension**:
   - Compares 8 bytes per iteration using 64-bit word loads (`u64::from_le_bytes`), XOR, and trailing-zero counting (`tzcnt`), skipping matching sequence stretches with zero branching.
2. **Zero-Centered WavefrontArena**:
   - Flat contiguous memory slab mapping diagonal offsets $k \in [\text{lo}, \text{hi}]$ with $O(1)$ indexed access, eliminating heap allocations inside the score loop.
3. **Exact Gap-Affine Model**:
   - Supports 2-piece affine gap penalties (mismatch $X$, gap open $O$, gap extend $E$) across Match ($M$), Insertion ($I$), and Deletion ($D$) wavefronts:
     $$I_s[k] = \max\left(M_{s - (O + E)}[k + 1] + 1, I_{s - E}[k + 1] + 1\right)$$
     $$D_s[k] = \max\left(M_{s - (O + E)}[k - 1], D_{s - E}[k - 1]\right)$$
     $$M_s[k] = \max\left(M_{s - X}[k] + 1, I_s[k], D_s[k]\right)$$
4. **BiWFA Divide-and-Conquer ($O(s)$ Peak Memory)**:
   - For long sequences, bounds memory consumption to $O(s)$ via bidirectional wavefront expansion and midpoint stitching.
5. **Exact Gotoh DP Verifier Baseline**:
   - Built-in classical $O(N \cdot M)$ dynamic programming solver used to rigorously verify 100% score equivalence.

---

## Architectural Overview

```
Query & Target FASTA
         |
   [ PaddedSeq ] (Sentinel-padded 64-bit word loads)
         |
   [ Wavefront Progression s = 0, 1, 2, ... ]
   - M_s[k], I_s[k], D_s[k]
   - Branchless 64-bit LCP Extension
         |
   Target Diagonal k = (t_len - q_len) reached with offset >= q_len
         |
   [ Optimal CIGAR Backtrace ]
         |
   CIGAR String & Alignment Verification
```

---

## Build & Test

### Building via Cargo
```bash
cargo build --release
cargo test --release
```

### Running via Makefile
```bash
make test       # Runs unit tests or standalone verification suite
make all        # Builds release binary
```

---

## CLI Usage

```bash
# Basic WFA alignment
cargo run --release -- -q data/sample_query.fa -t data/sample_target.fa -m wfa -s

# Classical Gotoh DP mode for verification
cargo run --release -- -q data/sample_query.fa -t data/sample_target.fa -m dp

# Custom affine gap parameters (X=4, O=6, E=2)
cargo run --release -- -q query.fa -t target.fa -x 4 -o 6 -e 2
```

---

## Verification & Benchmarks

Run automated verification:
```bash
./scripts/test_example.sh
./scripts/benchmark_compare.sh
```

### Performance Across Divergence (2% Mutations):
| Sequence Length | WFA Latency | Gotoh DP Latency | Speedup | Score Match |
| :--- | :--- | :--- | :--- | :--- |
| **200 bp** | **0.08 ms** | 12.96 ms | **164.3×** | 100% Exact (16 == 16) |
| **500 bp** | **0.43 ms** | 103.72 ms | **242.0×** | 100% Exact (40 == 40) |
| **1,000 bp** | **3.02 ms** | 454.94 ms | **150.6×** | 100% Exact (80 == 80) |

---

## License

MIT License. Copyright (c) 2026 Andreas Sagen.
