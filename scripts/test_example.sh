#!/usr/bin/env bash
set -euo pipefail

echo "=== 1. Verifying WFA Algorithmic Equivalence ==="
python3 scripts/verify_logic.py

if command -v cargo >/dev/null 2>&1; then
    echo -e "\n=== 2. Running Rust CLI ==="
    cargo run --release -- -q data/sample_query.fa -t data/sample_target.fa -m wfa -s
else
    echo -e "\n[Note: cargo toolchain not yet installed; algorithmic verification ran via Python suite]"
fi

echo -e "\nExample run completed successfully!"
