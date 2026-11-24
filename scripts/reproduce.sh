#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/erdos99-mpl
python -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python -m erdos99.experiment --n-min 4 --n-max 30 --seeds 99 2026 31415 --margins 0.001 0.01 --maxiter 600 --output "${1:-results/reproduction}"
.venv/bin/python -m erdos99.analyze "${1:-results/reproduction}"
