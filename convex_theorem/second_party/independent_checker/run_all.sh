#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MALLOC_ARENA_MAX=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 PYTHONOPTIMIZE=0
PYTHON="${PYTHON:-python3}"
a=()
if [[ "${1:-}" == --smoke ]]; then shift; a+=(--sample 5 --work "$PWD/smoke_work" --jobs 2); fi
exec "$PYTHON" code/driver.py run "${a[@]}" "$@"
