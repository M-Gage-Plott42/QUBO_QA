# Project Status

## Last Updated

- Date: 2026-02-15
- Branch: `main`
- Head commit: `c9afda4`
- Repo: `git@github.com:M-Gage-Plott42/QUBO_QA.git`

## Purpose

Benchmark digital adiabatic/QA behavior for three QUBO families using Qiskit Aer on classical hardware:
- Random symmetric QUBO (`x^T Q x`)
- MaxCut QUBO (Erdos-Renyi)
- MIS QUBO (Erdos-Renyi)

Primary script: `qa_adiabatic_steps_bench.py`

## Implemented

- Trotterized adiabatic evolution with `rx`, `rz`, `rzz` and `H(s) = (1-s) H_driver + s H_problem`
- `steps_to_opt` convergence metric over `K` with `t = K * delta_t`
- Optional optimum references: `--opt-ref {qa_best,exact}`
- Statistical testing modes:
  - `--stats-method mw` (one-sided Mann-Whitney U)
  - `--stats-method perm` (one-sided permutation test on median difference)
  - Holm correction + Cliff's delta reporting
- `--n-list` multi-size scan mode and `scan_summary.csv`
- Graph relabeling for MaxCut/MIS (BFS high-degree ordering) to improve linear locality for MPS
- Modern Aer primitives sampling via `SamplerV2.from_backend(...)`
- Transpile/template reuse across instances with cache modes:
  - `support` mode (keyed by nonzero support signature)
  - `full` mode (keyed by `n` and schedule, broader reuse)
  - auto-disable fallback for low measured-template hit rates
- Optional Estimator diagnostics via `--estimator-diagnostics` and `--estimator-precision`
  - Writes `expectation_energy.png`
  - Adds estimator diagnostics keys to `summary.json`

## Output Contract

Single-`n` mode outputs:
- `results.csv`
- `summary.json`
- `convergence_energy.png`
- `success_prob.png`
- `steps_boxplot.png`
- `expectation_energy.png` (only when `--estimator-diagnostics` is enabled)

`--n-list` mode additionally outputs:
- `scan_summary.csv`
- per-`n` subdirectories (`n_<value>/...`)

## Repository Output Policy

- Scratch outputs are local-only and ignored by git:
  - `qa_out*/`, `qa_scan*/`, `out/`, `diagnostics_local/`, and `.mplconfig` cache directories.
- Provenance runs should use explicit `--outdir` paths under `artifacts/`.
- Preferred layout:
  - single-`n`: `artifacts/runs/<YYYY-MM-DD>/<run_tag>/`
  - `--n-list`: `artifacts/scans/<YYYY-MM-DD>/<run_tag>/`
- Diagnostic/exploratory runs should use:
  - `diagnostics_local/<YYYY-MM-DD>/<run_tag>/`
  - these remain untracked unless explicitly promoted.
- Commit provenance artifacts in separate commits from code changes.

## Validation (Most Recent Pass)

Executed successfully in this repo:
- `.venv/bin/python -m py_compile qa_adiabatic_steps_bench.py`
- `make smoke`
- `make smoke-perm`
- `make scan-smoke`
- Estimator diagnostics smoke:
  - `.venv/bin/python qa_adiabatic_steps_bench.py -n 4 --instances 2 --t-max 1 --shots 16 --aer-method statevector --opt-ref exact --estimator-diagnostics --outdir /tmp/qa_estimator_diag`
- Cache benchmark matrix restart (for cache on/off characterization):
  - `OMP_NUM_THREADS=11` run set under `diagnostics_local/2026-02-15/cache_bench_s02/` with per-run logs and `timing.csv`

## Cache Benchmark Snapshot (2026-02-15)

- Bench roots:
  - `diagnostics_local/2026-02-15/cache_bench_s01/` (partial run under 1-core contention)
  - `diagnostics_local/2026-02-15/cache_bench_s02/` (full restart with `OMP_NUM_THREADS=11`)
- `s02` matrix:
  - `n in {6,10}`, `instances in {10,25}`, cache `on/off`, `--aer-method statevector --opt-ref exact --no-plots`
- Timing summary (`diagnostics_local/2026-02-15/cache_bench_s02/timing.csv`):
  - `n6_i10`: on `9.42s`, off `8.00s`
  - `n6_i25`: on `18.90s`, off `19.35s`
  - `n10_i10`: on `16.26s`, off `14.60s`
  - `n10_i25`: on `37.39s`, off `36.97s`
- Outcome:
  - No consistent speedup from transpile cache in this tested range.
  - Cache benchmarking is deprioritized for now.

## Workstation Handoff Checklist

On a new workstation/clone, run:
- `make install`
- `make smoke`
- `make smoke-perm`
- `make scan-smoke`

Then confirm policy behavior:
- commit-intended outputs go under `artifacts/...`
- exploratory outputs go under `diagnostics_local/...` (ignored)

## Patch Compliance Policy (Required)

After any patch to benchmark code or CLI, run before commit:
- `.venv/bin/python -m py_compile qa_adiabatic_steps_bench.py`
- `make smoke`
- `make smoke-perm`
- `make scan-smoke`
- A targeted audit for the modified subsystem (for example cache on/off artifact audit, estimator diagnostics audit, or stats audit)

Commit only after these checks pass and audit results are reviewed.

## Open Risks / Gaps

- Cache performance remains workload-dependent; current `s02` matrix shows mixed `on` vs `off` timing with no clear speedup.
- No auto-stop/easy-case threshold logic in `--n-list` scans yet.
- MPS scalability remains entanglement-dependent for dense couplings (expected limitation).

## Immediate Next Tasks

1. Add an optional "easy-case rate" metric and include it in `summary.json` and `scan_summary.csv`.
2. Add an optional scan stopping criterion based on adjusted p-value and easy-case-rate threshold.
3. Optional/low priority: revisit cache benchmarks only if larger-`n` runs make transpile cost dominant.

## Quick Commands

```bash
cd ~/git/QUBO_QA
git status
make smoke
make smoke-perm
make scan-smoke
```
