# Project Status

## Last Updated

- Date: 2026-02-15
- Branch: `main`
- Head commit: `057f8da`
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
- Transpile/template reuse across instances (enabled by default; disable via `--no-transpile-cache`)
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

## Validation (Most Recent Pass)

Executed successfully in this repo:
- `.venv/bin/python -m py_compile qa_adiabatic_steps_bench.py`
- `make smoke`
- `make smoke-perm`
- `make scan-smoke`
- Estimator diagnostics smoke:
  - `.venv/bin/python qa_adiabatic_steps_bench.py -n 4 --instances 2 --t-max 1 --shots 16 --aer-method statevector --opt-ref exact --estimator-diagnostics --outdir /tmp/qa_estimator_diag`

## Open Risks / Gaps

- Runtime/perf characterization for transpile-template cache (`on` vs `off`) is not yet benchmarked and documented.
- No auto-stop/easy-case threshold logic in `--n-list` scans yet.
- MPS scalability remains entanglement-dependent for dense couplings (expected limitation).

## Immediate Next Tasks

1. Run a controlled performance benchmark comparing cache enabled/disabled across representative `n` and instance counts.
2. Add an optional "easy-case rate" metric and include it in `summary.json` and `scan_summary.csv`.
3. Add an optional scan stopping criterion based on adjusted p-value and easy-case-rate threshold.

## Quick Commands

```bash
cd ~/git/QUBO_QA
git status
make smoke
make smoke-perm
make scan-smoke
```
