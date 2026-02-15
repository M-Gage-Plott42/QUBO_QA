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
- Easy-case-rate metrics in `summary.json` and `scan_summary.csv`:
  - per-family easy-case rate at `k <= 0`
  - overall mean easy-case rate per run
- Success-probability reporting split into:
  - cumulative-by-time success (`success_prob.png`)
  - instantaneous-at-time success (`success_prob_instantaneous.png`)
- Ising energy evaluation aligned with `dimod` Ising/BINARY conversion convention (`x=(s+1)/2`).
- Optional `--n-list` scan-stop criterion based on:
  - max Holm-adjusted p-value threshold
  - overall mean easy-case-rate threshold
  - configurable easy-case comparison operator (`le`/`ge`) and minimum evaluated `n` count
- Graph relabeling for MaxCut/MIS (BFS high-degree ordering) to improve linear locality for MPS
- Modern Aer primitives sampling via `SamplerV2.from_backend(...)`
- Transpile/template reuse across instances with cache modes:
  - `support` mode (keyed by nonzero support signature)
  - `full` mode (keyed by `n` and schedule, broader reuse)
  - auto-disable fallback for low measured-template hit rates
- Optional Estimator diagnostics via `--estimator-diagnostics` and `--estimator-precision`
  - Writes `expectation_energy.png`
  - Adds estimator diagnostics keys to `summary.json`
- Run-level timing metadata in `summary.json` under `run_timing`
  - UTC start/end timestamps
  - total walltime and per-family walltime seconds

## Output Contract

Single-`n` mode outputs:
- `results.csv`
- `summary.json`
- `convergence_energy.png`
- `success_prob.png`
- `success_prob_instantaneous.png`
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
- `.venv/bin/python -m unittest discover -s tests -p 'test_*.py' -v`
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

## Cache Revisit Snapshot (Larger n, 2026-02-15)

- Bench root:
  - `diagnostics_local/2026-02-15/cache_revisit_s03/`
- Matrix:
  - `n=12, instances=10`, cache on/off
  - `n=14, instances=6`, cache on/off
  - common flags: `--aer-method statevector --opt-ref exact --no-plots`, `OMP_NUM_THREADS=11`
- Timing summary:
  - `n12_i10`: on `25.87s`, off `25.67s`
  - `n14_i6`: on `35.69s`, off `35.30s`
- Cache-stat summary:
  - cache-on runs showed zero measured-template hits in this workload (`measured_hits=0`).
- Outcome:
  - No observed speedup at these larger-`n` points for the tested instance mix and cache mode.

## Scan-Stop Threshold Tuning Snapshot (2026-02-15)

- Baseline scan (no scan-stop):
  - `diagnostics_local/2026-02-15/scan_stop_tune_baseline/scan_summary.csv`
- Aggressive preset audit:
  - `--scan-stop-p-holm-threshold 0.9`
  - `--scan-stop-easy-case-threshold 0.9`
  - `--scan-stop-easy-case-op ge`
  - `--scan-stop-min-n-evals 2`
  - Result: early-stop triggered at `n=5` in
    `diagnostics_local/2026-02-15/scan_stop_tune_aggressive/scan_summary.csv`
- Conservative preset audit:
  - `--scan-stop-p-holm-threshold 0.95`
  - `--scan-stop-easy-case-threshold 0.7`
  - `--scan-stop-easy-case-op ge`
  - `--scan-stop-min-n-evals 3`
  - Result: early-stop triggered at `n=6` in
    `diagnostics_local/2026-02-15/scan_stop_tune_conservative_trigger/scan_summary.csv`

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
- Scan-stop criterion has recommended presets, but defaults are still user-provided (no auto-selected policy).
- MPS scalability remains entanglement-dependent for dense couplings (expected limitation).

## Immediate Next Tasks

1. No required follow-up items are currently open from the method-audit action list.
2. Optional future work: revisit cache only if workload/instance structure changes enough to improve template-hit rates.

## Method Audit Action List (2026-02-15)

Priority `P0`:
- Enforce `--mis-lambda > 1` for current MIS objective form `-sum x_i + lambda * sum x_i x_j`.
  - Rationale: this preserves independent-set encoding in the standard QUBO construction.
  - Acceptance: CLI rejects invalid values with a clear error and README documents the bound.
  - Status: complete (implemented in commit `37eadaf`).

Priority `P1`:
- Clarify success metric labeling.
  - Current behavior computes success on `best_so_far` trajectories, so the plot is cumulative by time.
  - Acceptance: title/legend/docs explicitly say "reached by time t" and a separate instantaneous variant is added.
  - Status: complete (implemented on `main`, see commit history after `37eadaf`).

Priority `P1`:
- Add missing parameter guardrails.
  - Validate `--graph-p`, `--random-density` in `[0,1]`.
  - Acceptance: invalid values fail fast with actionable messages.
  - Status: complete (implemented on `main`, see commit history after `c6a5dac`).

Priority `P1`:
- Add focused tests for the above.
  - Acceptance: tests fail before fix and pass after fix for each targeted behavior.
  - Status: complete (implemented on `main`, see commit history after `b85f94c`).

Reference basis used for the action list:
- A. Lucas (2014), "Ising formulations of many NP problems" (MIS penalty condition, Eq. 37 / `B < A`): https://arxiv.org/pdf/1302.5843.pdf
- Qiskit bit-ordering guide (measurement-string conventions): https://qiskit.qotlabs.org/docs/guides/bit-ordering
- D-Wave `dimod` BQM docs (Ising/QUBO relationship and conventions): https://docs.dwavequantum.com/en/latest/ocean/api_ref_dimod/bqm.html

## Quick Commands

```bash
cd ~/git/QUBO_QA
git status
make smoke
make smoke-perm
make scan-smoke
```
