# AGENTS.md

## Purpose

This repo benchmarks digital adiabatic/QA behavior for three QUBO families using Qiskit Aer:

- Random symmetric QUBO
- MaxCut QUBO
- MIS QUBO

Primary script: `qa_adiabatic_steps_bench.py`

## Environment

- Preferred development: WSL/Linux filesystem clone (for performance and fewer path/permission issues).
- Current mapped path may be `/mnt/c/git/QUBO_QA` (Windows `C:\git\QUBO_QA`).
- Python dependencies are in `requirements.txt`.

## Core Commands

- Install deps: `make install`
- Baseline smoke: `make smoke`
- Permutation-stat smoke: `make smoke-perm`
- n-scan smoke: `make scan-smoke`

Direct examples:

- Single n: `python qa_adiabatic_steps_bench.py -n 6 --instances 10 --t-max 5 --shots 64 --aer-method statevector --opt-ref exact`
- Scan mode: `python qa_adiabatic_steps_bench.py --n-list 4,5,6 --instances 50 --t-max 10 --shots 128 --opt-ref exact`

## Implementation Notes

- `num_processes=1` is set in `transpile` to avoid multiprocessing semaphore issues in constrained/sandboxed environments.
- Matplotlib config directory is set under output dirs to avoid cache permission issues.
- Graph-based families use BFS-high-degree relabeling to improve linear locality for MPS.
- Statistical modes:
  - `--stats-method mw` for one-sided Mann-Whitney U.
  - `--stats-method perm` for one-sided permutation test on median differences.

## Output Contract

Single-`n` run should produce:

- `results.csv`
- `summary.json`
- `convergence_energy.png`
- `success_prob.png`
- `steps_boxplot.png`

`--n-list` scan should additionally produce:

- `scan_summary.csv`
- Per-`n` subdirectories (`n_<value>/...`)

## Change Workflow

- Keep CLI flags backward-compatible unless explicitly changing interface.
- Prefer small, isolated edits with clear validation.
- Update `README.md` when adding/removing CLI options or output files.

## Validation Checklist

Before committing:

- `python -m py_compile qa_adiabatic_steps_bench.py`
- `make smoke`
- If stats logic changed: `make smoke-perm`
- If scan logic changed: `make scan-smoke`
- Confirm expected artifacts and keys in `summary.json`.

## Safety

- Do not run destructive git resets on user work.
- Do not remove user data/output directories unless explicitly requested.
- Treat WSL path differences as a first-class risk when reproducing commands.
