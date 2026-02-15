# QUBO QA Benchmark (Qiskit Aer, CPU)

Digital adiabatic/QA benchmark for three QUBO families:
- Random symmetric QUBO (`x^T Q x`)
- MaxCut QUBO on Erdos-Renyi graphs
- MIS QUBO on Erdos-Renyi graphs

The simulator uses Qiskit Aer with:
- `statevector` for small `n`
- `matrix_product_state` (MPS) for larger `n` when entanglement is moderate

## Current Features

- Trotterized adiabatic evolution with `delta_t` and sweep over step count `K`
- Convergence metric `steps_to_opt` (first `K` where best-so-far reaches reference energy)
- Three benchmark families: Random / MaxCut / MIS
- Success-probability curves vs time (`success_prob.png`)
- Modern Qiskit Aer primitives path for sampling (`SamplerV2.from_backend`)
- Optional modern diagnostics via `EstimatorV2` expectation-energy curves
- Graph relabeling (BFS from highest-degree node) for MaxCut/MIS to improve linear edge locality for MPS
- Parameterized transpile-template reuse across instances (can be disabled with `--no-transpile-cache`)
- Statistical testing modes:
  - `--stats-method mw` (one-sided Mann-Whitney U, random < comparator)
  - `--stats-method perm` (one-sided permutation test on median difference)
- Multi-`n` scan mode via `--n-list` with aggregated `scan_summary.csv`

## Quickstart

```bash
make install
make smoke
```

## Run Examples

Small exact sanity check:

```bash
.venv/bin/python qa_adiabatic_steps_bench.py -n 6 --instances 10 --t-max 5 --shots 64 --aer-method statevector --opt-ref exact
```

Permutation-statistics run:

```bash
.venv/bin/python qa_adiabatic_steps_bench.py -n 6 --instances 50 --t-max 10 --shots 128 \
  --aer-method statevector --opt-ref exact \
  --stats-method perm --perm-iterations 10000
```

`n`-scan run:

```bash
.venv/bin/python qa_adiabatic_steps_bench.py --n-list 4,5,6,7,8 \
  --instances 50 --t-max 10 --shots 128 \
  --aer-method statevector --opt-ref exact --outdir qa_scan
```

Push higher qubits with MPS:

```bash
.venv/bin/python qa_adiabatic_steps_bench.py -n 30 --instances 50 --t-max 10 --shots 64 \
  --aer-method matrix_product_state \
  --mps-max-bond-dimension 128 \
  --mps-truncation-threshold 1e-10
```

Try `n=42` (sparse random QUBO recommended):

```bash
.venv/bin/python qa_adiabatic_steps_bench.py -n 42 --instances 30 --t-max 8 --shots 32 \
  --aer-method matrix_product_state \
  --random-density 0.15 \
  --mps-max-bond-dimension 64 \
  --mps-truncation-threshold 1e-8
```

Enable estimator diagnostics (writes `expectation_energy.png`):

```bash
.venv/bin/python qa_adiabatic_steps_bench.py -n 6 --instances 10 --t-max 5 --shots 64 \
  --aer-method statevector --opt-ref exact --estimator-diagnostics
```

## Outputs

### Single-`n` mode

Default output directory: `qa_out/`

- `results.csv` (per-instance metrics)
- `summary.json` (aggregate metrics + test results)
- `convergence_energy.png` (median best-so-far energy vs time)
- `success_prob.png` (success probability vs time)
- `expectation_energy.png` (optional, only with `--estimator-diagnostics`)
- `steps_boxplot.png` (steps-to-opt distribution)

### `--n-list` scan mode

If `--outdir qa_scan` and `--n-list 4,5,6`, outputs are:
- `qa_scan/n_4/`, `qa_scan/n_5/`, `qa_scan/n_6/` (each with per-`n` results and summary)
- `qa_scan/scan_summary.csv` (one row per `n`, medians, p-values, Holm-adjusted p-values, Cliff's deltas)

## Output Policy

- Scratch/regression outputs remain local-only:
  - `qa_out*/`, `qa_scan*/`, and `out/` are git-ignored.
- Provenance outputs should be written under `artifacts/` and committed when needed.
- For reproducible provenance, always set an explicit run directory:
  - single-`n`: `artifacts/runs/<YYYY-MM-DD>/<run_tag>/`
  - `--n-list`: `artifacts/scans/<YYYY-MM-DD>/<run_tag>/`
- Recommended provenance files to commit:
  - `summary.json`
  - `scan_summary.csv` (for scans)
  - plots (`convergence_energy.png`, `success_prob.png`, `steps_boxplot.png`, optional `expectation_energy.png`)
  - `results.csv` when file size is reasonable for the repo.

## Method Notes

- This is digital adiabatic simulation (Trotterized circuit evolution), not analog hardware quantum annealing.
- `--opt-ref exact` uses brute-force exact solving (`dimod.ExactSolver`) and is intended for small `n`.
- `--opt-ref qa_best` uses the best observed final energy as the reference optimum.
- Statevector scaling is exponential in `n`.
- MPS scaling depends on entanglement growth and bond-dimension/truncation settings.
- For dense random couplings, MPS bond growth can still make large `n` expensive.
- Default transpile-template reuse reduces repeated compile overhead across instances; use `--no-transpile-cache` to disable.
- `--estimator-diagnostics` enables EstimatorV2 expectation-value tracking and adds estimator keys to `summary.json`.

## Qiskit v2 Alignment Notes (2026-02-15)

Validated against `M-Gage-Plott42/qiskit-v2-guide` (Qiskit 2.3.x patterns):
- Hamiltonian-to-gate mapping is consistent:
  - driver term `-sum X` is implemented as `rx(-2 t)`,
  - problem terms are implemented as `rz(2 t h)` and `rzz(2 t J)`.
- Measurement/energy path is consistent:
  - Qiskit counts are treated as MSB-first classical strings,
  - keys are reversed to qubit-0-first before Ising energy evaluation.
- Aer backend usage remains valid for Qiskit 2.x:
  - execution now uses `SamplerV2.from_backend(AerSimulator(...))` on transpiled circuits,
  - MPS options used by this script are recognized by Aer (`qiskit-aer 0.17.2`),
  - auto-switch to MPS at `n >= --mps-auto-threshold` is functioning.
- CLI and output contracts remain unchanged after this modernization.

## Makefile Shortcuts

- `make smoke` for quick baseline check
- `make smoke-perm` for permutation-test quick check
- `make scan-smoke` for quick `--n-list` check
