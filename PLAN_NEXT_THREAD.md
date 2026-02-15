# Plan: Next Thread Execution

## Current State

- Canonical repo for ongoing work: `~/git/QUBO_QA`
- Remote: `git@github.com:M-Gage-Plott42/QUBO_QA.git`
- Baseline benchmark script exists and runs:
  - `qa_adiabatic_steps_bench.py`
- Implemented features already in repo:
  - Success probability curve output (`success_prob.png`)
  - Graph relabeling for MaxCut/MIS (BFS high-degree ordering)
  - `--stats-method {mw,perm}` with Holm correction
  - `--n-list` scan mode + `scan_summary.csv`

## Primary Objective (Next)

Integrate and validate against your Qiskit guide repo to ensure modeling and API usage are aligned.

## Step Plan

1. Locate/access Qiskit guide repo.
- Needed input: local path or GitHub URL (and access details if private).

2. Perform alignment audit against this script.
- Validate Hamiltonian-to-gate mapping/sign conventions.
- Validate `rx/rz/rzz` angle factors.
- Validate Aer backend/MPS option usage and defaults.
- Validate measurement bit ordering and energy evaluation path.

3. Apply corrections if needed.
- Patch only minimal required code.
- Keep CLI and output contract stable unless change is required.

4. Re-validate behavior.
- `python -m py_compile qa_adiabatic_steps_bench.py`
- `make smoke`
- `make smoke-perm`
- `make scan-smoke`

5. Document outcomes.
- Update `README.md` with any changed assumptions/options.
- Add brief “Qiskit alignment notes” section summarizing what changed and why.

## Acceptance Criteria

- Script behavior remains reproducible for smoke runs.
- Any Qiskit correctness mismatches are fixed and documented.
- Repo is clean after commit(s), with clear commit messages.

## Command Reminders

```bash
cd ~/git/QUBO_QA
git status
make smoke
make smoke-perm
make scan-smoke
```

