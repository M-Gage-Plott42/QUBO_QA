#!/usr/bin/env python3
"""
Prototype PRR-style analog local-detuning optimization.

This script is intentionally separate from qa_adiabatic_steps_bench.py.
It performs piecewise-constant pulse optimization using a BFGS/NM/BFGS
sequence for an analog Ising evolution model with:
  H(t) = Omega(t) * H_driver + g(t) * H_linear + H_quadratic

Use this for exploratory PRR-style studies; keep the digital K-sweep
benchmark in qa_adiabatic_steps_bench.py unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
from typing import Any, Dict, List, Tuple

import dimod
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from qa_adiabatic_steps_bench import (
    approximation_ratio_minimization,
    bqm_to_ising_arrays,
    ising_energy_from_bitstring,
    make_random_qubo_bqm,
    make_maxcut_qubo_bqm,
    make_mis_qubo_bqm,
)


def kron_all(mats: List[np.ndarray]) -> np.ndarray:
    out = np.array([[1.0 + 0.0j]], dtype=complex)
    for mat in mats:
        out = np.kron(out, mat)
    return out


def build_pauli_ops(n: int) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    i2 = np.eye(2, dtype=complex)
    x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

    x_ops: List[np.ndarray] = []
    z_ops: List[np.ndarray] = []
    for q in range(n):
        mats_x = [i2] * n
        mats_z = [i2] * n
        mats_x[q] = x
        mats_z[q] = z
        x_ops.append(kron_all(mats_x))
        z_ops.append(kron_all(mats_z))
    return x_ops, z_ops


def build_hamiltonians(
    n: int,
    h: np.ndarray,
    j_terms: List[Tuple[int, int, float]],
    offset: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_ops, z_ops = build_pauli_ops(n)
    dim = 2**n
    identity = np.eye(dim, dtype=complex)

    h_driver = np.zeros((dim, dim), dtype=complex)
    for op in x_ops:
        h_driver += -op

    h_linear = np.zeros((dim, dim), dtype=complex)
    for i in range(n):
        if h[i] != 0.0:
            h_linear += float(h[i]) * z_ops[i]

    h_quadratic = float(offset) * identity
    for (i, j, jij) in j_terms:
        if jij != 0.0:
            h_quadratic += float(jij) * (z_ops[i] @ z_ops[j])

    return h_driver, h_linear, h_quadratic


def evolve_piecewise_constant(
    omega_vals: np.ndarray,
    g_vals: np.ndarray,
    total_time: float,
    h_driver: np.ndarray,
    h_linear: np.ndarray,
    h_quadratic: np.ndarray,
    psi0: np.ndarray,
) -> np.ndarray:
    segments = int(len(omega_vals))
    dt = float(total_time) / float(segments)
    psi = psi0.copy()
    for k in range(segments):
        h_seg = float(omega_vals[k]) * h_driver + float(g_vals[k]) * h_linear + h_quadratic
        psi = expm(-1.0j * dt * h_seg) @ psi
    return psi


def basis_index_to_bitstring_lsb0(index: int, n: int) -> str:
    bits_msb = format(int(index), f"0{n}b")
    return bits_msb[::-1]


def optimize_local_detuning_schedule(
    h_driver: np.ndarray,
    h_linear: np.ndarray,
    h_quadratic: np.ndarray,
    target_hamiltonian: np.ndarray,
    n: int,
    total_time: float,
    segments: int,
    omega_max: float,
    g_min: float,
    g_max: float,
    g_start: float,
    maxiter_bfgs: int,
    maxiter_nm: int,
) -> Dict[str, Any]:
    dim = 2**n
    psi0 = np.zeros(dim, dtype=complex)
    psi0[0] = 1.0 + 0.0j

    if segments < 3:
        raise ValueError("segments must be at least 3.")

    # Optimized variables:
    # - omega for interior segments (endpoints fixed to 0)
    # - g for all but last segment (last fixed to 1.0)
    n_omega = segments - 2
    n_g = segments - 1

    def unpack(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        omega = np.zeros(segments, dtype=float)
        gvals = np.zeros(segments, dtype=float)
        omega[1:-1] = np.clip(x[:n_omega], 0.0, omega_max)
        gvals[:-1] = np.clip(x[n_omega:], g_min, g_max)
        gvals[-1] = 1.0
        return omega, gvals

    omega_guess = omega_max * np.sin(np.linspace(0.0, np.pi, segments))[1:-1]
    g_guess_full = np.linspace(float(g_start), 1.0, segments, dtype=float)
    x0 = np.concatenate([omega_guess, g_guess_full[:-1]])

    bounds = [(0.0, float(omega_max))] * n_omega + [(float(g_min), float(g_max))] * n_g
    eval_history: List[float] = []

    def objective(x: np.ndarray) -> float:
        omega, gvals = unpack(x)
        psi = evolve_piecewise_constant(
            omega_vals=omega,
            g_vals=gvals,
            total_time=total_time,
            h_driver=h_driver,
            h_linear=h_linear,
            h_quadratic=h_quadratic,
            psi0=psi0,
        )
        e = float(np.real(np.vdot(psi, target_hamiltonian @ psi)))
        eval_history.append(e)
        return e

    t0 = time.perf_counter()
    eval_start_stage1 = int(len(eval_history))
    stage1 = minimize(
        objective,
        x0=x0,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": int(maxiter_bfgs)},
    )
    eval_end_stage1 = int(len(eval_history))
    t1 = time.perf_counter()
    eval_start_stage2 = int(len(eval_history))
    stage2 = minimize(
        objective,
        x0=stage1.x,
        method="Nelder-Mead",
        options={"maxiter": int(maxiter_nm), "xatol": 1e-4, "fatol": 1e-4},
    )
    eval_end_stage2 = int(len(eval_history))
    t2 = time.perf_counter()
    eval_start_stage3 = int(len(eval_history))
    stage3 = minimize(
        objective,
        x0=stage2.x,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": int(maxiter_bfgs)},
    )
    eval_end_stage3 = int(len(eval_history))
    t3 = time.perf_counter()

    omega_opt, g_opt = unpack(stage3.x)
    psi_final = evolve_piecewise_constant(
        omega_vals=omega_opt,
        g_vals=g_opt,
        total_time=total_time,
        h_driver=h_driver,
        h_linear=h_linear,
        h_quadratic=h_quadratic,
        psi0=psi0,
    )
    e_final = float(np.real(np.vdot(psi_final, target_hamiltonian @ psi_final)))

    return {
        "omega_opt": omega_opt,
        "g_opt": g_opt,
        "psi_final": psi_final,
        "objective_final": e_final,
        "eval_history": eval_history,
        "stage_eval_ranges": {
            "stage1": {"start": eval_start_stage1, "end": eval_end_stage1},
            "stage2": {"start": eval_start_stage2, "end": eval_end_stage2},
            "stage3": {"start": eval_start_stage3, "end": eval_end_stage3},
        },
        "stage": {
            "stage1": {
                "success": bool(stage1.success),
                "nit": int(stage1.nit) if stage1.nit is not None else None,
                "fun": float(stage1.fun),
                "message": str(stage1.message),
                "walltime_seconds": float(t1 - t0),
            },
            "stage2": {
                "success": bool(stage2.success),
                "nit": int(stage2.nit) if stage2.nit is not None else None,
                "fun": float(stage2.fun),
                "message": str(stage2.message),
                "walltime_seconds": float(t2 - t1),
            },
            "stage3": {
                "success": bool(stage3.success),
                "nit": int(stage3.nit) if stage3.nit is not None else None,
                "fun": float(stage3.fun),
                "message": str(stage3.message),
                "walltime_seconds": float(t3 - t2),
            },
            "total_walltime_seconds": float(t3 - t0),
        },
    }


def protocol_mode_for_problem(problem: str) -> str:
    if str(problem) in {"maxcut", "mis"}:
        return "paper_aligned"
    return "exploratory_generalized"


def _stage_label_curve(num_points: int, stage_eval_ranges: Dict[str, Dict[str, int]]) -> List[str]:
    labels = ["unknown"] * int(num_points)
    for stage_name in ("stage1", "stage2", "stage3"):
        section = stage_eval_ranges.get(stage_name, {})
        start = int(section.get("start", 0))
        end = int(section.get("end", 0))
        for idx in range(max(0, start), min(int(num_points), end)):
            labels[idx] = stage_name
    return labels


def build_prr_objective_trace(
    eval_history: List[float],
    stage_eval_ranges: Dict[str, Dict[str, int]],
) -> Dict[str, List[Any]]:
    n = int(len(eval_history))
    if n <= 0:
        return {
            "eval_index": [],
            "budget_fraction": [],
            "objective_energy": [],
            "best_objective_so_far": [],
            "stage": [],
        }

    objective = np.asarray(eval_history, dtype=float)
    best_so_far = np.minimum.accumulate(objective)
    budget = (np.arange(n, dtype=float) + 1.0) / float(n)
    stage_labels = _stage_label_curve(num_points=n, stage_eval_ranges=stage_eval_ranges)
    return {
        "eval_index": [int(v) for v in np.arange(n, dtype=int)],
        "budget_fraction": [float(v) for v in budget],
        "objective_energy": [float(v) for v in objective],
        "best_objective_so_far": [float(v) for v in best_so_far],
        "stage": stage_labels,
    }


def make_problem_bqm(
    *,
    problem: str,
    rng: np.random.Generator,
    n: int,
    graph_p: float,
    random_low: int,
    random_high: int,
    random_density: float,
    maxcut_weight_low: int,
    maxcut_weight_high: int,
    mis_node_weight_low: int,
    mis_node_weight_high: int,
    mis_lambda: float,
) -> Tuple[dimod.BinaryQuadraticModel, Dict[str, Any]]:
    if str(problem) == "maxcut":
        return make_maxcut_qubo_bqm(
            rng=rng,
            n=n,
            p=float(graph_p),
            weight_low=int(maxcut_weight_low),
            weight_high=int(maxcut_weight_high),
        )
    if str(problem) == "mis":
        return make_mis_qubo_bqm(
            rng=rng,
            n=n,
            p=float(graph_p),
            penalty_lambda=float(mis_lambda),
            node_weight_low=int(mis_node_weight_low),
            node_weight_high=int(mis_node_weight_high),
        )
    if str(problem) == "random":
        return make_random_qubo_bqm(
            rng=rng,
            n=n,
            low=int(random_low),
            high=int(random_high),
            density=float(random_density),
        )
    raise ValueError(f"Unsupported problem: {problem}")


def run_prr_local_detuning_for_bqm(
    *,
    bqm: dimod.BinaryQuadraticModel,
    problem: str,
    meta: Dict[str, Any],
    n: int,
    graph_p: float,
    seed: int,
    total_time: float,
    segments: int,
    omega_max: float,
    g_min: float,
    g_max: float,
    g_start: float,
    maxiter_bfgs: int,
    maxiter_nm: int,
    exact_max_n: int,
) -> Dict[str, Any]:
    model = bqm_to_ising_arrays(bqm, n=n)
    h_driver, h_linear, h_quadratic = build_hamiltonians(
        n=n,
        h=model.h,
        j_terms=model.j_terms,
        offset=model.offset,
    )
    target_hamiltonian = h_linear + h_quadratic

    opt_result = optimize_local_detuning_schedule(
        h_driver=h_driver,
        h_linear=h_linear,
        h_quadratic=h_quadratic,
        target_hamiltonian=target_hamiltonian,
        n=n,
        total_time=float(total_time),
        segments=int(segments),
        omega_max=float(omega_max),
        g_min=float(g_min),
        g_max=float(g_max),
        g_start=float(g_start),
        maxiter_bfgs=int(maxiter_bfgs),
        maxiter_nm=int(maxiter_nm),
    )

    psi_final = opt_result["psi_final"]
    probs = np.abs(psi_final) ** 2
    top_idx = int(np.argmax(probs))
    top_bitstring_lsb0 = basis_index_to_bitstring_lsb0(top_idx, n=n)
    sampled_best_energy = ising_energy_from_bitstring(top_bitstring_lsb0, model)
    expectation_energy_final = float(opt_result["objective_final"])

    exact_opt_energy: float | None = None
    approx_ratio: float | None = None
    opt_subspace_probability: float | None = None
    if n <= int(exact_max_n):
        exact = dimod.ExactSolver().sample(bqm)
        exact_opt_energy = float(exact.first.energy)
        approx_ratio = approximation_ratio_minimization(
            obtained=sampled_best_energy,
            optimal=exact_opt_energy,
        )
        opt_prob = 0.0
        lowest = exact.lowest(atol=1e-9, rtol=0.0)
        for rec in lowest.data(fields=["sample"]):
            sample = rec.sample
            bits_msb = "".join("1" if int(sample[i]) == 1 else "0" for i in range(n))
            opt_prob += float(probs[int(bits_msb, 2)])
        opt_subspace_probability = float(opt_prob)

    objective_trace = build_prr_objective_trace(
        eval_history=opt_result["eval_history"],
        stage_eval_ranges=opt_result["stage_eval_ranges"],
    )
    mode = protocol_mode_for_problem(problem)
    summary: Dict[str, Any] = {
        "problem": str(problem),
        "protocol_mode": mode,
        "mode_note": (
            "paper-aligned family (PRR-style MaxCut/MIS local-detuning protocol)."
            if mode == "paper_aligned"
            else "exploratory extension beyond PRR paper scope."
        ),
        "n": int(n),
        "graph_p": float(graph_p),
        "seed": int(seed),
        "segments": int(segments),
        "total_time": float(total_time),
        "omega_max": float(omega_max),
        "g_min": float(g_min),
        "g_max": float(g_max),
        "g_start": float(g_start),
        "meta": meta,
        "sampled_best_bitstring_lsb0": top_bitstring_lsb0,
        "sampled_best_energy": float(sampled_best_energy),
        "expectation_energy_final": float(expectation_energy_final),
        "exact_opt_energy": exact_opt_energy,
        "approximation_ratio_vs_exact": approx_ratio,
        "opt_subspace_probability": opt_subspace_probability,
        "schedule": {
            "omega": [float(v) for v in opt_result["omega_opt"]],
            "detuning_scale_g": [float(v) for v in opt_result["g_opt"]],
        },
        "optimization": {
            "stage": opt_result["stage"],
            "stage_eval_ranges": opt_result["stage_eval_ranges"],
            "num_objective_evals": int(len(opt_result["eval_history"])),
            "best_objective_seen": float(min(opt_result["eval_history"])) if opt_result["eval_history"] else None,
            "trace_budget_fraction": objective_trace["budget_fraction"],
            "trace_objective_energy": objective_trace["objective_energy"],
            "trace_best_objective_so_far": objective_trace["best_objective_so_far"],
            "trace_stage": objective_trace["stage"],
        },
    }
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", type=str, required=True, choices=["random", "maxcut", "mis"])
    ap.add_argument("-n", "--n", type=int, required=True)
    ap.add_argument("--graph-p", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)

    ap.add_argument("--random-low", type=int, default=-5)
    ap.add_argument("--random-high", type=int, default=5)
    ap.add_argument("--random-density", type=float, default=1.0)
    ap.add_argument("--maxcut-weight-low", type=int, default=1)
    ap.add_argument("--maxcut-weight-high", type=int, default=5)
    ap.add_argument("--mis-node-weight-low", type=int, default=1)
    ap.add_argument("--mis-node-weight-high", type=int, default=4)
    ap.add_argument("--mis-lambda", type=float, default=5.0)

    ap.add_argument("--total-time", type=float, default=3.5)
    ap.add_argument("--segments", type=int, default=8)
    ap.add_argument("--omega-max", type=float, default=5.0)
    ap.add_argument("--g-min", type=float, default=-8.0)
    ap.add_argument("--g-max", type=float, default=4.0)
    ap.add_argument("--g-start", type=float, default=-3.0)
    ap.add_argument("--maxiter-bfgs", type=int, default=120)
    ap.add_argument("--maxiter-nm", type=int, default=180)

    ap.add_argument("--exact-max-n", type=int, default=18)
    ap.add_argument("--outdir", type=str, default="diagnostics_local/prr_local_detuning_run")
    args = ap.parse_args()

    if int(args.n) <= 0:
        raise SystemExit("--n must be positive.")
    if not (0.0 <= float(args.graph_p) <= 1.0):
        raise SystemExit("--graph-p must be in [0, 1].")
    if int(args.segments) < 3:
        raise SystemExit("--segments must be at least 3.")
    if float(args.total_time) <= 0.0:
        raise SystemExit("--total-time must be positive.")
    if float(args.omega_max) <= 0.0:
        raise SystemExit("--omega-max must be positive.")
    if float(args.g_min) >= float(args.g_max):
        raise SystemExit("--g-min must be < --g-max.")
    if int(args.maxiter_bfgs) <= 0 or int(args.maxiter_nm) <= 0:
        raise SystemExit("--maxiter-bfgs and --maxiter-nm must be positive.")
    if int(args.random_low) > int(args.random_high):
        raise SystemExit("--random-low must be <= --random-high.")
    if not (0.0 <= float(args.random_density) <= 1.0):
        raise SystemExit("--random-density must be in [0, 1].")
    if int(args.maxcut_weight_low) > int(args.maxcut_weight_high):
        raise SystemExit("--maxcut-weight-low must be <= --maxcut-weight-high.")
    if int(args.mis_node_weight_low) > int(args.mis_node_weight_high):
        raise SystemExit("--mis-node-weight-low must be <= --mis-node-weight-high.")
    if float(args.mis_lambda) <= max(1.0, float(args.mis_node_weight_high)):
        raise SystemExit("--mis-lambda must be > max(1.0, --mis-node-weight-high).")

    rng = np.random.default_rng(int(args.seed))
    n = int(args.n)

    bqm, meta = make_problem_bqm(
        problem=str(args.problem),
        rng=rng,
        n=n,
        graph_p=float(args.graph_p),
        random_low=int(args.random_low),
        random_high=int(args.random_high),
        random_density=float(args.random_density),
        maxcut_weight_low=int(args.maxcut_weight_low),
        maxcut_weight_high=int(args.maxcut_weight_high),
        mis_node_weight_low=int(args.mis_node_weight_low),
        mis_node_weight_high=int(args.mis_node_weight_high),
        mis_lambda=float(args.mis_lambda),
    )

    summary = run_prr_local_detuning_for_bqm(
        bqm=bqm,
        problem=str(args.problem),
        meta=meta,
        n=n,
        graph_p=float(args.graph_p),
        seed=int(args.seed),
        total_time=float(args.total_time),
        segments=int(args.segments),
        omega_max=float(args.omega_max),
        g_min=float(args.g_min),
        g_max=float(args.g_max),
        g_start=float(args.g_start),
        maxiter_bfgs=int(args.maxiter_bfgs),
        maxiter_nm=int(args.maxiter_nm),
        exact_max_n=int(args.exact_max_n),
    )

    os.makedirs(args.outdir, exist_ok=True)

    summary_path = os.path.join(args.outdir, "summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote: {summary_path}")

    trace_path = os.path.join(args.outdir, "optimization_trace.csv")
    with open(trace_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["eval_index", "budget_fraction", "objective_energy", "best_objective_so_far", "stage"])
        opt = summary["optimization"]
        for idx in range(int(opt["num_objective_evals"])):
            w.writerow(
                [
                    int(idx),
                    float(opt["trace_budget_fraction"][idx]),
                    float(opt["trace_objective_energy"][idx]),
                    float(opt["trace_best_objective_so_far"][idx]),
                    str(opt["trace_stage"][idx]),
                ]
            )
    print(f"Wrote: {trace_path}")


if __name__ == "__main__":
    main()
