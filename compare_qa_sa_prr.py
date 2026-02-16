#!/usr/bin/env python3
"""
Run side-by-side QA, SA, and PRR comparisons on shared QUBO instances.

Outputs:
- instance_bank.json
- comparison_results.csv
- comparison_curves.csv
- comparison_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
from typing import Any, Dict, List, Optional

import dimod
import numpy as np

from comparison_common import (
    COMPARISON_SCHEMA_VERSION,
    FAMILY_ORDER,
    build_instance_bank,
    write_instance_bank_manifest,
)
from prr_local_detuning_opt import run_prr_local_detuning_for_bqm
from qa_adiabatic_steps_bench import (
    approximation_ratio_minimization,
    bqm_to_ising_arrays,
    build_backend,
    run_instance_sweep,
    run_sa_convergence_trace,
    scale_ising_for_dynamics,
)


def parse_csv_items(csv_text: str) -> List[str]:
    return [token.strip() for token in str(csv_text).split(",") if token.strip()]


def parse_sa_sweep_checkpoints(csv_text: str) -> List[int]:
    values: List[int] = []
    for item in parse_csv_items(csv_text):
        ivalue = int(item)
        if ivalue <= 0:
            raise ValueError("All SA sweep checkpoints must be positive integers.")
        values.append(ivalue)
    if not values:
        raise ValueError("At least one SA sweep checkpoint is required.")
    return sorted(set(values))


def resolve_family_order(csv_text: str) -> List[str]:
    requested = parse_csv_items(csv_text)
    if not requested:
        raise ValueError("At least one family must be selected.")
    valid = set(FAMILY_ORDER)
    if any(fam not in valid for fam in requested):
        raise ValueError(f"Unsupported family in --families. Allowed: {','.join(FAMILY_ORDER)}")
    return requested


def _metric_or_none(values: List[Optional[float]], fn: Any) -> Optional[float]:
    filtered = [float(v) for v in values if v is not None]
    if not filtered:
        return None
    return float(fn(np.asarray(filtered, dtype=float)))


def _bool_rate(values: List[Optional[bool]]) -> Optional[float]:
    filtered = [bool(v) for v in values if v is not None]
    if not filtered:
        return None
    return float(np.mean(np.asarray(filtered, dtype=float)))


def build_summary_aggregates(result_rows: List[Dict[str, Any]], families: List[str]) -> Dict[str, Any]:
    summary: Dict[str, Any] = {}
    algorithms = ["qa", "sa", "prr"]

    for family in families:
        summary[family] = {}
        for algo in algorithms:
            rows = [row for row in result_rows if row["family"] == family and row["algorithm"] == algo]
            approx_values = [row["approx_ratio"] for row in rows]
            ratio_error_values = [row["approx_ratio_error"] for row in rows]
            gap_values = [row["gap_to_opt"] for row in rows]
            runtime_values = [row["runtime_seconds"] for row in rows]
            success_values = [row["success_vs_opt"] for row in rows]

            summary[family][algo] = {
                "instances": int(len(rows)),
                "mean_final_energy": _metric_or_none([row["final_energy"] for row in rows], np.mean),
                "median_final_energy": _metric_or_none([row["final_energy"] for row in rows], np.median),
                "approx_ratio_mean": _metric_or_none(approx_values, np.mean),
                "approx_ratio_median": _metric_or_none(approx_values, np.median),
                "approx_ratio_error_mean_1_minus_r": _metric_or_none(ratio_error_values, np.mean),
                "gap_to_opt_mean": _metric_or_none(gap_values, np.mean),
                "gap_to_opt_median": _metric_or_none(gap_values, np.median),
                "runtime_seconds_mean": _metric_or_none(runtime_values, np.mean),
                "runtime_seconds_median": _metric_or_none(runtime_values, np.median),
                "success_rate_vs_opt": _bool_rate(success_values),
            }
    return summary


def _rows_to_curve_index(curve_rows: List[Dict[str, Any]]) -> Dict[tuple, Dict[str, np.ndarray]]:
    grouped: Dict[tuple, List[Dict[str, Any]]] = {}
    for row in curve_rows:
        key = (str(row["family"]), int(row["instance_id"]), str(row["algorithm"]))
        grouped.setdefault(key, []).append(row)

    out: Dict[tuple, Dict[str, np.ndarray]] = {}
    for key, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda r: float(r["budget_fraction"]))
        budget = np.asarray([float(r["budget_fraction"]) for r in rows_sorted], dtype=float)
        energy = np.asarray([float(r["best_energy_so_far"]) for r in rows_sorted], dtype=float)
        out[key] = {"budget": budget, "energy": energy}
    return out


def _write_plots(
    *,
    outdir: str,
    families: List[str],
    result_rows: List[Dict[str, Any]],
    curve_rows: List[Dict[str, Any]],
    plot_grid_points: int,
) -> List[str]:
    mpl_config_dir = os.path.join(outdir, ".mplconfig")
    os.makedirs(mpl_config_dir, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = mpl_config_dir

    import matplotlib.pyplot as plt

    saved: List[str] = []
    algorithms = ["qa", "sa", "prr"]
    algo_labels = {"qa": "QA", "sa": "SA", "prr": "PRR"}
    algo_colors = {"qa": "#1f77b4", "sa": "#2ca02c", "prr": "#d62728"}
    tol = 1e-9

    opt_by_instance: Dict[tuple, Optional[float]] = {}
    for row in result_rows:
        key = (str(row["family"]), int(row["instance_id"]))
        if key not in opt_by_instance:
            opt_by_instance[key] = row["opt_energy"]

    curve_index = _rows_to_curve_index(curve_rows)
    grid = np.linspace(0.0, 1.0, int(max(8, plot_grid_points)))

    # Convergence ratio R by family.
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.2), squeeze=False)
    axes_list = axes[0]
    ratio_plot_has_data = False
    for ax, family in zip(axes_list, families):
        ax.set_title(f"{family} (R)")
        for algo in algorithms:
            curves: List[np.ndarray] = []
            for inst_id in sorted({int(r["instance_id"]) for r in result_rows if r["family"] == family}):
                curve_key = (family, inst_id, algo)
                inst_key = (family, inst_id)
                if curve_key not in curve_index:
                    continue
                opt = opt_by_instance.get(inst_key)
                if opt is None:
                    continue
                budget = curve_index[curve_key]["budget"]
                energy = curve_index[curve_key]["energy"]
                interp = np.interp(grid, budget, energy)
                ratio_curve = np.asarray(
                    [approximation_ratio_minimization(obtained=float(e), optimal=float(opt)) for e in interp],
                    dtype=float,
                )
                curves.append(ratio_curve)
            if curves:
                ratio_plot_has_data = True
                stacked = np.vstack(curves)
                ax.plot(grid, np.median(stacked, axis=0), label=algo_labels[algo], color=algo_colors[algo])
        ax.set_xlabel("budget fraction")
        ax.set_ylabel("median approximation ratio R")
        ax.set_ylim(0.0, 1.01)
        ax.grid(alpha=0.25)
        ax.legend()
    fig.tight_layout()
    ratio_path = os.path.join(outdir, "convergence_ratio_compare.png")
    fig.savefig(ratio_path, dpi=200)
    plt.close(fig)
    if ratio_plot_has_data:
        saved.append(ratio_path)

    # Success probability by family.
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.2), squeeze=False)
    axes_list = axes[0]
    success_plot_has_data = False
    for ax, family in zip(axes_list, families):
        ax.set_title(f"{family} success")
        for algo in algorithms:
            curves: List[np.ndarray] = []
            for inst_id in sorted({int(r["instance_id"]) for r in result_rows if r["family"] == family}):
                curve_key = (family, inst_id, algo)
                inst_key = (family, inst_id)
                if curve_key not in curve_index:
                    continue
                opt = opt_by_instance.get(inst_key)
                if opt is None:
                    continue
                budget = curve_index[curve_key]["budget"]
                energy = curve_index[curve_key]["energy"]
                interp = np.interp(grid, budget, energy)
                curves.append((interp <= float(opt) + tol).astype(float))
            if curves:
                success_plot_has_data = True
                stacked = np.vstack(curves)
                ax.plot(grid, np.mean(stacked, axis=0), label=algo_labels[algo], color=algo_colors[algo])
        ax.set_xlabel("budget fraction")
        ax.set_ylabel("success probability")
        ax.set_ylim(0.0, 1.01)
        ax.grid(alpha=0.25)
        ax.legend()
    fig.tight_layout()
    success_path = os.path.join(outdir, "success_prob_compare.png")
    fig.savefig(success_path, dpi=200)
    plt.close(fig)
    if success_plot_has_data:
        saved.append(success_path)

    # Histogram: final approximation ratio by family.
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.2), squeeze=False)
    axes_list = axes[0]
    ratio_hist_has_data = False
    for ax, family in zip(axes_list, families):
        ax.set_title(f"{family} final R")
        for algo in algorithms:
            vals = [
                float(r["approx_ratio"])
                for r in result_rows
                if r["family"] == family and r["algorithm"] == algo and r["approx_ratio"] is not None
            ]
            if vals:
                ratio_hist_has_data = True
                ax.hist(vals, bins=10, alpha=0.5, label=algo_labels[algo], color=algo_colors[algo])
        ax.set_xlabel("approximation ratio R")
        ax.set_ylabel("count")
        ax.grid(alpha=0.2)
        ax.legend()
    fig.tight_layout()
    ratio_hist_path = os.path.join(outdir, "hist_final_ratio.png")
    fig.savefig(ratio_hist_path, dpi=200)
    plt.close(fig)
    if ratio_hist_has_data:
        saved.append(ratio_hist_path)

    # Histogram: final gap to opt by family.
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.2), squeeze=False)
    axes_list = axes[0]
    gap_hist_has_data = False
    for ax, family in zip(axes_list, families):
        ax.set_title(f"{family} final gap")
        for algo in algorithms:
            vals = [
                float(r["gap_to_opt"])
                for r in result_rows
                if r["family"] == family and r["algorithm"] == algo and r["gap_to_opt"] is not None
            ]
            if vals:
                gap_hist_has_data = True
                ax.hist(vals, bins=10, alpha=0.5, label=algo_labels[algo], color=algo_colors[algo])
        ax.set_xlabel("gap to optimum")
        ax.set_ylabel("count")
        ax.grid(alpha=0.2)
        ax.legend()
    fig.tight_layout()
    gap_hist_path = os.path.join(outdir, "hist_final_gap.png")
    fig.savefig(gap_hist_path, dpi=200)
    plt.close(fig)
    if gap_hist_has_data:
        saved.append(gap_hist_path)

    # Histogram: runtime by family.
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.2), squeeze=False)
    axes_list = axes[0]
    runtime_hist_has_data = False
    for ax, family in zip(axes_list, families):
        ax.set_title(f"{family} runtime")
        for algo in algorithms:
            vals = [float(r["runtime_seconds"]) for r in result_rows if r["family"] == family and r["algorithm"] == algo]
            if vals:
                runtime_hist_has_data = True
                ax.hist(vals, bins=10, alpha=0.5, label=algo_labels[algo], color=algo_colors[algo])
        ax.set_xlabel("runtime (s)")
        ax.set_ylabel("count")
        ax.grid(alpha=0.2)
        ax.legend()
    fig.tight_layout()
    runtime_hist_path = os.path.join(outdir, "hist_runtime_seconds.png")
    fig.savefig(runtime_hist_path, dpi=200)
    plt.close(fig)
    if runtime_hist_has_data:
        saved.append(runtime_hist_path)

    return saved


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n", type=int, required=True)
    ap.add_argument("--instances", type=int, default=10, help="Instances per family.")
    ap.add_argument("--families", type=str, default="random,maxcut,mis")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", type=str, default="diagnostics_local/compare")

    # Shared QUBO generation config.
    ap.add_argument("--random-low", type=int, default=-5)
    ap.add_argument("--random-high", type=int, default=5)
    ap.add_argument("--random-density", type=float, default=1.0)
    ap.add_argument("--graph-p", type=float, default=0.3)
    ap.add_argument("--maxcut-weight-low", type=int, default=1)
    ap.add_argument("--maxcut-weight-high", type=int, default=1)
    ap.add_argument("--mis-node-weight-low", type=int, default=1)
    ap.add_argument("--mis-node-weight-high", type=int, default=1)
    ap.add_argument("--mis-lambda", type=float, default=2.0)

    # QA config.
    ap.add_argument("--delta-t", type=float, default=0.25)
    ap.add_argument("--t-max", type=float, default=10.0)
    ap.add_argument("--shots", type=int, default=128)
    ap.add_argument("--trotter-order", type=int, default=2, choices=[1, 2])
    ap.add_argument("--scale-dynamics", type=str, default="maxabs", choices=["none", "maxabs"])
    ap.add_argument("--aer-method", type=str, default="statevector", choices=["auto", "statevector", "matrix_product_state"])
    ap.add_argument("--mps-auto-threshold", type=int, default=26)
    ap.add_argument("--mps-max-bond-dimension", type=int, default=256)
    ap.add_argument("--mps-truncation-threshold", type=float, default=1e-12)
    ap.add_argument(
        "--mps-sample-measure-algorithm",
        type=str,
        default="mps_apply_measure",
        choices=["mps_apply_measure", "mps_probabilities"],
    )
    ap.add_argument(
        "--mps-swap-direction",
        type=str,
        default="mps_swap_left",
        choices=["mps_swap_left", "mps_swap_right"],
    )
    ap.add_argument("--mps-omp-threads", type=int, default=1)
    ap.add_argument("--mps-log-data", action="store_true")
    ap.add_argument("--mps-lapack", action="store_true")

    # SA config.
    ap.add_argument("--sa-reads", type=int, default=128)
    ap.add_argument("--sa-sweep-checkpoints", type=str, default="64,256,1024,4096")

    # PRR config.
    ap.add_argument("--prr-total-time", type=float, default=3.5)
    ap.add_argument("--prr-segments", type=int, default=8)
    ap.add_argument("--prr-omega-max", type=float, default=5.0)
    ap.add_argument("--prr-g-min", type=float, default=-8.0)
    ap.add_argument("--prr-g-max", type=float, default=4.0)
    ap.add_argument("--prr-g-start", type=float, default=-3.0)
    ap.add_argument("--prr-maxiter-bfgs", type=int, default=120)
    ap.add_argument("--prr-maxiter-nm", type=int, default=180)

    # Opt reference.
    ap.add_argument("--opt-ref", type=str, default="exact", choices=["exact", "none"])
    ap.add_argument("--exact-max-n", type=int, default=16)
    ap.add_argument("--plot-grid-points", type=int, default=101)
    ap.add_argument("--no-plots", action="store_true")

    args = ap.parse_args()

    if int(args.n) <= 0:
        raise SystemExit("--n must be positive.")
    if int(args.instances) <= 0:
        raise SystemExit("--instances must be positive.")
    if int(args.shots) <= 0:
        raise SystemExit("--shots must be positive.")
    if float(args.delta_t) <= 0.0:
        raise SystemExit("--delta-t must be positive.")
    if float(args.t_max) < 0.0:
        raise SystemExit("--t-max must be nonnegative.")
    if int(args.random_low) > int(args.random_high):
        raise SystemExit("--random-low must be <= --random-high.")
    if not (0.0 <= float(args.random_density) <= 1.0):
        raise SystemExit("--random-density must be in [0,1].")
    if not (0.0 <= float(args.graph_p) <= 1.0):
        raise SystemExit("--graph-p must be in [0,1].")
    if int(args.maxcut_weight_low) <= 0 or int(args.maxcut_weight_low) > int(args.maxcut_weight_high):
        raise SystemExit("--maxcut weights require positive low <= high.")
    if int(args.mis_node_weight_low) <= 0 or int(args.mis_node_weight_low) > int(args.mis_node_weight_high):
        raise SystemExit("--mis node weights require positive low <= high.")
    if float(args.mis_lambda) <= max(1.0, float(args.mis_node_weight_high)):
        raise SystemExit("--mis-lambda must be > max(1.0, --mis-node-weight-high).")
    if int(args.sa_reads) <= 0:
        raise SystemExit("--sa-reads must be positive.")
    if int(args.prr_segments) < 3:
        raise SystemExit("--prr-segments must be >= 3.")
    if float(args.prr_total_time) <= 0.0:
        raise SystemExit("--prr-total-time must be positive.")

    families = resolve_family_order(args.families)
    sa_checkpoints = parse_sa_sweep_checkpoints(args.sa_sweep_checkpoints)
    k_max = int(round(float(args.t_max) / float(args.delta_t)))
    qa_budget = (
        np.zeros(k_max + 1, dtype=float)
        if k_max <= 0
        else np.arange(k_max + 1, dtype=float) / float(k_max)
    )
    sa_budget = np.asarray(sa_checkpoints, dtype=float) / float(max(sa_checkpoints))

    os.makedirs(args.outdir, exist_ok=True)
    generation_config = {
        "random_low": int(args.random_low),
        "random_high": int(args.random_high),
        "random_density": float(args.random_density),
        "graph_p": float(args.graph_p),
        "maxcut_weight_low": int(args.maxcut_weight_low),
        "maxcut_weight_high": int(args.maxcut_weight_high),
        "mis_node_weight_low": int(args.mis_node_weight_low),
        "mis_node_weight_high": int(args.mis_node_weight_high),
        "mis_lambda": float(args.mis_lambda),
    }

    bank = build_instance_bank(
        n=int(args.n),
        instances_per_family=int(args.instances),
        seed_base=int(args.seed),
        random_low=int(args.random_low),
        random_high=int(args.random_high),
        random_density=float(args.random_density),
        graph_p=float(args.graph_p),
        maxcut_weight_low=int(args.maxcut_weight_low),
        maxcut_weight_high=int(args.maxcut_weight_high),
        mis_lambda=float(args.mis_lambda),
        mis_node_weight_low=int(args.mis_node_weight_low),
        mis_node_weight_high=int(args.mis_node_weight_high),
    )

    manifest_path = write_instance_bank_manifest(
        outdir=args.outdir,
        n=int(args.n),
        instances_per_family=int(args.instances),
        seed_base=int(args.seed),
        generation_config=generation_config,
        bank=bank,
    )
    print(f"Wrote: {manifest_path}")

    backend = build_backend(args, n=int(args.n))
    print(f"[compare] [backend] Aer method={backend.options.method}")

    result_rows: List[Dict[str, Any]] = []
    curve_rows: List[Dict[str, Any]] = []
    exact_reference_enabled = bool(str(args.opt_ref) == "exact" and int(args.n) <= int(args.exact_max_n))
    sa_engine_name = "unknown"

    for family in families:
        print(f"\n=== compare n={args.n} family={family} ===")
        for entry in bank[family]:
            bqm = entry.bqm
            model_eval = bqm_to_ising_arrays(bqm, n=int(args.n))
            model_dyn, _dyn_scale = scale_ising_for_dynamics(model_eval, mode=str(args.scale_dynamics))

            opt_energy: Optional[float] = None
            if exact_reference_enabled:
                exact = dimod.ExactSolver().sample(bqm)
                opt_energy = float(exact.first.energy)

            # QA run.
            qa_start = time.perf_counter()
            qa_curve_raw, _qa_bits, _qa_expect = run_instance_sweep(
                backend=backend,
                model_eval=model_eval,
                model_dyn=model_dyn,
                k_max=k_max,
                delta_t=float(args.delta_t),
                trotter_order=int(args.trotter_order),
                shots=int(args.shots),
                seed=int(entry.seed) + 11,
                seed_transpiler=int(entry.seed) + 101,
                cache_mode="support",
                estimator=None,
                use_transpile_cache=False,
            )
            qa_runtime = float(time.perf_counter() - qa_start)
            qa_curve = np.minimum.accumulate(qa_curve_raw)
            qa_final = float(qa_curve[-1])

            # SA run.
            sa_curve, sa_runtime_curve, sa_engine = run_sa_convergence_trace(
                bqm=bqm,
                num_reads=int(args.sa_reads),
                sweep_checkpoints=sa_checkpoints,
                seed=int(entry.seed) + 22,
            )
            sa_engine_name = str(sa_engine)
            sa_final = float(sa_curve[-1])
            sa_runtime = float(sa_runtime_curve[-1])

            # PRR run.
            prr_summary = run_prr_local_detuning_for_bqm(
                bqm=bqm,
                problem=family,
                meta=entry.meta,
                n=int(args.n),
                graph_p=float(args.graph_p),
                seed=int(entry.seed) + 33,
                total_time=float(args.prr_total_time),
                segments=int(args.prr_segments),
                omega_max=float(args.prr_omega_max),
                g_min=float(args.prr_g_min),
                g_max=float(args.prr_g_max),
                g_start=float(args.prr_g_start),
                maxiter_bfgs=int(args.prr_maxiter_bfgs),
                maxiter_nm=int(args.prr_maxiter_nm),
                exact_max_n=int(args.exact_max_n),
            )
            prr_final = float(prr_summary["sampled_best_energy"])
            prr_runtime = float(prr_summary["optimization"]["stage"]["total_walltime_seconds"])
            prr_curve = np.asarray(prr_summary["optimization"]["trace_best_objective_so_far"], dtype=float)
            prr_budget = np.asarray(prr_summary["optimization"]["trace_budget_fraction"], dtype=float)

            for idx in range(int(len(qa_curve))):
                curve_rows.append(
                    {
                        "family": family,
                        "instance_id": int(entry.instance_id),
                        "algorithm": "qa",
                        "point_index": int(idx),
                        "budget_fraction": float(qa_budget[idx]),
                        "best_energy_so_far": float(qa_curve[idx]),
                        "energy_definition": "best_sampled_energy",
                    }
                )
            for idx in range(int(len(sa_curve))):
                curve_rows.append(
                    {
                        "family": family,
                        "instance_id": int(entry.instance_id),
                        "algorithm": "sa",
                        "point_index": int(idx),
                        "budget_fraction": float(sa_budget[idx]),
                        "best_energy_so_far": float(sa_curve[idx]),
                        "energy_definition": "best_sampled_energy",
                    }
                )
            for idx in range(int(len(prr_curve))):
                curve_rows.append(
                    {
                        "family": family,
                        "instance_id": int(entry.instance_id),
                        "algorithm": "prr",
                        "point_index": int(idx),
                        "budget_fraction": float(prr_budget[idx]),
                        "best_energy_so_far": float(prr_curve[idx]),
                        "energy_definition": "best_expectation_energy",
                    }
                )

            per_algo = [
                ("qa", qa_final, qa_runtime, "digital_qa"),
                ("sa", sa_final, sa_runtime, f"sa_{sa_engine_name}"),
                ("prr", prr_final, prr_runtime, str(prr_summary["protocol_mode"])),
            ]
            for algo, final_energy, runtime_seconds, mode_tag in per_algo:
                ratio: Optional[float] = None
                ratio_error: Optional[float] = None
                gap_to_opt: Optional[float] = None
                success: Optional[bool] = None
                if opt_energy is not None:
                    ratio = float(approximation_ratio_minimization(obtained=final_energy, optimal=float(opt_energy)))
                    ratio_error = float(1.0 - ratio)
                    gap_to_opt = float(final_energy - float(opt_energy))
                    success = bool(final_energy <= float(opt_energy) + 1e-9)

                result_rows.append(
                    {
                        "family": family,
                        "instance_id": int(entry.instance_id),
                        "bqm_digest": str(entry.bqm_digest),
                        "algorithm": algo,
                        "mode_tag": mode_tag,
                        "final_energy": float(final_energy),
                        "opt_energy": float(opt_energy) if opt_energy is not None else None,
                        "gap_to_opt": gap_to_opt,
                        "approx_ratio": ratio,
                        "approx_ratio_error": ratio_error,
                        "success_vs_opt": success,
                        "runtime_seconds": float(runtime_seconds),
                    }
                )

            if (int(entry.instance_id) + 1) % max(1, int(args.instances) // 5) == 0:
                print(f"  instance {int(entry.instance_id) + 1}/{int(args.instances)} done")

    results_csv_path = os.path.join(args.outdir, "comparison_results.csv")
    with open(results_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "family",
                "instance_id",
                "bqm_digest",
                "algorithm",
                "mode_tag",
                "final_energy",
                "opt_energy",
                "gap_to_opt",
                "approx_ratio",
                "approx_ratio_error",
                "success_vs_opt",
                "runtime_seconds",
            ]
        )
        for row in result_rows:
            w.writerow(
                [
                    row["family"],
                    row["instance_id"],
                    row["bqm_digest"],
                    row["algorithm"],
                    row["mode_tag"],
                    row["final_energy"],
                    row["opt_energy"],
                    row["gap_to_opt"],
                    row["approx_ratio"],
                    row["approx_ratio_error"],
                    int(row["success_vs_opt"]) if row["success_vs_opt"] is not None else "",
                    row["runtime_seconds"],
                ]
            )
    print(f"Wrote: {results_csv_path}")

    curves_csv_path = os.path.join(args.outdir, "comparison_curves.csv")
    with open(curves_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "family",
                "instance_id",
                "algorithm",
                "point_index",
                "budget_fraction",
                "best_energy_so_far",
                "energy_definition",
            ]
        )
        for row in curve_rows:
            w.writerow(
                [
                    row["family"],
                    row["instance_id"],
                    row["algorithm"],
                    row["point_index"],
                    row["budget_fraction"],
                    row["best_energy_so_far"],
                    row["energy_definition"],
                ]
            )
    print(f"Wrote: {curves_csv_path}")

    summary = {
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "n": int(args.n),
        "instances_per_family": int(args.instances),
        "families": families,
        "opt_ref": str(args.opt_ref),
        "exact_reference_enabled": bool(exact_reference_enabled),
        "exact_max_n": int(args.exact_max_n),
        "aer_method": str(backend.options.method),
        "qa": {
            "delta_t": float(args.delta_t),
            "t_max": float(args.t_max),
            "k_max": int(k_max),
            "shots": int(args.shots),
            "trotter_order": int(args.trotter_order),
            "scale_dynamics": str(args.scale_dynamics),
        },
        "sa": {
            "reads": int(args.sa_reads),
            "sweep_checkpoints": [int(v) for v in sa_checkpoints],
            "engine": str(sa_engine_name),
        },
        "prr": {
            "total_time": float(args.prr_total_time),
            "segments": int(args.prr_segments),
            "omega_max": float(args.prr_omega_max),
            "g_min": float(args.prr_g_min),
            "g_max": float(args.prr_g_max),
            "g_start": float(args.prr_g_start),
            "maxiter_bfgs": int(args.prr_maxiter_bfgs),
            "maxiter_nm": int(args.prr_maxiter_nm),
        },
        "generation_config": generation_config,
        "aggregates": build_summary_aggregates(result_rows=result_rows, families=families),
        "note": {
            "ratio_error_definition": "approx_ratio_error = 1 - R",
            "energy_curve_definitions": {
                "qa": "best_sampled_energy",
                "sa": "best_sampled_energy",
                "prr": "best_expectation_energy",
            },
        },
    }

    plot_paths: List[str] = []
    if not bool(args.no_plots):
        plot_paths = _write_plots(
            outdir=args.outdir,
            families=families,
            result_rows=result_rows,
            curve_rows=curve_rows,
            plot_grid_points=int(args.plot_grid_points),
        )
        for path in plot_paths:
            print(f"Wrote: {path}")
    summary["plots"] = {
        "enabled": not bool(args.no_plots),
        "generated": [os.path.basename(path) for path in plot_paths],
    }

    summary_json_path = os.path.join(args.outdir, "comparison_summary.json")
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote: {summary_json_path}")


if __name__ == "__main__":
    main()
