#!/usr/bin/env python3
"""
QA/adiabatic (digital) annealing benchmark for three QUBO families:
  - Random symmetric Q in [-5, 5] (optionally sparse), objective x^T Q x
  - MaxCut QUBO on an Erdos-Renyi graph
  - MIS QUBO on an Erdos-Renyi graph, objective -sum x_i + lambda * sum_{(i,j)} x_i x_j

We simulate adiabatic evolution under:
  H(s) = (1-s) * H_driver + s * H_problem
  H_driver = - sum_i X_i
  H_problem = sum_i h_i Z_i + sum_{i<j} J_ij Z_i Z_j

Dynamics are approximated by first- or second-order Trotterization with time slices of width delta_t.

Convergence metric:
  delta_t fixed (default 0.25), start at t=0
  Sweep K = 0..Kmax where total time t = K * delta_t
  For each K, run the circuit, compute best energy among shots.
  steps_to_opt is the smallest K where best-so-far energy equals a reference optimum:
    - opt_ref="qa_best": reference optimum is the best energy observed over the sweep
    - opt_ref="exact": for small n only, compute exact optimum using dimod.ExactSolver

Outputs:
  - convergence_energy.png : one plot with 3 curves (Random/MaxCut/MIS), median best-so-far energy vs time
  - success_prob.png : one plot with 3 curves, success probability vs time
  - expectation_energy.png : optional estimator diagnostics plot (when --estimator-diagnostics is set)
  - steps_boxplot.png : distribution of steps_to_opt
  - results.csv : per-instance summary
  - summary.json : statistics (Mann-Whitney or permutation one-sided + Holm correction), effect sizes
  - scan_summary.csv : when --n-list is used, one row per n with medians, p-values, and deltas

Notes:
  - For large n, use Aer method="matrix_product_state" (MPS). This can scale to more qubits if
    entanglement remains moderate; otherwise it can still blow up.
  - Transpile cache supports `support` and `full` modes; low hit-rate runs can auto-disable cache.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import deque
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import dimod
import networkx as nx
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import EstimatorV2, SamplerV2

try:
    from scipy.stats import mannwhitneyu
except Exception:  # pragma: no cover
    mannwhitneyu = None


def _rand_symmetric_int_matrix(
    rng: np.random.Generator,
    n: int,
    low: int = -5,
    high: int = 5,
    density: float = 1.0,
) -> np.ndarray:
    """Random symmetric integer matrix Q with entries in [low, high]."""
    q = rng.integers(low, high + 1, size=(n, n), dtype=np.int32)
    q = np.triu(q)
    q = q + q.T - np.diag(np.diag(q))
    if density < 1.0:
        mask = rng.random((n, n)) < density
        mask = np.triu(mask, 1)
        q_off = np.triu(q, 1)
        q_off[~mask] = 0
        q = np.diag(np.diag(q)) + q_off + q_off.T
    return q


def qubo_dict_from_symmetric_q_for_x_t_q_x(q_sym: np.ndarray) -> Dict[Tuple[int, int], float]:
    """
    Build an upper-triangular QUBO dict for objective x^T Qsym x with x in {0,1}^n.

    dimod expects QUBO energy: sum_i Q_ii x_i + sum_{i<j} Q_ij x_i x_j
    For symmetric Qsym:
        x^T Qsym x = sum_i Qsym_ii x_i + 2 * sum_{i<j} Qsym_ij x_i x_j
    """
    n = q_sym.shape[0]
    q: Dict[Tuple[int, int], float] = {}
    for i in range(n):
        q[(i, i)] = float(q_sym[i, i])
    for i in range(n):
        for j in range(i + 1, n):
            v = q_sym[i, j]
            if v != 0:
                q[(i, j)] = float(2 * v)
    return q


def make_random_qubo_bqm(
    rng: np.random.Generator,
    n: int,
    low: int,
    high: int,
    density: float,
) -> Tuple[dimod.BinaryQuadraticModel, Dict[str, Any]]:
    q_sym = _rand_symmetric_int_matrix(rng, n, low=low, high=high, density=density)
    qubo = qubo_dict_from_symmetric_q_for_x_t_q_x(q_sym)
    bqm = dimod.BinaryQuadraticModel.from_qubo(qubo)
    meta = {"Qsym": q_sym.tolist(), "density": density}
    return bqm, meta


def _average_edge_span(g: nx.Graph) -> float:
    if g.number_of_edges() == 0:
        return 0.0
    return float(np.mean([abs(int(i) - int(j)) for (i, j) in g.edges()]))


def _bfs_high_degree_order(g: nx.Graph) -> List[int]:
    remaining = set(g.nodes())
    order: List[int] = []

    while remaining:
        start = max(remaining, key=lambda node: (g.degree[node], -int(node)))
        q: deque[int] = deque([int(start)])
        remaining.remove(start)
        while q:
            u = q.popleft()
            order.append(u)
            nbrs = [int(v) for v in g.neighbors(u) if v in remaining]
            nbrs.sort(key=lambda node: (g.degree[node], -int(node)), reverse=True)
            for v in nbrs:
                remaining.remove(v)
                q.append(v)
    return order


def relabel_graph_bfs_high_degree(g: nx.Graph) -> Tuple[nx.Graph, Dict[int, int], Dict[str, float]]:
    """
    Relabel graph nodes so BFS traversal from highest-degree seeds maps to nearby indices.
    This tends to reduce edge distances in the linear qubit ordering used by MPS.
    """
    order_old = _bfs_high_degree_order(g)
    mapping_old_to_new = {old: new for new, old in enumerate(order_old)}
    g_relabel = nx.relabel_nodes(g, mapping_old_to_new, copy=True)
    metrics = {
        "edge_span_before": _average_edge_span(g),
        "edge_span_after": _average_edge_span(g_relabel),
    }
    return g_relabel, mapping_old_to_new, metrics


def make_maxcut_qubo_bqm(
    rng: np.random.Generator,
    n: int,
    p: float,
    weight_low: int = 1,
    weight_high: int = 1,
) -> Tuple[dimod.BinaryQuadraticModel, Dict[str, Any]]:
    """
    MaxCut objective (maximize):
        sum_{(i,j)} w_ij * (x_i XOR x_j)
    Turn into minimization by negating.
    """
    g_raw = nx.erdos_renyi_graph(n, p, seed=int(rng.integers(0, 2**31 - 1)))
    g, mapping_old_to_new, metrics = relabel_graph_bfs_high_degree(g_raw)
    lin = {i: 0.0 for i in range(n)}
    quad: Dict[Tuple[int, int], float] = {}
    for (i, j) in g.edges():
        w = float(rng.integers(weight_low, weight_high + 1))
        lin[i] += -w
        lin[j] += -w
        ii, jj = (i, j) if i < j else (j, i)
        quad[(ii, jj)] = quad.get((ii, jj), 0.0) + 2.0 * w
    bqm = dimod.BinaryQuadraticModel(lin, quad, 0.0, vartype=dimod.BINARY)
    meta = {
        "p": p,
        "m": g.number_of_edges(),
        "node_order_old_to_new": mapping_old_to_new,
        **metrics,
    }
    return bqm, meta


def make_mis_qubo_bqm(
    rng: np.random.Generator,
    n: int,
    p: float,
    penalty_lambda: float,
) -> Tuple[dimod.BinaryQuadraticModel, Dict[str, Any]]:
    """
    MIS QUBO (minimization):
        minimize -sum_i x_i + lambda * sum_{(i,j) in E} x_i x_j
    """
    g_raw = nx.erdos_renyi_graph(n, p, seed=int(rng.integers(0, 2**31 - 1)))
    g, mapping_old_to_new, metrics = relabel_graph_bfs_high_degree(g_raw)
    lin = {i: -1.0 for i in range(n)}
    quad: Dict[Tuple[int, int], float] = {}
    for (i, j) in g.edges():
        ii, jj = (i, j) if i < j else (j, i)
        quad[(ii, jj)] = quad.get((ii, jj), 0.0) + float(penalty_lambda)
    bqm = dimod.BinaryQuadraticModel(lin, quad, 0.0, vartype=dimod.BINARY)
    meta = {
        "p": p,
        "m": g.number_of_edges(),
        "lambda": float(penalty_lambda),
        "node_order_old_to_new": mapping_old_to_new,
        **metrics,
    }
    return bqm, meta


@dataclass(frozen=True)
class IsingModel:
    n: int
    h: np.ndarray
    j_terms: List[Tuple[int, int, float]]
    offset: float


def bqm_to_ising_arrays(bqm: dimod.BinaryQuadraticModel, n: int) -> IsingModel:
    h_dict, j_dict, offset = bqm.to_ising()
    h = np.zeros(n, dtype=float)
    for i, v in h_dict.items():
        h[int(i)] = float(v)
    j_list: List[Tuple[int, int, float]] = []
    for (i, j), v in j_dict.items():
        ii = int(i)
        jj = int(j)
        if ii == jj:
            continue
        if ii < jj:
            j_list.append((ii, jj, float(v)))
        else:
            j_list.append((jj, ii, float(v)))
    return IsingModel(n=n, h=h, j_terms=j_list, offset=float(offset))


def scale_ising_for_dynamics(model: IsingModel, mode: str = "maxabs") -> Tuple[IsingModel, float]:
    """
    Scale h,J for dynamics only (angles). Positive scaling preserves argmin state.
    """
    if mode == "none":
        return model, 1.0
    if mode != "maxabs":
        raise ValueError(f"Unknown scale mode: {mode}")

    max_h = float(np.max(np.abs(model.h))) if model.h.size else 0.0
    max_j = max((abs(v) for _, _, v in model.j_terms), default=0.0)
    denom = max(max_h, max_j, 1e-12)
    scale = 1.0 / denom
    h2 = model.h * scale
    j2 = [(i, j, v * scale) for (i, j, v) in model.j_terms]
    return IsingModel(n=model.n, h=h2, j_terms=j2, offset=model.offset * scale), scale


@dataclass
class TranspiledTemplate:
    h_indices: Tuple[int, ...]
    j_pairs: Tuple[Tuple[int, int], ...]
    h_params: Dict[int, Parameter]
    j_params: Dict[Tuple[int, int], Parameter]
    circuits: List[QuantumCircuit]


_FULL_SUPPORT_CACHE: Dict[int, Tuple[Tuple[int, ...], Tuple[Tuple[int, int], ...]]] = {}


def support_signature_sparse(model_dyn: IsingModel, tol: float = 1e-15) -> Tuple[Tuple[int, ...], Tuple[Tuple[int, int], ...]]:
    h_indices = tuple(i for i in range(model_dyn.n) if abs(float(model_dyn.h[i])) > tol)
    j_pairs = sorted(
        (
            (int(i), int(j)) if int(i) < int(j) else (int(j), int(i))
            for (i, j, v) in model_dyn.j_terms
            if abs(float(v)) > tol
        )
    )
    return h_indices, tuple(j_pairs)


def support_signature_full(n: int) -> Tuple[Tuple[int, ...], Tuple[Tuple[int, int], ...]]:
    if n in _FULL_SUPPORT_CACHE:
        return _FULL_SUPPORT_CACHE[n]
    h_indices = tuple(range(int(n)))
    j_pairs = tuple((i, j) for i in range(int(n)) for j in range(i + 1, int(n)))
    _FULL_SUPPORT_CACHE[n] = (h_indices, j_pairs)
    return _FULL_SUPPORT_CACHE[n]


def resolve_support_signature(model_dyn: IsingModel, cache_mode: str) -> Tuple[Tuple[int, ...], Tuple[Tuple[int, int], ...]]:
    if cache_mode == "support":
        return support_signature_sparse(model_dyn)
    if cache_mode == "full":
        return support_signature_full(model_dyn.n)
    raise ValueError(f"Unknown transpile cache mode: {cache_mode}")


def bump_cache_stat(cache_stats: Optional[Dict[str, int]], key: str) -> None:
    if cache_stats is None:
        return
    cache_stats[key] = cache_stats.get(key, 0) + 1


def build_parameterized_anneal_circuits(
    n: int,
    h_indices: Tuple[int, ...],
    j_pairs: Tuple[Tuple[int, int], ...],
    k_max: int,
    delta_t: float,
    trotter_order: int,
    include_measurements: bool,
) -> Tuple[List[QuantumCircuit], Dict[int, Parameter], Dict[Tuple[int, int], Parameter]]:
    h_params = {i: Parameter(f"h_{i}") for i in h_indices}
    j_params = {(i, j): Parameter(f"j_{i}_{j}") for (i, j) in j_pairs}

    circuits: List[QuantumCircuit] = []
    for k in range(0, k_max + 1):
        if include_measurements:
            qc = QuantumCircuit(n, n)
        else:
            qc = QuantumCircuit(n)
        qc.h(range(n))

        if k > 0:
            for step in range(k):
                s = (step + 0.5) / float(k)
                t_d = delta_t * (1.0 - s)
                t_p = delta_t * s

                if trotter_order == 1:
                    _apply_driver_layer(qc, n, t_d)
                    for i in h_indices:
                        qc.rz(2.0 * t_p * h_params[i], i)
                    for (i, j) in j_pairs:
                        qc.rzz(2.0 * t_p * j_params[(i, j)], i, j)
                elif trotter_order == 2:
                    _apply_driver_layer(qc, n, 0.5 * t_d)
                    for i in h_indices:
                        qc.rz(2.0 * t_p * h_params[i], i)
                    for (i, j) in j_pairs:
                        qc.rzz(2.0 * t_p * j_params[(i, j)], i, j)
                    _apply_driver_layer(qc, n, 0.5 * t_d)
                else:
                    raise ValueError("trotter_order must be 1 or 2")

        if include_measurements:
            qc.measure(range(n), range(n))
        qc.metadata = {"k": k}
        circuits.append(qc)

    return circuits, h_params, j_params


def get_transpiled_template(
    backend: AerSimulator,
    model_dyn: IsingModel,
    k_max: int,
    delta_t: float,
    trotter_order: int,
    seed_transpiler: int,
    include_measurements: bool,
    cache_mode: str,
    cache: Optional[Dict[Tuple[Any, ...], TranspiledTemplate]],
    cache_stats: Optional[Dict[str, int]],
) -> TranspiledTemplate:
    h_indices, j_pairs = resolve_support_signature(model_dyn, cache_mode=cache_mode)

    if cache_mode == "support":
        key = (
            model_dyn.n,
            cache_mode,
            h_indices,
            j_pairs,
            int(k_max),
            float(delta_t),
            int(trotter_order),
            bool(include_measurements),
        )
    else:
        key = (
            model_dyn.n,
            cache_mode,
            int(k_max),
            float(delta_t),
            int(trotter_order),
            bool(include_measurements),
        )

    if cache is not None and key in cache:
        bump_cache_stat(cache_stats, "measured_hits" if include_measurements else "unmeasured_hits")
        return cache[key]

    if cache is None:
        bump_cache_stat(cache_stats, "measured_bypassed" if include_measurements else "unmeasured_bypassed")
    else:
        bump_cache_stat(cache_stats, "measured_misses" if include_measurements else "unmeasured_misses")

    circuits, h_params, j_params = build_parameterized_anneal_circuits(
        n=model_dyn.n,
        h_indices=h_indices,
        j_pairs=j_pairs,
        k_max=k_max,
        delta_t=delta_t,
        trotter_order=trotter_order,
        include_measurements=include_measurements,
    )
    tqc = transpile(
        circuits,
        backend=backend,
        optimization_level=0,
        seed_transpiler=seed_transpiler,
        num_processes=1,
    )
    template = TranspiledTemplate(
        h_indices=h_indices,
        j_pairs=j_pairs,
        h_params=h_params,
        j_params=j_params,
        circuits=tqc,
    )
    if cache is not None:
        cache[key] = template
    return template


def bind_template_circuits(template: TranspiledTemplate, model_dyn: IsingModel) -> List[QuantumCircuit]:
    j_lookup = {(int(i), int(j)): float(v) for (i, j, v) in model_dyn.j_terms}
    bind_map: Dict[Parameter, float] = {}
    for i in template.h_indices:
        bind_map[template.h_params[i]] = float(model_dyn.h[i])
    for pair in template.j_pairs:
        bind_map[template.j_params[pair]] = float(j_lookup.get(pair, 0.0))
    return [qc.assign_parameters(bind_map, inplace=False, strict=False) for qc in template.circuits]


def ising_model_to_sparse_pauli_op(model_eval: IsingModel) -> SparsePauliOp:
    pauli_terms: List[Tuple[str, float]] = []
    n = model_eval.n

    if model_eval.offset != 0.0:
        pauli_terms.append(("I" * n, float(model_eval.offset)))

    for i, h_i in enumerate(model_eval.h):
        if h_i == 0.0:
            continue
        chars = ["I"] * n
        chars[n - 1 - i] = "Z"
        pauli_terms.append(("".join(chars), float(h_i)))

    for (i, j, jij) in model_eval.j_terms:
        if jij == 0.0:
            continue
        chars = ["I"] * n
        chars[n - 1 - i] = "Z"
        chars[n - 1 - j] = "Z"
        pauli_terms.append(("".join(chars), float(jij)))

    if not pauli_terms:
        pauli_terms = [("I" * n, 0.0)]
    return SparsePauliOp.from_list(pauli_terms)


def _apply_driver_layer(qc: QuantumCircuit, n: int, t: float) -> None:
    # H_driver = -sum X, so exp(-i t H_driver) = exp(+i t sum X) = RX(-2t) per qubit
    theta = -2.0 * t
    for q in range(n):
        qc.rx(theta, q)


def _apply_problem_layer(qc: QuantumCircuit, model_dyn: IsingModel, t: float) -> None:
    # exp(-i t h Z) -> RZ(2 t h), exp(-i t J ZZ) -> RZZ(2 t J)
    n = model_dyn.n
    for i in range(n):
        hi = model_dyn.h[i]
        if hi != 0.0:
            qc.rz(2.0 * t * hi, i)
    for (i, j, jij) in model_dyn.j_terms:
        if jij != 0.0:
            qc.rzz(2.0 * t * jij, i, j)


def build_anneal_circuit(
    model_dyn: IsingModel,
    num_steps: int,
    delta_t: float,
    trotter_order: int = 2,
) -> QuantumCircuit:
    """
    Build adiabatic evolution circuit with num_steps slices of width delta_t.
    Schedule: s_k = (k + 0.5) / num_steps (midpoint).
    """
    n = model_dyn.n
    qc = QuantumCircuit(n, n)
    qc.h(range(n))

    if num_steps <= 0:
        qc.measure(range(n), range(n))
        return qc

    for k in range(num_steps):
        s = (k + 0.5) / float(num_steps)
        t_d = delta_t * (1.0 - s)
        t_p = delta_t * s

        if trotter_order == 1:
            _apply_driver_layer(qc, n, t_d)
            _apply_problem_layer(qc, model_dyn, t_p)
        elif trotter_order == 2:
            _apply_driver_layer(qc, n, 0.5 * t_d)
            _apply_problem_layer(qc, model_dyn, t_p)
            _apply_driver_layer(qc, n, 0.5 * t_d)
        else:
            raise ValueError("trotter_order must be 1 or 2")

    qc.measure(range(n), range(n))
    return qc


def ising_energy_from_bitstring(bitstring_lsb0: str, model_eval: IsingModel) -> float:
    """
    Compute Ising energy for measured computational basis state.
    bitstring_lsb0 is qubit-0 first (LSB-first).
    """
    n = model_eval.n
    # dimod's to_ising uses x = (s + 1) / 2 for binary variables, so measured bits map as:
    # bit 0 -> s=-1, bit 1 -> s=+1.
    s = np.fromiter((-1.0 if ch == "0" else 1.0 for ch in bitstring_lsb0[:n]), dtype=float, count=n)
    e = model_eval.offset + float(np.dot(model_eval.h, s))
    for (i, j, jij) in model_eval.j_terms:
        e += jij * s[i] * s[j]
    return float(e)


def best_energy_from_counts(counts: Dict[str, int], model_eval: IsingModel) -> Tuple[float, str]:
    """
    Return (best_energy, bitstring_lsb0) for the best sample present in counts.
    Qiskit count keys are MSB-first by classical register order, so reverse for qubit-0 first.
    """
    best_e = float("inf")
    best_b = "0" * model_eval.n
    for b_msb, _freq in counts.items():
        b = b_msb[::-1]
        e = ising_energy_from_bitstring(b, model_eval)
        if e < best_e:
            best_e = e
            best_b = b
    return best_e, best_b


def counts_from_sampler_pub(pub_result: Any) -> Dict[str, int]:
    """
    Extract counts from a SamplerV2 pub result without assuming classical-register name.
    """
    for value in pub_result.data.values():
        if hasattr(value, "get_counts"):
            return value.get_counts()
    raise RuntimeError("SamplerV2 result did not contain a classical register with counts.")


def cliffs_delta_smaller_is_better(x: np.ndarray, y: np.ndarray) -> float:
    """
    Cliff's delta for "smaller is better":
      delta = P(x < y) - P(x > y)
    Positive means x tends to be smaller.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    n = x.size
    m = y.size
    less = 0
    greater = 0
    for xi in x:
        less += int(np.sum(xi < y))
        greater += int(np.sum(xi > y))
    return (less - greater) / float(n * m)


def holm_adjust(pvals: Dict[str, float]) -> Dict[str, float]:
    """Holm-Bonferroni adjusted p-values."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    adj: Dict[str, float] = {}
    running_max = 0.0
    for i, (name, p) in enumerate(items):
        mult = m - i
        p_adj = min(1.0, p * mult)
        running_max = max(running_max, p_adj)
        adj[name] = running_max
    return adj


def permutation_test_median_diff_less(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    iterations: int,
) -> float:
    """
    One-sided permutation test for median(x) - median(y) under alternative "x < y".
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    observed = float(np.median(x) - np.median(y))
    combined = np.concatenate([x, y])
    n_x = x.size

    num_less_equal = 0
    for _ in range(iterations):
        perm = rng.permutation(combined)
        diff = float(np.median(perm[:n_x]) - np.median(perm[n_x:]))
        if diff <= observed + 1e-12:
            num_less_equal += 1

    # Add-one smoothing to avoid exact zero p-values in finite permutations.
    return (num_less_equal + 1.0) / (iterations + 1.0)


def parse_n_list(n_single: Optional[int], n_list: str) -> List[int]:
    if n_list:
        values: List[int] = []
        for chunk in n_list.split(","):
            token = chunk.strip()
            if not token:
                continue
            if "-" in token:
                parts = token.split("-", 1)
                if len(parts) != 2:
                    raise ValueError(f"Invalid n-list token: {token}")
                lo = int(parts[0].strip())
                hi = int(parts[1].strip())
                if lo > hi:
                    raise ValueError(f"Invalid descending range in n-list: {token}")
                values.extend(range(lo, hi + 1))
            else:
                values.append(int(token))
    else:
        if n_single is None:
            raise ValueError("Specify either --n or --n-list.")
        values = [int(n_single)]

    if not values:
        raise ValueError("No n values were parsed.")
    if any(v <= 0 for v in values):
        raise ValueError("All n values must be positive integers.")
    return values


@dataclass
class InstanceResult:
    family: str
    instance_id: int
    n: int
    k_max: int
    delta_t: float
    opt_ref: str
    opt_energy: float
    steps_to_opt: int
    best_final_energy: float
    reached_opt: bool


def build_backend(args: argparse.Namespace, n: int) -> AerSimulator:
    method = args.aer_method
    if method == "auto":
        method = "matrix_product_state" if n >= args.mps_auto_threshold else "statevector"

    backend_kwargs: Dict[str, Any] = {"method": method}

    if method == "matrix_product_state":
        backend_kwargs.update(
            {
                "matrix_product_state_max_bond_dimension": args.mps_max_bond_dimension,
                "matrix_product_state_truncation_threshold": args.mps_truncation_threshold,
                "mps_sample_measure_algorithm": args.mps_sample_measure_algorithm,
                "mps_log_data": bool(args.mps_log_data),
                "mps_swap_direction": args.mps_swap_direction,
                "mps_omp_threads": int(args.mps_omp_threads),
                "mps_lapack": bool(args.mps_lapack),
            }
        )

    return AerSimulator(**backend_kwargs)


def run_instance_sweep(
    backend: AerSimulator,
    model_eval: IsingModel,
    model_dyn: IsingModel,
    k_max: int,
    delta_t: float,
    trotter_order: int,
    shots: int,
    seed: int,
    seed_transpiler: int,
    cache_mode: str,
    estimator: Optional[EstimatorV2] = None,
    estimator_precision: Optional[float] = None,
    use_transpile_cache: bool = True,
    measured_template_cache: Optional[Dict[Tuple[Any, ...], TranspiledTemplate]] = None,
    unmeasured_template_cache: Optional[Dict[Tuple[Any, ...], TranspiledTemplate]] = None,
    cache_stats: Optional[Dict[str, int]] = None,
) -> Tuple[np.ndarray, List[str], Optional[np.ndarray]]:
    """
    Returns:
      energies[k] = best energy among shots at step-count k
      best_bits[k] = bitstring (qubit0-first) achieving energies[k]
    """
    measured_template = get_transpiled_template(
        backend=backend,
        model_dyn=model_dyn,
        k_max=k_max,
        delta_t=delta_t,
        trotter_order=trotter_order,
        seed_transpiler=seed_transpiler,
        include_measurements=True,
        cache_mode=cache_mode,
        cache=measured_template_cache if use_transpile_cache else None,
        cache_stats=cache_stats,
    )
    measured_circuits = bind_template_circuits(measured_template, model_dyn)

    sampler = SamplerV2.from_backend(backend, seed=seed)
    job = sampler.run(measured_circuits, shots=shots)
    result = job.result()

    energies = np.zeros(k_max + 1, dtype=float)
    best_bits: List[str] = []
    for idx in range(k_max + 1):
        counts = counts_from_sampler_pub(result[idx])
        e, b = best_energy_from_counts(counts, model_eval)
        energies[idx] = e
        best_bits.append(b)

    expectation_curve: Optional[np.ndarray] = None
    if estimator is not None:
        unmeasured_template = get_transpiled_template(
            backend=backend,
            model_dyn=model_dyn,
            k_max=k_max,
            delta_t=delta_t,
            trotter_order=trotter_order,
            seed_transpiler=seed_transpiler,
            include_measurements=False,
            cache_mode=cache_mode,
            cache=unmeasured_template_cache if use_transpile_cache else None,
            cache_stats=cache_stats,
        )
        unmeasured_circuits = bind_template_circuits(unmeasured_template, model_dyn)
        observable = ising_model_to_sparse_pauli_op(model_eval)
        pubs = [(qc, observable.apply_layout(qc.layout)) for qc in unmeasured_circuits]
        if estimator_precision is None:
            est_result = estimator.run(pubs).result()
        else:
            est_result = estimator.run(pubs, precision=estimator_precision).result()
        expectation_curve = np.asarray([float(est_result[idx].data.evs) for idx in range(k_max + 1)], dtype=float)

    return energies, best_bits, expectation_curve


def run_benchmark_for_n(
    args: argparse.Namespace,
    n: int,
    outdir: str,
    seed_base: int,
) -> Dict[str, Any]:
    if args.opt_ref == "exact" and n > int(args.exact_max_n):
        raise SystemExit(
            f"--opt-ref=exact requested but n={n} > exact-max-n={args.exact_max_n}. "
            "Use --opt-ref=qa_best for larger n."
        )

    rng = np.random.default_rng(seed_base)
    k_max = int(round(float(args.t_max) / float(args.delta_t)))
    t_axis = np.arange(k_max + 1, dtype=float) * float(args.delta_t)

    os.makedirs(outdir, exist_ok=True)
    backend = build_backend(args, n=n)
    print(f"[n={n}] [backend] Aer method={backend.options.method}")
    estimator = EstimatorV2.from_backend(backend) if bool(args.estimator_diagnostics) else None
    cache_requested = not bool(args.no_transpile_cache)
    use_transpile_cache = bool(cache_requested)
    cache_mode = str(args.transpile_cache_mode)
    cache_autodisable_enabled = bool(cache_requested and not bool(args.no_cache_autodisable))
    cache_autodisable_triggered = False
    cache_autodisable_reason: Optional[str] = None
    measured_template_cache: Dict[Tuple[Any, ...], TranspiledTemplate] = {}
    unmeasured_template_cache: Dict[Tuple[Any, ...], TranspiledTemplate] = {}
    cache_stats: Dict[str, int] = {
        "measured_hits": 0,
        "measured_misses": 0,
        "measured_bypassed": 0,
        "unmeasured_hits": 0,
        "unmeasured_misses": 0,
        "unmeasured_bypassed": 0,
    }

    families = ["random", "maxcut", "mis"]
    curves_best_so_far: Dict[str, List[np.ndarray]] = {f: [] for f in families}
    curves_raw_energy: Dict[str, List[np.ndarray]] = {f: [] for f in families}
    curves_expectation: Optional[Dict[str, List[np.ndarray]]] = (
        {f: [] for f in families} if bool(args.estimator_diagnostics) else None
    )
    opt_energies_by_family: Dict[str, List[float]] = {f: [] for f in families}
    steps_to_opt: Dict[str, List[int]] = {f: [] for f in families}
    all_rows: List[InstanceResult] = []

    for fam_idx, fam in enumerate(families):
        print(f"\n=== n={n} Family: {fam} ===")
        for inst_id in range(int(args.instances)):
            if fam == "random":
                bqm, _meta = make_random_qubo_bqm(
                    rng,
                    n,
                    low=int(args.random_low),
                    high=int(args.random_high),
                    density=float(args.random_density),
                )
            elif fam == "maxcut":
                bqm, _meta = make_maxcut_qubo_bqm(rng, n, p=float(args.graph_p))
            elif fam == "mis":
                bqm, _meta = make_mis_qubo_bqm(
                    rng,
                    n,
                    p=float(args.graph_p),
                    penalty_lambda=float(args.mis_lambda),
                )
            else:
                raise RuntimeError("Unknown family")

            model_eval = bqm_to_ising_arrays(bqm, n)
            model_dyn, _dyn_scale = scale_ising_for_dynamics(model_eval, mode=str(args.scale_dynamics))

            run_seed = int(seed_base) + fam_idx * 100000 + inst_id
            energies, _best_bits, exp_curve = run_instance_sweep(
                backend=backend,
                model_eval=model_eval,
                model_dyn=model_dyn,
                k_max=k_max,
                delta_t=float(args.delta_t),
                trotter_order=int(args.trotter_order),
                shots=int(args.shots),
                seed=run_seed,
                seed_transpiler=int(seed_base) + fam_idx * 100000,
                cache_mode=cache_mode,
                estimator=estimator,
                estimator_precision=args.estimator_precision,
                use_transpile_cache=use_transpile_cache,
                measured_template_cache=measured_template_cache,
                unmeasured_template_cache=unmeasured_template_cache,
                cache_stats=cache_stats,
            )

            best_so_far = np.minimum.accumulate(energies)
            curves_raw_energy[fam].append(energies)
            curves_best_so_far[fam].append(best_so_far)
            if curves_expectation is not None:
                if exp_curve is None:
                    raise RuntimeError("Estimator diagnostics enabled but no expectation curve was produced.")
                curves_expectation[fam].append(exp_curve)

            if args.opt_ref == "exact":
                exact = dimod.ExactSolver().sample(bqm)
                opt_energy = float(exact.first.energy)
            else:
                opt_energy = float(best_so_far[-1])
            opt_energies_by_family[fam].append(opt_energy)

            tol = 1e-9
            hit = np.where(best_so_far <= opt_energy + tol)[0]
            if hit.size > 0:
                k_star = int(hit[0])
                reached = True
            else:
                k_star = int(k_max + 1)
                reached = False

            steps_to_opt[fam].append(k_star)
            row = InstanceResult(
                family=fam,
                instance_id=inst_id,
                n=n,
                k_max=k_max,
                delta_t=float(args.delta_t),
                opt_ref=str(args.opt_ref),
                opt_energy=opt_energy,
                steps_to_opt=k_star,
                best_final_energy=float(best_so_far[-1]),
                reached_opt=reached,
            )
            all_rows.append(row)

            if cache_autodisable_enabled and use_transpile_cache:
                attempts = int(cache_stats.get("measured_hits", 0) + cache_stats.get("measured_misses", 0))
                if attempts >= int(args.cache_autodisable_min_attempts):
                    hit_rate = float(cache_stats.get("measured_hits", 0)) / float(max(1, attempts))
                    if hit_rate < float(args.cache_autodisable_min_hit_rate):
                        use_transpile_cache = False
                        cache_autodisable_enabled = False
                        cache_autodisable_triggered = True
                        cache_autodisable_reason = (
                            f"auto-disabled cache after {attempts} measured template requests: "
                            f"hit_rate={hit_rate:.3f} < threshold={float(args.cache_autodisable_min_hit_rate):.3f}"
                        )
                        print(f"[n={n}] [cache] {cache_autodisable_reason}")

            if (inst_id + 1) % max(1, (int(args.instances) // 10)) == 0:
                print(f"  instance {inst_id + 1}/{args.instances} done")

    def agg_curve(curves: List[np.ndarray], agg: str = "median") -> np.ndarray:
        m = np.vstack(curves)
        if agg == "mean":
            return np.mean(m, axis=0)
        return np.median(m, axis=0)

    def success_curve(curves: List[np.ndarray], opt_energies: List[float], tol: float = 1e-9) -> np.ndarray:
        m = np.vstack(curves)
        refs = np.asarray(opt_energies, dtype=float)[:, None]
        return np.mean(m <= refs + tol, axis=0)

    agg = "median"
    curve_random = agg_curve(curves_best_so_far["random"], agg=agg)
    curve_maxcut = agg_curve(curves_best_so_far["maxcut"], agg=agg)
    curve_mis = agg_curve(curves_best_so_far["mis"], agg=agg)

    success_cumulative_random = success_curve(curves_best_so_far["random"], opt_energies_by_family["random"])
    success_cumulative_maxcut = success_curve(curves_best_so_far["maxcut"], opt_energies_by_family["maxcut"])
    success_cumulative_mis = success_curve(curves_best_so_far["mis"], opt_energies_by_family["mis"])
    success_instant_random = success_curve(curves_raw_energy["random"], opt_energies_by_family["random"])
    success_instant_maxcut = success_curve(curves_raw_energy["maxcut"], opt_energies_by_family["maxcut"])
    success_instant_mis = success_curve(curves_raw_energy["mis"], opt_energies_by_family["mis"])
    exp_random: Optional[np.ndarray] = None
    exp_maxcut: Optional[np.ndarray] = None
    exp_mis: Optional[np.ndarray] = None
    if curves_expectation is not None:
        exp_random = agg_curve(curves_expectation["random"], agg=agg)
        exp_maxcut = agg_curve(curves_expectation["maxcut"], agg=agg)
        exp_mis = agg_curve(curves_expectation["mis"], agg=agg)

    steps_r = np.asarray(steps_to_opt["random"], dtype=float)
    steps_c = np.asarray(steps_to_opt["maxcut"], dtype=float)
    steps_m = np.asarray(steps_to_opt["mis"], dtype=float)
    easy_case_k_threshold = 0

    stats: Dict[str, Any] = {
        "n": n,
        "instances_per_family": int(args.instances),
        "delta_t": float(args.delta_t),
        "t_max": float(args.t_max),
        "k_max": k_max,
        "shots": int(args.shots),
        "opt_ref": str(args.opt_ref),
        "aer_method": str(backend.options.method),
        "scale_dynamics": str(args.scale_dynamics),
        "stats_method": str(args.stats_method),
        "transpile_cache": {
            "enabled": bool(cache_requested),
            "effective_enabled_final": bool(use_transpile_cache),
            "mode": str(cache_mode),
            "auto_disable": {
                "enabled": bool(cache_requested and not bool(args.no_cache_autodisable)),
                "min_attempts": int(args.cache_autodisable_min_attempts),
                "min_hit_rate": float(args.cache_autodisable_min_hit_rate),
                "triggered": bool(cache_autodisable_triggered),
                "reason": cache_autodisable_reason,
            },
            **cache_stats,
        },
        "estimator_diagnostics": {
            "enabled": bool(args.estimator_diagnostics),
            "precision": args.estimator_precision,
        },
    }

    stats["steps_summary"] = {
        "random": {"median": float(np.median(steps_r)), "mean": float(np.mean(steps_r))},
        "maxcut": {"median": float(np.median(steps_c)), "mean": float(np.mean(steps_c))},
        "mis": {"median": float(np.median(steps_m)), "mean": float(np.mean(steps_m))},
    }
    easy_case_rates = {
        "random": float(np.mean(steps_r <= easy_case_k_threshold)),
        "maxcut": float(np.mean(steps_c <= easy_case_k_threshold)),
        "mis": float(np.mean(steps_m <= easy_case_k_threshold)),
    }
    stats["easy_case"] = {
        "k_threshold": int(easy_case_k_threshold),
        "t_threshold": float(easy_case_k_threshold * float(args.delta_t)),
        "rate": easy_case_rates,
        "overall_mean_rate": float(np.mean(list(easy_case_rates.values()))),
    }
    stats["success_probability_cumulative_at_tmax"] = {
        "random": float(success_cumulative_random[-1]),
        "maxcut": float(success_cumulative_maxcut[-1]),
        "mis": float(success_cumulative_mis[-1]),
    }
    # Backward-compatible alias retained for existing consumers.
    stats["success_probability_at_tmax"] = dict(stats["success_probability_cumulative_at_tmax"])
    stats["success_probability_instantaneous_at_tmax"] = {
        "random": float(success_instant_random[-1]),
        "maxcut": float(success_instant_maxcut[-1]),
        "mis": float(success_instant_mis[-1]),
    }
    if exp_random is not None and exp_maxcut is not None and exp_mis is not None:
        stats["estimator_expectation_at_tmax"] = {
            "random": float(exp_random[-1]),
            "maxcut": float(exp_maxcut[-1]),
            "mis": float(exp_mis[-1]),
        }

    effects = {
        "delta_random_vs_maxcut": float(cliffs_delta_smaller_is_better(steps_r, steps_c)),
        "delta_random_vs_mis": float(cliffs_delta_smaller_is_better(steps_r, steps_m)),
    }
    stats["cliffs_delta_smaller_is_better"] = effects

    pvals: Optional[Dict[str, float]] = None
    if args.stats_method == "mw":
        if mannwhitneyu is None:
            stats["note"] = "SciPy not installed; cannot run Mann-Whitney tests."
        else:
            pvals = {
                "random_vs_maxcut": float(mannwhitneyu(steps_r, steps_c, alternative="less", method="auto").pvalue),
                "random_vs_mis": float(mannwhitneyu(steps_r, steps_m, alternative="less", method="auto").pvalue),
            }
            stats["pvals_method"] = "mannwhitneyu_one_sided_less"
    elif args.stats_method == "perm":
        perm_rng = np.random.default_rng(int(seed_base) + 777)
        pvals = {
            "random_vs_maxcut": float(
                permutation_test_median_diff_less(steps_r, steps_c, rng=perm_rng, iterations=int(args.perm_iterations))
            ),
            "random_vs_mis": float(
                permutation_test_median_diff_less(steps_r, steps_m, rng=perm_rng, iterations=int(args.perm_iterations))
            ),
        }
        stats["pvals_method"] = "permutation_median_diff_one_sided_less"
        stats["perm_iterations"] = int(args.perm_iterations)
    else:
        raise RuntimeError(f"Unknown stats method: {args.stats_method}")

    stats["pvals_one_sided_less"] = pvals
    stats["pvals_holm_adjusted"] = holm_adjust(pvals) if pvals is not None else None

    csv_path = os.path.join(outdir, "results.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "family",
                "instance_id",
                "n",
                "k_max",
                "delta_t",
                "opt_ref",
                "opt_energy",
                "steps_to_opt",
                "best_final_energy",
                "reached_opt",
            ]
        )
        for r in all_rows:
            w.writerow(
                [
                    r.family,
                    r.instance_id,
                    r.n,
                    r.k_max,
                    r.delta_t,
                    r.opt_ref,
                    r.opt_energy,
                    r.steps_to_opt,
                    r.best_final_energy,
                    int(r.reached_opt),
                ]
            )

    json_path = os.path.join(outdir, "summary.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print("\n=== Summary ===")
    print(json.dumps(stats["steps_summary"], indent=2))
    print("Easy-case rates (k <= 0):", stats["easy_case"]["rate"])
    if stats.get("pvals_one_sided_less") is not None:
        print(f"p-values ({stats['pvals_method']}):", stats["pvals_one_sided_less"])
        print("Holm-adjusted p-values:", stats["pvals_holm_adjusted"])
        print("Cliff's deltas:", stats["cliffs_delta_smaller_is_better"])
    print(f"\nWrote: {csv_path}")
    print(f"Wrote: {json_path}")

    if args.no_plots:
        return stats

    mpl_config_dir = os.path.join(outdir, ".mplconfig")
    os.makedirs(mpl_config_dir, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = mpl_config_dir

    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(t_axis, curve_random, label="Random QUBO")
    plt.plot(t_axis, curve_maxcut, label="MaxCut QUBO")
    plt.plot(t_axis, curve_mis, label="MIS QUBO")
    plt.xlabel("t (dimensionless), with delta_t = %.2f" % float(args.delta_t))
    plt.ylabel(f"{agg} best-so-far energy")
    plt.title(f"QA convergence (n={n}, N={args.instances}, shots={args.shots})")
    plt.legend()
    plt.tight_layout()
    out_plot = os.path.join(outdir, "convergence_energy.png")
    plt.savefig(out_plot, dpi=200)
    plt.close()
    print(f"Wrote: {out_plot}")

    plt.figure()
    plt.plot(t_axis, success_cumulative_random, label="Random QUBO")
    plt.plot(t_axis, success_cumulative_maxcut, label="MaxCut QUBO")
    plt.plot(t_axis, success_cumulative_mis, label="MIS QUBO")
    plt.xlabel("t (dimensionless), with delta_t = %.2f" % float(args.delta_t))
    plt.ylabel("cumulative success probability (reached by time t)")
    plt.ylim(0.0, 1.01)
    plt.title(f"Cumulative success probability vs time (n={n}, N={args.instances})")
    plt.legend()
    plt.tight_layout()
    out_plot_success = os.path.join(outdir, "success_prob.png")
    plt.savefig(out_plot_success, dpi=200)
    plt.close()
    print(f"Wrote: {out_plot_success}")

    plt.figure()
    plt.plot(t_axis, success_instant_random, label="Random QUBO")
    plt.plot(t_axis, success_instant_maxcut, label="MaxCut QUBO")
    plt.plot(t_axis, success_instant_mis, label="MIS QUBO")
    plt.xlabel("t (dimensionless), with delta_t = %.2f" % float(args.delta_t))
    plt.ylabel("instantaneous success probability (at time t)")
    plt.ylim(0.0, 1.01)
    plt.title(f"Instantaneous success probability vs time (n={n}, N={args.instances})")
    plt.legend()
    plt.tight_layout()
    out_plot_success_inst = os.path.join(outdir, "success_prob_instantaneous.png")
    plt.savefig(out_plot_success_inst, dpi=200)
    plt.close()
    print(f"Wrote: {out_plot_success_inst}")

    if exp_random is not None and exp_maxcut is not None and exp_mis is not None:
        plt.figure()
        plt.plot(t_axis, exp_random, label="Random QUBO")
        plt.plot(t_axis, exp_maxcut, label="MaxCut QUBO")
        plt.plot(t_axis, exp_mis, label="MIS QUBO")
        plt.xlabel("t (dimensionless), with delta_t = %.2f" % float(args.delta_t))
        plt.ylabel(f"{agg} estimator expectation energy")
        plt.title(f"Estimator diagnostics vs time (n={n}, N={args.instances})")
        plt.legend()
        plt.tight_layout()
        out_plot_expectation = os.path.join(outdir, "expectation_energy.png")
        plt.savefig(out_plot_expectation, dpi=200)
        plt.close()
        print(f"Wrote: {out_plot_expectation}")

    plt.figure()
    plt.boxplot([steps_r, steps_c, steps_m], tick_labels=["Random", "MaxCut", "MIS"])
    plt.ylabel("steps_to_opt (K), where t = K * delta_t")
    plt.title(f"Steps to optimum (opt_ref={args.opt_ref})")
    plt.tight_layout()
    out_plot2 = os.path.join(outdir, "steps_boxplot.png")
    plt.savefig(out_plot2, dpi=200)
    plt.close()
    print(f"Wrote: {out_plot2}")

    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--n", type=int, help="Number of qubits / QUBO dimension.")
    ap.add_argument("--n-list", type=str, default="", help="Comma-separated n values or ranges (e.g. 4,5,6 or 4-8).")
    ap.add_argument("--instances", type=int, default=100, help="Number of instances per family.")
    ap.add_argument("--delta-t", type=float, default=0.25, help="Time step delta_t.")
    ap.add_argument("--t-max", type=float, default=10.0, help="Max total time t to sweep (inclusive).")
    ap.add_argument("--shots", type=int, default=128, help="Shots per (instance, step-count).")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed.")

    ap.add_argument("--random-low", type=int, default=-5)
    ap.add_argument("--random-high", type=int, default=5)
    ap.add_argument("--random-density", type=float, default=1.0, help="Off-diagonal density for random Q.")
    ap.add_argument("--graph-p", type=float, default=0.3, help="Edge probability for MaxCut and MIS graphs.")
    ap.add_argument("--mis-lambda", type=float, default=2.0, help="Penalty lambda for MIS QUBO.")

    ap.add_argument("--trotter-order", type=int, default=2, choices=[1, 2])
    ap.add_argument("--scale-dynamics", type=str, default="maxabs", choices=["none", "maxabs"])

    ap.add_argument("--opt-ref", type=str, default="qa_best", choices=["qa_best", "exact"])
    ap.add_argument("--exact-max-n", type=int, default=16, help="Max n for opt_ref=exact (brute force).")

    ap.add_argument(
        "--aer-method",
        type=str,
        default="auto",
        choices=["auto", "statevector", "matrix_product_state"],
    )
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

    ap.add_argument(
        "--no-transpile-cache",
        action="store_true",
        help="Disable transpile-template reuse across instances.",
    )
    ap.add_argument(
        "--transpile-cache-mode",
        type=str,
        default="support",
        choices=["support", "full"],
        help="Cache template mode: 'support' keys by nonzero term pattern, 'full' keys by n/schedule only.",
    )
    ap.add_argument(
        "--no-cache-autodisable",
        action="store_true",
        help="Do not auto-disable cache when measured-template hit rate stays below threshold.",
    )
    ap.add_argument(
        "--cache-autodisable-min-attempts",
        type=int,
        default=12,
        help="Minimum measured-template requests before auto-disable check.",
    )
    ap.add_argument(
        "--cache-autodisable-min-hit-rate",
        type=float,
        default=0.05,
        help="Auto-disable cache when measured hit-rate falls below this threshold.",
    )
    ap.add_argument(
        "--estimator-diagnostics",
        action="store_true",
        help="Compute EstimatorV2 expectation-energy diagnostics alongside shot-based curves.",
    )
    ap.add_argument(
        "--estimator-precision",
        type=float,
        default=None,
        help="Optional precision target passed to EstimatorV2.run.",
    )

    ap.add_argument("--stats-method", type=str, default="mw", choices=["mw", "perm"])
    ap.add_argument("--perm-iterations", type=int, default=10000)
    ap.add_argument(
        "--scan-stop-p-holm-threshold",
        type=float,
        default=None,
        help=(
            "Optional --n-list early stop: require max Holm-adjusted p-value to be >= this threshold "
            "alongside the easy-case threshold criterion."
        ),
    )
    ap.add_argument(
        "--scan-stop-easy-case-threshold",
        type=float,
        default=None,
        help="Optional --n-list early stop threshold for overall mean easy-case rate.",
    )
    ap.add_argument(
        "--scan-stop-easy-case-op",
        type=str,
        default="le",
        choices=["le", "ge"],
        help="Comparison operator for easy-case threshold in scan-stop criterion: le means <=, ge means >=.",
    )
    ap.add_argument(
        "--scan-stop-min-n-evals",
        type=int,
        default=1,
        help="Minimum number of evaluated n values before applying scan-stop checks.",
    )

    ap.add_argument("--outdir", type=str, default="qa_out")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    if args.instances <= 0:
        raise SystemExit("--instances must be positive.")
    if args.shots <= 0:
        raise SystemExit("--shots must be positive.")
    if args.delta_t <= 0.0:
        raise SystemExit("--delta-t must be positive.")
    if args.t_max < 0.0:
        raise SystemExit("--t-max must be nonnegative.")
    if int(args.random_low) > int(args.random_high):
        raise SystemExit("--random-low must be <= --random-high.")
    if not (0.0 <= float(args.random_density) <= 1.0):
        raise SystemExit("--random-density must be in [0, 1].")
    if not (0.0 <= float(args.graph_p) <= 1.0):
        raise SystemExit("--graph-p must be in [0, 1].")
    if int(args.exact_max_n) <= 0:
        raise SystemExit("--exact-max-n must be positive.")
    if int(args.mps_auto_threshold) <= 0:
        raise SystemExit("--mps-auto-threshold must be positive.")
    if int(args.mps_max_bond_dimension) <= 0:
        raise SystemExit("--mps-max-bond-dimension must be positive.")
    if float(args.mps_truncation_threshold) < 0.0:
        raise SystemExit("--mps-truncation-threshold must be nonnegative.")
    if int(args.mps_omp_threads) <= 0:
        raise SystemExit("--mps-omp-threads must be positive.")
    if float(args.mis_lambda) <= 1.0:
        raise SystemExit(
            "--mis-lambda must be > 1.0 for the current MIS QUBO formulation "
            "(-sum x_i + lambda * sum x_i x_j)."
        )
    if args.stats_method == "perm" and int(args.perm_iterations) <= 0:
        raise SystemExit("--perm-iterations must be positive.")
    if args.estimator_precision is not None and float(args.estimator_precision) <= 0.0:
        raise SystemExit("--estimator-precision must be positive.")
    if int(args.cache_autodisable_min_attempts) <= 0:
        raise SystemExit("--cache-autodisable-min-attempts must be positive.")
    if not (0.0 <= float(args.cache_autodisable_min_hit_rate) <= 1.0):
        raise SystemExit("--cache-autodisable-min-hit-rate must be in [0, 1].")
    if int(args.scan_stop_min_n_evals) <= 0:
        raise SystemExit("--scan-stop-min-n-evals must be positive.")
    scan_stop_configured = (
        args.scan_stop_p_holm_threshold is not None or args.scan_stop_easy_case_threshold is not None
    )
    if scan_stop_configured and (
        args.scan_stop_p_holm_threshold is None or args.scan_stop_easy_case_threshold is None
    ):
        raise SystemExit(
            "--scan-stop-p-holm-threshold and --scan-stop-easy-case-threshold must be provided together."
        )
    if args.scan_stop_p_holm_threshold is not None and not (0.0 <= float(args.scan_stop_p_holm_threshold) <= 1.0):
        raise SystemExit("--scan-stop-p-holm-threshold must be in [0, 1].")
    if args.scan_stop_easy_case_threshold is not None and not (
        0.0 <= float(args.scan_stop_easy_case_threshold) <= 1.0
    ):
        raise SystemExit("--scan-stop-easy-case-threshold must be in [0, 1].")

    try:
        n_values = parse_n_list(args.n, args.n_list)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    scan_mode = bool(args.n_list)
    os.makedirs(args.outdir, exist_ok=True)

    if not scan_mode:
        run_benchmark_for_n(args, n=n_values[0], outdir=args.outdir, seed_base=int(args.seed))
        return

    scan_rows: List[Dict[str, Any]] = []
    scan_stop_enabled = bool(
        args.scan_stop_p_holm_threshold is not None and args.scan_stop_easy_case_threshold is not None
    )
    for idx, n in enumerate(n_values):
        n_outdir = os.path.join(args.outdir, f"n_{n}")
        stats = run_benchmark_for_n(args, n=n, outdir=n_outdir, seed_base=int(args.seed) + idx * 1_000_000)
        pvals = stats.get("pvals_one_sided_less") or {}
        pvals_holm = stats.get("pvals_holm_adjusted") or {}
        deltas = stats.get("cliffs_delta_smaller_is_better") or {}
        steps_summary = stats["steps_summary"]
        easy_case = stats.get("easy_case") or {}
        easy_case_rates = easy_case.get("rate") or {}
        p_holm_values = [float(v) for v in pvals_holm.values() if isinstance(v, (int, float))]
        p_holm_max = float(max(p_holm_values)) if p_holm_values else None
        easy_case_overall_mean_rate = easy_case.get("overall_mean_rate", "")
        easy_case_overall_mean_rate_num = (
            float(easy_case_overall_mean_rate)
            if isinstance(easy_case_overall_mean_rate, (int, float))
            else None
        )
        scan_stop_triggered = False
        scan_stop_reason = ""
        if scan_stop_enabled and (idx + 1) >= int(args.scan_stop_min_n_evals):
            p_stop = float(args.scan_stop_p_holm_threshold)
            easy_stop = float(args.scan_stop_easy_case_threshold)
            if p_holm_max is not None and easy_case_overall_mean_rate_num is not None:
                easy_condition = (
                    easy_case_overall_mean_rate_num <= easy_stop
                    if str(args.scan_stop_easy_case_op) == "le"
                    else easy_case_overall_mean_rate_num >= easy_stop
                )
                if p_holm_max >= p_stop and easy_condition:
                    scan_stop_triggered = True
                    scan_stop_reason = (
                        f"scan-stop triggered at n={n}: "
                        f"max_holm_p={p_holm_max:.4g} >= {p_stop:.4g} and "
                        f"overall_mean_easy_case_rate={easy_case_overall_mean_rate_num:.4g} "
                        f"{args.scan_stop_easy_case_op} {easy_stop:.4g}"
                    )

        scan_rows.append(
            {
                "n": int(n),
                "stats_method": str(args.stats_method),
                "random_median_steps": float(steps_summary["random"]["median"]),
                "maxcut_median_steps": float(steps_summary["maxcut"]["median"]),
                "mis_median_steps": float(steps_summary["mis"]["median"]),
                "easy_case_k_threshold": easy_case.get("k_threshold", ""),
                "random_easy_case_rate": easy_case_rates.get("random", ""),
                "maxcut_easy_case_rate": easy_case_rates.get("maxcut", ""),
                "mis_easy_case_rate": easy_case_rates.get("mis", ""),
                "overall_mean_easy_case_rate": easy_case.get("overall_mean_rate", ""),
                "scan_stop_triggered": int(scan_stop_triggered),
                "scan_stop_reason": scan_stop_reason,
                "scan_stop_p_holm_max": p_holm_max if p_holm_max is not None else "",
                "scan_stop_p_holm_threshold": args.scan_stop_p_holm_threshold if scan_stop_enabled else "",
                "scan_stop_easy_case_threshold": args.scan_stop_easy_case_threshold if scan_stop_enabled else "",
                "scan_stop_easy_case_op": args.scan_stop_easy_case_op if scan_stop_enabled else "",
                "p_random_vs_maxcut": pvals.get("random_vs_maxcut", ""),
                "p_random_vs_mis": pvals.get("random_vs_mis", ""),
                "p_holm_random_vs_maxcut": pvals_holm.get("random_vs_maxcut", ""),
                "p_holm_random_vs_mis": pvals_holm.get("random_vs_mis", ""),
                "delta_random_vs_maxcut": deltas.get("delta_random_vs_maxcut", ""),
                "delta_random_vs_mis": deltas.get("delta_random_vs_mis", ""),
                "n_outdir": n_outdir,
            }
        )
        if scan_stop_triggered:
            print(f"[scan] {scan_stop_reason}")
            break

    scan_csv_path = os.path.join(args.outdir, "scan_summary.csv")
    with open(scan_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "n",
                "stats_method",
                "random_median_steps",
                "maxcut_median_steps",
                "mis_median_steps",
                "easy_case_k_threshold",
                "random_easy_case_rate",
                "maxcut_easy_case_rate",
                "mis_easy_case_rate",
                "overall_mean_easy_case_rate",
                "scan_stop_triggered",
                "scan_stop_reason",
                "scan_stop_p_holm_max",
                "scan_stop_p_holm_threshold",
                "scan_stop_easy_case_threshold",
                "scan_stop_easy_case_op",
                "p_random_vs_maxcut",
                "p_random_vs_mis",
                "p_holm_random_vs_maxcut",
                "p_holm_random_vs_mis",
                "delta_random_vs_maxcut",
                "delta_random_vs_mis",
                "n_outdir",
            ],
        )
        w.writeheader()
        for row in scan_rows:
            w.writerow(row)
    print(f"\nWrote: {scan_csv_path}")


if __name__ == "__main__":
    main()
