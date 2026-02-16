#!/usr/bin/env python3
"""
Shared utilities for QA/SA/PRR comparison workflows.

This module centralizes:
- deterministic instance-bank generation across families,
- instance manifest writing with stable digests,
- common schema/version constants used by comparison artifacts.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from typing import Any, Dict, List, Tuple

import dimod
import numpy as np

from qa_adiabatic_steps_bench import (
    make_maxcut_qubo_bqm,
    make_mis_qubo_bqm,
    make_random_qubo_bqm,
)


COMPARISON_SCHEMA_VERSION = "qa_sa_prr_compare_v1"
FAMILY_ORDER: Tuple[str, str, str] = ("random", "maxcut", "mis")


@dataclass
class InstanceEntry:
    family: str
    instance_id: int
    seed: int
    bqm: dimod.BinaryQuadraticModel
    meta: Dict[str, Any]
    bqm_digest: str


def _canonical_bqm_lines(bqm: dimod.BinaryQuadraticModel) -> List[str]:
    lines: List[str] = []

    for var, bias in sorted(((int(k), float(v)) for k, v in bqm.linear.items()), key=lambda kv: kv[0]):
        lines.append(f"L,{var},{bias:.17g}")

    quad_items: List[Tuple[int, int, float]] = []
    for (u, v), bias in bqm.quadratic.items():
        i = int(u)
        j = int(v)
        if i > j:
            i, j = j, i
        quad_items.append((i, j, float(bias)))
    for i, j, bias in sorted(quad_items):
        lines.append(f"Q,{i},{j},{bias:.17g}")

    lines.append(f"O,{float(bqm.offset):.17g}")
    return lines


def bqm_digest_sha256(bqm: dimod.BinaryQuadraticModel) -> str:
    payload = "\n".join(_canonical_bqm_lines(bqm)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def summarize_instance_meta(family: str, meta: Dict[str, Any]) -> Dict[str, Any]:
    if family == "random":
        qsym = np.asarray(meta.get("Qsym", []), dtype=float)
        if qsym.size == 0:
            return {"density": meta.get("density")}
        return {
            "density": float(meta.get("density", 1.0)),
            "shape": [int(v) for v in qsym.shape],
            "min": float(np.min(qsym)),
            "max": float(np.max(qsym)),
            "mean": float(np.mean(qsym)),
            "std": float(np.std(qsym)),
            "trace": float(np.trace(qsym)),
        }

    out: Dict[str, Any] = {}
    for key in (
        "p",
        "m",
        "lambda",
        "node_weight_low",
        "node_weight_high",
        "edge_span_before",
        "edge_span_after",
    ):
        if key in meta:
            out[key] = meta[key]
    return out


def build_instance_bank(
    *,
    n: int,
    instances_per_family: int,
    seed_base: int,
    random_low: int,
    random_high: int,
    random_density: float,
    graph_p: float,
    maxcut_weight_low: int,
    maxcut_weight_high: int,
    mis_lambda: float,
    mis_node_weight_low: int,
    mis_node_weight_high: int,
) -> Dict[str, List[InstanceEntry]]:
    bank: Dict[str, List[InstanceEntry]] = {family: [] for family in FAMILY_ORDER}

    for fam_idx, family in enumerate(FAMILY_ORDER):
        for inst_id in range(int(instances_per_family)):
            entry_seed = int(seed_base) + fam_idx * 100_000 + inst_id
            rng = np.random.default_rng(entry_seed)

            if family == "random":
                bqm, meta = make_random_qubo_bqm(
                    rng=rng,
                    n=int(n),
                    low=int(random_low),
                    high=int(random_high),
                    density=float(random_density),
                )
            elif family == "maxcut":
                bqm, meta = make_maxcut_qubo_bqm(
                    rng=rng,
                    n=int(n),
                    p=float(graph_p),
                    weight_low=int(maxcut_weight_low),
                    weight_high=int(maxcut_weight_high),
                )
            elif family == "mis":
                bqm, meta = make_mis_qubo_bqm(
                    rng=rng,
                    n=int(n),
                    p=float(graph_p),
                    penalty_lambda=float(mis_lambda),
                    node_weight_low=int(mis_node_weight_low),
                    node_weight_high=int(mis_node_weight_high),
                )
            else:
                raise RuntimeError(f"Unknown family: {family}")

            bank[family].append(
                InstanceEntry(
                    family=family,
                    instance_id=int(inst_id),
                    seed=int(entry_seed),
                    bqm=bqm,
                    meta=meta,
                    bqm_digest=bqm_digest_sha256(bqm),
                )
            )

    return bank


def write_instance_bank_manifest(
    *,
    outdir: str,
    n: int,
    instances_per_family: int,
    seed_base: int,
    generation_config: Dict[str, Any],
    bank: Dict[str, List[InstanceEntry]],
) -> str:
    os.makedirs(outdir, exist_ok=True)

    payload: Dict[str, Any] = {
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "n": int(n),
        "instances_per_family": int(instances_per_family),
        "seed_base": int(seed_base),
        "generation_config": generation_config,
        "families": {},
    }

    for family in FAMILY_ORDER:
        payload["families"][family] = []
        for entry in bank[family]:
            payload["families"][family].append(
                {
                    "instance_id": int(entry.instance_id),
                    "seed": int(entry.seed),
                    "bqm_digest": str(entry.bqm_digest),
                    "meta": summarize_instance_meta(family, entry.meta),
                }
            )

    manifest_path = os.path.join(outdir, "instance_bank.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return manifest_path
