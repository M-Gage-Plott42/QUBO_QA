import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

import qa_adiabatic_steps_bench as bench


class TestEnergyConsistency(unittest.TestCase):
    def test_ising_energy_matches_dimod_bqm_energy_for_sampled_bits(self) -> None:
        rng = np.random.default_rng(12345)
        n = 6

        makers = [
            lambda: bench.make_random_qubo_bqm(rng, n=n, low=-5, high=5, density=1.0)[0],
            lambda: bench.make_maxcut_qubo_bqm(rng, n=n, p=0.3)[0],
            lambda: bench.make_mis_qubo_bqm(rng, n=n, p=0.3, penalty_lambda=2.0)[0],
        ]

        for make_bqm in makers:
            for _ in range(4):
                bqm = make_bqm()
                model_eval = bench.bqm_to_ising_arrays(bqm, n=n)
                for _ in range(16):
                    bits = rng.integers(0, 2, size=n)
                    bitstring_lsb0 = "".join(str(int(b)) for b in bits)
                    e_ising = bench.ising_energy_from_bitstring(bitstring_lsb0, model_eval)
                    sample = {i: int(bits[i]) for i in range(n)}
                    e_bqm = float(bqm.energy(sample))
                    self.assertAlmostEqual(e_ising, e_bqm, places=9)


class TestCliGuards(unittest.TestCase):
    def test_mis_lambda_guard_rejects_invalid_value(self) -> None:
        proc = subprocess.run(
            [
                sys.executable,
                "qa_adiabatic_steps_bench.py",
                "-n",
                "4",
                "--instances",
                "1",
                "--t-max",
                "0.5",
                "--shots",
                "8",
                "--mis-lambda",
                "1.0",
                "--no-plots",
                "--outdir",
                "/tmp/qa_test_invalid_mis_lambda",
            ],
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("--mis-lambda must be greater than max(1.0, --mis-node-weight-high)", proc.stderr)


class TestSuccessSemantics(unittest.TestCase):
    def test_summary_has_cumulative_and_instantaneous_success_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "qa_semantics"
            proc = subprocess.run(
                [
                    sys.executable,
                    "qa_adiabatic_steps_bench.py",
                    "-n",
                    "4",
                    "--instances",
                    "3",
                    "--t-max",
                    "1.0",
                    "--shots",
                    "16",
                    "--aer-method",
                    "statevector",
                    "--opt-ref",
                    "exact",
                    "--no-plots",
                    "--outdir",
                    str(outdir),
                ],
                cwd=Path(__file__).resolve().parents[1],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stdout + "\n" + proc.stderr)

            summary = json.loads((outdir / "summary.json").read_text(encoding="utf-8"))
            cumulative = summary["success_probability_cumulative_at_tmax"]
            instantaneous = summary["success_probability_instantaneous_at_tmax"]
            legacy = summary["success_probability_at_tmax"]

            self.assertEqual(legacy, cumulative)
            for fam in ("random", "maxcut", "mis"):
                self.assertIn(fam, cumulative)
                self.assertIn(fam, instantaneous)
                self.assertGreaterEqual(float(cumulative[fam]), float(instantaneous[fam]))


if __name__ == "__main__":
    unittest.main()
