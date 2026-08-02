"""Run with: sage -python -m unittest discover -v"""

import json
from random import Random
import unittest
from unittest.mock import patch

from core import Slipway
from reference import Poseidon
from search import construct, is_mds, search
from verify_paper_rf8_rp20 import INSTANCE, verify_paper
from sage.all import GF, matrix, vector


class SlipwayTests(unittest.TestCase):
    def setUp(self):
        self.instance = Slipway(json.loads(INSTANCE.read_text()))

    def test_upstream_verifier_uses_our_matrix_and_constants(self):
        upstream = self.instance.reference
        self.assertIsInstance(upstream, Poseidon)
        self.assertEqual(upstream.mds, [list(row) for row in self.instance.matrix.rows()])
        self.assertEqual(upstream.round_constants, [list(row) for row in self.instance.constants])
        # The normal Poseidon permutation matches the paper's round schedule.
        inputs = [0] * self.instance.t
        with patch.object(upstream, "permutation", wraps=upstream.permutation) as evaluate:
            actual = self.instance.reference_permutation(inputs)
            evaluate.assert_called_once_with(inputs)
        state = vector(self.instance.field, inputs)
        for index in range(self.instance.rf + self.instance.rp):
            state = self.instance.round(state, index)
        self.assertEqual(actual, list(state))

    def test_fixed_control_root_extraction(self):
        self.instance.certify()
        result = self.instance.solve(self.instance.control)
        self.assertEqual(result["inactive_partial_rounds"], 14)
        self.assertEqual(result["output_degrees"], [59049, 59049])
        self.assertEqual(result["candidates"], [
            {"root": 18129142, "outputs": [0, 2007855200]},
            {"root": 1404166674, "outputs": [0, 0]},
        ])
        self.assertEqual(len(result["solutions"]), 1)
        solution = result["solutions"][0]
        self.assertEqual(solution["input"], [
            0, 0, 1994407524, 1518406169, 1869697938, 1680758970,
            1974616136, 2130375035, 5166025, 909317996, 1154992185,
            1870851706, 82420390, 1232840978, 1169612330, 1690942889,
        ])
        self.assertEqual(solution["output"], [
            0, 0, 1080431931, 1291539939, 842366831, 410517411,
            390312851, 1938087994, 1744025294, 601741681, 1822652989,
            1131652874, 158028812, 266015611, 365198253, 1770326166,
        ])

    def test_a_root_of_only_one_output_is_not_a_solution(self):
        # At w=1, a root of output 0 fails the second output constraint.
        result = self.instance.solve(1)
        self.assertIn({"root": 1108395457, "outputs": [0, 1772989505]},
                      result["candidates"])
        self.assertEqual(result["solutions"], [])

    def test_collapsed_control_has_one_more_inactive_sbox(self):
        f0, f1, degrees = self.instance.output_polynomials(244601253)
        self.assertEqual(degrees[:15], [0] * 15)
        self.assertEqual(degrees[15], 1)
        self.assertEqual([f0.degree(), f1.degree()], [19683, 19683])

    def test_changed_prefix_constant_fails_certificate(self):
        self.instance.constants[1][0] += 1
        with self.assertRaisesRegex(ValueError, "Prefix identity failed"):
            self.instance.certify()

    def test_out_of_field_controls_are_rejected(self):
        for control in (-1, self.instance.p):
            with self.subTest(control=control):
                with self.assertRaisesRegex(ValueError, "canonical field element"):
                    self.instance.solve(control)

    def test_inverse_prefix_recovers_paper_controls(self):
        for w in (0, 1, 244601253, self.instance.control, self.instance.p - 1):
            g, h = self.instance.field(1), self.instance.field(w)
            for lg, lh in self.instance.ratios:
                g, h = (g + lg*h)**3, (g + lh*h)**3
            for X in (1, self.instance.p - 1):
                self.assertEqual(self.instance.recover_control(X*g, X*(h-g)), (X, w))
        self.assertEqual(self.instance.recover_control(0, 0), (0, 0))


class PaperVerifierTests(unittest.TestCase):
    def test_supplied_witness_is_verified_without_search(self):
        data = json.loads(INSTANCE.read_text())
        with patch.object(Slipway, "solve", side_effect=AssertionError("No root search")), \
             patch.object(Slipway, "solve_joint", side_effect=AssertionError("No joint search")):
            result = verify_paper(data)
        self.assertEqual(result["root"], data["witness"]["root"])

    def test_changed_witness_is_rejected(self):
        for field in ("root", "output"):
            with self.subTest(field=field):
                data = json.loads(INSTANCE.read_text())
                if field == "root":
                    data["witness"][field] += 1
                else:
                    data["witness"][field][2] += 1
                with self.assertRaisesRegex(ValueError, "Supplied"):
                    verify_paper(data)


class SearchTests(unittest.TestCase):
    def test_fresh_matrix_and_control_search(self):
        reports = [search(seed, p=257) for seed in (0, 1)]
        self.assertNotEqual(reports[0]["instance"]["matrix"],
                            reports[1]["instance"]["matrix"])
        self.assertEqual(reports[0]["instance"]["round_constants"],
                         reports[1]["instance"]["round_constants"])
        for report in reports:
            instance = Slipway(report["instance"])
            instance.certify()
            self.assertTrue(is_mds(instance.matrix))
            result = report["result"]
            self.assertEqual(result["output_degrees"], [27, 27])
            self.assertEqual(result["inactive_partial_rounds"], 4)
            self.assertTrue(result["solutions"])
            for solution in result["solutions"]:
                self.assertEqual(solution["input"][:2], [0, 0])
                self.assertEqual(solution["output"][:2], [0, 0])
                self.assertEqual(instance.reference_permutation(solution["input"]),
                                 solution["output"])
            # Exhaustively check the polynomial representation on this small field.
            f0, f1, _ = instance.output_polynomials(result["control"])
            for root in range(instance.p):
                inputs = instance.recover_input(root, result["control"])
                self.assertEqual([int(f0(root)), int(f1(root))],
                                 instance.reference_permutation(inputs)[:2])

    def test_configurable_koalabear_rounds_and_upstream_constants(self):
        for rf, rp in ((4, 4), (4, 5), (6, 7)):
            with self.subTest(rf=rf, rp=rp):
                report = search(rf=rf, rp=rp)
                data, result = report["instance"], report["result"]
                t, p = data["width"], data["prime"]
                self.assertEqual((data["full_rounds"], data["partial_rounds"]), (rf, rp))
                self.assertEqual(data["round_constants"], Poseidon(p, 3, t, rf, rp).round_constants)
                expected = 3**(rf//2 + max(0, rp-(t-2)))
                self.assertEqual(result["output_degrees"], [expected, expected])
                instance = Slipway(data)
                for solution in result["solutions"]:
                    self.assertEqual(solution["input"][:2], [0, 0])
                    self.assertEqual(solution["output"][:2], [0, 0])
                    self.assertEqual(instance.reference_permutation(solution["input"]),
                                     solution["output"])

    def test_unsupported_parameters_fail_before_search(self):
        for params in ({"rf": 3}, {"rp": -1}, {"rf": 6, "t": 6},
                       {"p": 259}, {"rp": 7}, {"max_attempts": 0}):
            with self.subTest(params=params):
                with self.assertRaises(ValueError):
                    search(**params)

    def test_four_round_prefix_construction(self):
        p, t, rf, rp = 2130706433, 12, 8, 10
        constants = Poseidon(p, 3, t, rf, rp).round_constants
        rng, candidate = Random(0), None
        # Test the longer algebraic construction without millions of minor
        # checks in the routine suite; MDS is tested separately at small widths.
        with patch("search.is_mds", return_value=True):
            for _ in range(100):
                candidate = construct(rng, p, t, rf, rp, constants)
                if candidate is not None:
                    break
        self.assertIsNotNone(candidate)
        _, instance = candidate  # construct() already proved the prefix identity.
        f0, f1, degrees = instance.output_polynomials(1)
        self.assertEqual(degrees, [0]*rp)
        self.assertEqual((f0.degree(), f1.degree()), (81, 81))

    def test_mds_requires_smaller_minors_as_well_as_invertibility(self):
        # Invertible, all entries nonzero, but its top-left 2x2 minor vanishes.
        M = matrix(GF(257), [[1, 1, 1], [1, 1, 2], [1, 2, 1]])
        self.assertTrue(M.is_invertible())
        self.assertFalse(is_mds(M))


if __name__ == "__main__":
    unittest.main()
