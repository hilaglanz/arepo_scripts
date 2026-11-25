"""Synthetic physical checks; no AREPO dependencies or snapshots are loaded."""
import copy
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyze_bowen_launch as diagnostic


LAW = {"temperature": 1100.0, "width": 50.0, "kappa_max": 8.0, "rho_min": 1.0}


def gas_fixture():
    """Unsorted gas: one interior cell, three selected, one exterior."""
    center = np.array([10.0, 20.0, 30.0])
    reference_velocity = np.array([100.0, -200.0, 300.0])
    displacement = np.array([[3, 0, 0], [.5, 0, 0], [0, 2, 0],
                             [0, 0, 4], [6, 0, 0]], dtype=float)
    velocity = np.array([[-2, 8, 9], [0, 0, 0], [10, 5, 9],
                         [8, -7, 0], [99, 0, 0]], dtype=float)
    data = {
        "pos": center + displacement,
        "vel": reference_velocity + velocity,
        "mass": np.array([3, 11, 2, 4, 13], dtype=float),
        "rho": np.full(5, 2.0),
        "temp": np.full(5, LAW["temperature"]),
        "ka_r": np.full(5, 0.5),
        "ka_s": np.full(5, 1.0),
        # Large tangential flux components must not enter radial acceleration.
        "frad": np.array([[-20, 800, 900], [0, 0, 0], [1000, 50, 900],
                          [800, -700, 0], [99, 0, 0]], dtype=float),
    }
    return data, center, reference_velocity


def force_fixture(masses, velocities, densities=None, temperatures=None):
    """Use a flux giving Gamma = opacity for straightforward physical cases."""
    mass = np.asarray(masses, dtype=float)
    count = len(mass)
    radius = np.arange(2, count + 2, dtype=float)
    enclosed_mass = 10 + np.cumsum(mass)
    gravity = diagnostic.GRAVITY * enclosed_mass / radius**2
    return {
        "radius": radius,
        "mass": mass,
        "enclosed_mass": enclosed_mass,
        "vr": np.asarray(velocities, dtype=float),
        "rho": np.full(count, 2.0) if densities is None else np.asarray(densities, dtype=float),
        "temp": np.full(count, 1100.0) if temperatures is None else np.asarray(temperatures, dtype=float),
        "flux_radial": diagnostic.LIGHT_SPEED * gravity,
        "ka_r": np.full(count, 0.5),
        "ka_s": np.full(count, 1.0),
    }


class RegionTests(unittest.TestCase):
    def test_common_boost_and_translation_leave_shell_kinematics_unchanged(self):
        original, center, reference = gas_fixture()
        transformed = copy.deepcopy(original)
        shift = np.array([83.0, -11.0, 27.0])
        boost = np.array([44.0, 52.0, -18.0])
        transformed["pos"] += shift
        transformed["vel"] += boost
        before = diagnostic.take_region(original, center, 7, reference, 1, 5, block_size=2)
        after = diagnostic.take_region(transformed, center + shift, 7,
                                       reference + boost, 1, 5, block_size=2)
        for name in before:
            np.testing.assert_allclose(before[name], after[name], err_msg=name)
        np.testing.assert_allclose(before["vr"], [5, -2, 0])

    def test_gravity_includes_unplotted_interior_gas_and_sorts_radii(self):
        data, center, reference = gas_fixture()
        cells = diagnostic.take_region(data, center, 7, reference, 1, 5)
        np.testing.assert_array_equal(cells["radius"], [2, 3, 4])
        np.testing.assert_array_equal(cells["mass"], [2, 3, 4])
        # Core 7 + interior gas 11, then each successive selected mass.
        np.testing.assert_array_equal(cells["enclosed_mass"], [20, 23, 27])

    def test_vector_flux_preserves_sign_and_ignores_tangential_components(self):
        data, center, reference = gas_fixture()
        cells = diagnostic.take_region(data, center, 7, reference, 1, 5)
        np.testing.assert_allclose(cells["flux_radial"], [50, -20, 0])
        diagnostic.add_force_diagnostics(cells, LAW)
        self.assertGreater(cells["gamma_total"][0], 0)
        self.assertLess(cells["gamma_total"][1], 0)
        self.assertEqual(cells["gamma_total"][2], 0)

    def test_missing_required_fields_raise_instead_of_inventing_values(self):
        for missing in ("pos", "mass", "vel", "rho", "temp", "ka_r", "frad", "ka_s"):
            with self.subTest(field=missing):
                data, center, reference = gas_fixture()
                del data[missing]
                data["fradr"] = np.ones(5)  # Scalar flux cannot replace required vector silently.
                with self.assertRaisesRegex(ValueError, missing):
                    diagnostic.take_region(data, center, 7, reference, 1, 5)

    def test_explicit_non_bowen_missing_scatter_permission_only_fills_scatter(self):
        data, center, reference = gas_fixture()
        del data["ka_s"]
        cells = diagnostic.take_region(data, center, 7, reference, 1, 5,
                                       allow_zero_scatter=True)
        np.testing.assert_array_equal(cells["ka_s"], np.zeros(3))
        diagnostic.add_force_diagnostics(cells, LAW)
        np.testing.assert_array_equal(cells["gamma_added"], np.zeros(3))
        np.testing.assert_array_equal(cells["gamma_total"], cells["gamma_base"])
        data, center, reference = gas_fixture()
        del data["ka_s"], data["ka_r"]
        with self.assertRaisesRegex(ValueError, "ka_r"):
            diagnostic.take_region(data, center, 7, reference, 1, 5,
                                   allow_zero_scatter=True)

    def test_empty_selected_region_is_valid_and_has_no_force_measurement(self):
        data, center, reference = gas_fixture()
        cells = diagnostic.take_region(data, center, 7, reference, 7, 8)
        diagnostic.add_force_diagnostics(cells, LAW)
        row = diagnostic.summarize(cells, np.array([7, 8]), LAW)[0]
        self.assertEqual(row["gas_mass_Msun"], 0)
        self.assertEqual(row["net_Msun_per_yr"], 0)
        self.assertTrue(np.isnan(row["gamma_total_mass_mean"]))
        self.assertTrue(np.isnan(row["fraction_of_inward_mass_gamma_total_gt1"]))


class OpacityAndSummaryTests(unittest.TestCase):
    def test_strict_density_gate_and_half_max_at_condensation_temperature(self):
        cells = force_fixture([1, 2, 3, 4], [-1, 2, -3, 4],
                              densities=[1, 1.1, 2, 2], temperatures=[900, 1100, 1150, 1050])
        diagnostic.add_force_diagnostics(cells, LAW)
        expected = np.array([0, 4, 8 / (1 + np.exp(1)), 8 / (1 + np.exp(-1))])
        np.testing.assert_allclose(cells["kappa_trial"], expected)
        cells["ka_s"] = expected.copy()
        diagnostic.add_force_diagnostics(cells, LAW)
        row = diagnostic.summarize(cells, np.array([1, 6]), LAW)[0]
        self.assertAlmostEqual(row["density_eligible_mass_fraction"], .9)
        self.assertAlmostEqual(row["cool_and_dense_mass_fraction"], .6)
        self.assertAlmostEqual(row["halfmax_added_opacity_mass_fraction"], .6)

    def test_trial_replaces_existing_bowen_instead_of_counting_it_twice(self):
        cells = force_fixture([1], [-1])
        cells["ka_r"][:] = 2
        cells["ka_s"][:] = 20
        diagnostic.add_force_diagnostics(cells, LAW)
        np.testing.assert_allclose(cells["gamma_total"], [22])
        np.testing.assert_allclose(cells["gamma_trial_total"], [6])
        # A matching stored law makes the counterfactual identical to actual.
        cells["ka_s"][:] = cells["kappa_trial"]
        diagnostic.add_force_diagnostics(cells, LAW)
        np.testing.assert_allclose(cells["gamma_trial_total"], cells["gamma_total"])

    def test_outward_inward_net_flux_and_mass_weighted_velocity(self):
        data, center, reference = gas_fixture()
        cells = diagnostic.take_region(data, center, 7, reference, 1, 5)
        diagnostic.add_force_diagnostics(cells, LAW)
        row = diagnostic.summarize(cells, np.array([1, 5]), LAW)[0]
        conversion = diagnostic.YEAR / diagnostic.MSUN / 4
        self.assertAlmostEqual(row["outward_Msun_per_yr"] / conversion, 10)
        self.assertAlmostEqual(row["inward_Msun_per_yr"] / conversion, 6)
        self.assertAlmostEqual(row["net_Msun_per_yr"] / conversion, 4)
        self.assertAlmostEqual(row["radial_velocity_mass_mean_km_s"] * 1e5, 4 / 9)
        self.assertAlmostEqual(row["inward_mass_fraction"], 3 / 9)
        self.assertAlmostEqual(row["density_volume_mean_g_cm3"], 2)

    def test_infalling_gamma_fraction_uses_only_infalling_mass_denominator(self):
        cells = force_fixture([1, 3, 6, 2], [-1, -1, 2, -1])
        diagnostic.add_force_diagnostics(cells, LAW)
        # Explicit exact values also check strict >1, rather than >=1.
        cells["gamma_total"] = np.array([2.0, .5, 1.2, 1.0])
        cells["gamma_trial_total"] = np.array([2.0, 2.0, 1.2, 1.0])
        row = diagnostic.summarize(cells, np.array([1, 6]), LAW)[0]
        self.assertAlmostEqual(row["gamma_total_gt1_mass_fraction"], 7 / 12)
        self.assertAlmostEqual(row["fraction_of_inward_mass_gamma_total_gt1"], 1 / 6)
        self.assertAlmostEqual(row["fraction_of_inward_mass_gamma_trial_total_gt1"], 4 / 6)
        cells["vr"][:] = 2
        no_infall = diagnostic.summarize(cells, np.array([1, 6]), LAW)[0]
        self.assertTrue(np.isnan(no_infall["fraction_of_inward_mass_gamma_total_gt1"]))

    def test_adjacent_shells_count_boundary_cells_once_and_omit_upper_edge(self):
        cells = force_fixture([1, 2, 3], [1, 1, 1])  # radii 2, 3, 4
        diagnostic.add_force_diagnostics(cells, LAW)
        rows = diagnostic.summarize(cells, np.array([1, 3, 4]), LAW)
        np.testing.assert_array_equal([row["cell_count"] for row in rows], [1, 1])
        np.testing.assert_allclose([row["gas_mass_Msun"] * diagnostic.MSUN for row in rows], [1, 2])


if __name__ == "__main__":
    unittest.main()
