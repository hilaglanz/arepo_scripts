"""Small synthetic corotation checks; no snapshot imports, IC writes, or simulations.

Run with: python -B tests/test_giant_binary_corotation.py
The orbital checks use the sibling arepo-snap-util source when it is available.
"""
import argparse
import ast
import contextlib
import io
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "GiantBinaries.py"
LIBRARY = ROOT.parent / "arepo-snap-util/stellar_ics/multiple_star.py"
G, RSUN, MSUN = 6.67430e-8, 6.957e10, 1.98847e33


def execute_nodes(nodes, namespace, filename):
    tree = ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[]))
    exec(compile(tree, str(filename), "exec"), namespace)


class Component:
    """In-memory data stand-in, including an off-center core as the last row."""

    def __init__(self, data, radius=1e13 / RSUN, **kwargs):
        self.data = {key: value.copy() for key, value in data.items()}
        self.mass = self.data["mass"].sum()
        self.radius = radius
        self.offset = np.zeros(3)

    def get_radius(self):
        return self.radius * RSUN

    def get_position(self):
        return self.data["pos"]

    def get_shifted_position(self):
        return self.data["pos"] + self.offset

    def get_velocity(self):
        return self.data["vel"]

    def set_velocity(self, velocity):
        self.data["vel"] = velocity


def donor_data():
    mass = np.array([2., 3., 4., 5.]) * (2e33 / 14.)
    position = np.array([[1., 2., .3], [-2., 1., .8],
                         [.2, -3., -1.], [.1, .4, .2]]) * 1e12
    position -= np.average(position, axis=0, weights=mass)
    velocity = np.array([[20., 30., 9.], [-60., 55., -8.],
                         [5., -8., 2.], [4., -5., 6.]]) * 1e3
    velocity -= np.average(velocity, axis=0, weights=mass)
    return dict(pos=position, vel=velocity, mass=mass, type=np.array([0, 0, 0, 1]),
                rho=np.ones(3) * 1e-8, u=np.ones(3) * 1e12,
                xnuc=np.tile([.7, .28, .02], (3, 1)))


def point_mass(mass):
    return Component(dict(pos=np.zeros((1, 3)), vel=np.zeros((1, 3)),
                          mass=np.array([mass]), type=np.array([1])))


class BuilderForwardingTests(unittest.TestCase):
    def setUp(self):
        self.binary = Mock()
        self.binary.orbital_separation = 1e14
        self.reader = Mock(return_value=SimpleNamespace(data=donor_data(), boxsize=1e15))
        self.namespace = dict(np=np, argparse=argparse, rsol=RSUN, msol=MSUN, G=G,
                              gadget_readsnapname=self.reader,
                              SnapshotComponent=Component, PointMassComponent=point_mass,
                              filter_donor_background=lambda *args: None,
                              enclosing_boxsize=lambda *args: 1e15,
                              roche_distance=lambda q: .4,
                              BinarySystemWithLogIDs=lambda **kwargs: self.binary)
        self.tree = ast.parse(BUILDER.read_text())
        execute_nodes([node for node in self.tree.body
                       if isinstance(node, ast.FunctionDef)],
                      self.namespace, BUILDER)
        self.spin = Mock()
        self.namespace["apply_donor_rotation"] = self.spin
        self.namespace["filter_donor_background"] = lambda *args: None
        self.namespace["enclosing_boxsize"] = lambda *args: 1e15

    def assert_forwarded(self, expected):
        call = self.binary.add_components_as_binary.call_args
        self.assertEqual(call.kwargs["corotation_factor"], 0.)  # Library adds only orbital motion.
        self.assertIs(call.kwargs["corotating_at_rlof"], False)
        self.assertEqual(call.kwargs["e"], 0.)
        primary, secondary = call.args
        omega = expected * np.sqrt(G * (primary.mass + secondary.mass) / 1e14**3)
        if expected == 0.:
            self.spin.assert_not_called()
        else:
            self.spin.assert_called_once_with(primary, omega, 250. * RSUN, 500. * RSUN,
                                              profile="wind", spin_mode="synchronize")
        self.binary.create_ics.assert_called_once()

    def test_builder_default_and_explicit_factor(self):
        for options, expected in (({}, 0.), ({"corotation_factor": 1.,
                                             "corotation_outer_radius_rsol": 500.}, 1.)):
            with self.subTest(expected=expected), contextlib.redirect_stdout(io.StringIO()):
                self.binary.reset_mock()
                self.spin.reset_mock()
                self.namespace["AddPointMassToFile"](
                    "unused_snapshot", "unused_ic", 1e33, None,
                    giant_radius_rsol=250., **options)
                self.assert_forwarded(expected)

    def test_cli_default_and_explicit_factor_reach_builder(self):
        for options, expected in (([], 0.), (["--corotation_factor", "1",
                                            "--corotation_outer_radius", "500"], 1.)):
            with self.subTest(expected=expected), contextlib.redirect_stdout(io.StringIO()):
                self.binary.reset_mock()
                self.spin.reset_mock()
                self.namespace.update(__name__="__main__", sys=sys)
                with patch.object(sys, "argv", [str(BUILDER), "--giant_radius", "250"] + options):
                    execute_nodes([node for node in self.tree.body if isinstance(node, ast.If)],
                                  self.namespace, BUILDER)
                self.assert_forwarded(expected)

    def test_nonfinite_factor_fails_before_loading_snapshot(self):
        for factor in (np.nan, np.inf, -np.inf):
            with self.subTest(factor=factor), self.assertRaisesRegex(ValueError, "finite"):
                self.namespace["AddPointMassToFile"](
                    "unused_snapshot", "unused_ic", 1e33, None, corotation_factor=factor)
        self.reader.assert_not_called()

    def test_missing_outer_radius_fails_before_loading_snapshot(self):
        with self.assertRaisesRegex(ValueError, "outer_radius"):
            self.namespace["AddPointMassToFile"](
                "unused_snapshot", "unused_ic", 1e33, None, corotation_factor=1.)
        self.reader.assert_not_called()

    def test_cli_density_cut_reaches_boundary_finder(self):
        finder = Mock(return_value=500. * RSUN)
        self.namespace.update(__name__="__main__", sys=sys, corotation_radius_from_density=finder)
        with patch.object(sys, "argv", [str(BUILDER), "--giant_radius", "250",
                                        "--corotation_factor", "1", "--corotation_density_cut", "1e-10"]), \
                contextlib.redirect_stdout(io.StringIO()):
            execute_nodes([node for node in self.tree.body if isinstance(node, ast.If)],
                          self.namespace, BUILDER)
        self.assert_forwarded(1.)
        primary = self.binary.add_components_as_binary.call_args.args[0]
        finder.assert_called_once_with(primary, 250. * RSUN, 1e-10)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            self.namespace["InitParser"]().parse_args(
                ["--corotation_outer_radius", "500", "--corotation_density_cut", "1e-10"])

    def test_invalid_density_cut_fails_before_loading_snapshot(self):
        for cut in (0., -1., np.nan, np.inf):
            with self.subTest(cut=cut), self.assertRaisesRegex(ValueError, "density_cut"):
                self.namespace["AddPointMassToFile"](
                    "unused_snapshot", "unused_ic", 1e33, None,
                    corotation_factor=1., corotation_density_cut=cut)
        self.reader.assert_not_called()

    def test_cli_explicit_legacy_profile_and_add_mode(self):
        self.namespace.update(__name__="__main__", sys=sys)
        arguments = [str(BUILDER), "--giant_radius", "250", "--corotation_factor", ".4",
                     "--corotation_outer_radius", "500", "--corotation_profile", "taper",
                     "--corotation_spin_mode", "add"]
        with patch.object(sys, "argv", arguments), contextlib.redirect_stdout(io.StringIO()):
            execute_nodes([node for node in self.tree.body if isinstance(node, ast.If)],
                          self.namespace, BUILDER)
        self.assertEqual(self.spin.call_args.kwargs, {"profile": "taper", "spin_mode": "add"})

    def test_density_crossing_ignores_later_dense_shell(self):
        edges = np.geomspace(1., 32., 129)
        radius = np.r_[np.sqrt(edges[:-1] * edges[1:]), 32.]
        density = 10. / radius**2
        density[radius > 8.] = 100.  # Detached, dense outer shell must not set boundary.
        position = np.column_stack((radius, np.zeros_like(radius), np.zeros_like(radius)))
        data = dict(pos=np.r_[position, -position, [[0., 0., 0.]]],
                    mass=np.r_[density, density, 1e6],
                    rho=np.r_[density, density],
                    vel=np.zeros((2 * len(radius) + 1, 3)),
                    type=np.r_[np.zeros(2 * len(radius)), 1])
        primary = Component(data)
        with contextlib.redirect_stdout(io.StringIO()):
            boundary = self.namespace["corotation_radius_from_density"](primary, 1., 1.)
        self.assertAlmostEqual(boundary, np.sqrt(10.), places=12)
        for key in data:
            np.testing.assert_array_equal(primary.data[key], data[key])

    def test_density_crossing_uses_volume_average(self):
        edges = np.geomspace(1., 16., 129)
        radius = np.r_[np.sqrt(edges[:-1] * edges[1:]), 16.]
        density = np.where(radius < 4., 4., .1)
        # Dense, tiny-volume clumps should not outweigh the diffuse shell volume.
        position = np.column_stack((np.r_[radius, radius], np.zeros(2 * len(radius)),
                                    np.zeros(2 * len(radius))))
        data = dict(pos=position, rho=np.r_[density, np.full(len(radius), 100.)],
                    mass=np.r_[100. * density, np.full(len(radius), .1)],
                    vel=np.zeros_like(position), type=np.zeros(len(position)))
        primary = Component(data)
        with contextlib.redirect_stdout(io.StringIO()):
            boundary = self.namespace["corotation_radius_from_density"](primary, 1., 1.)
        self.assertGreater(boundary, 3.9)
        self.assertLess(boundary, 4.1)
        for cut, message in ((1e-4, "No outward"), (1e3, "already below")):
            with self.assertRaisesRegex(ValueError, message):
                self.namespace["corotation_radius_from_density"](primary, 1., cut)

    def test_taper_preserves_background_and_has_no_linear_kick(self):
        data = dict(
            pos=np.array([[.5, 0., 0.], [0., 1., 0.], [-1.5, 0., 0.],
                          [0., -2., 0.], [.5, 5., 0.]]),
            mass=np.array([2., 3., 1., 4., 1.]),
            vel=np.arange(15, dtype=float).reshape(5, 3), type=np.array([0, 1, 0, 0, 0]))
        primary = Component(data)
        inherited = primary.get_velocity().copy()
        omega = 7.
        real_namespace = dict(np=np, rsol=RSUN)
        execute_nodes([node for node in self.tree.body if isinstance(node, ast.FunctionDef)],
                      real_namespace, BUILDER)
        with contextlib.redirect_stdout(io.StringIO()):
            real_namespace["apply_donor_rotation"](primary, omega, 1., 2.,
                                                  profile="taper", spin_mode="add")
        addition = primary.get_velocity() - inherited
        weights = np.array([1., 1., .5, 0., 0.])
        center = np.average(data["pos"], axis=0, weights=data["mass"] * weights)
        np.testing.assert_allclose(
            addition, weights[:, None] * np.cross([0., 0., omega], data["pos"] - center), atol=1e-14)
        np.testing.assert_allclose(np.einsum("i,ij->j", data["mass"], addition), 0., atol=1e-14)
        np.testing.assert_array_equal(primary.get_velocity()[3:], inherited[3:])
        np.testing.assert_array_equal(primary.get_position(), data["pos"])
        # Core and inner gas follow the same rigid rotation, retaining their residual flow.
        np.testing.assert_allclose(addition[1] - addition[0],
                                   np.cross([0., 0., omega], data["pos"][1] - data["pos"][0]))

    def test_invalid_taper_radii_do_not_change_velocities(self):
        for inner, outer in ((-1., 2.), (1., 1.), (2., 1.), (np.nan, 2.), (1., np.inf)):
            with self.subTest(inner=inner, outer=outer):
                primary = Component(donor_data())
                inherited = primary.get_velocity().copy()
                real_namespace = dict(np=np, rsol=RSUN)
                execute_nodes([node for node in self.tree.body if isinstance(node, ast.FunctionDef)],
                              real_namespace, BUILDER)
                with self.assertRaisesRegex(ValueError, "radii"):
                    real_namespace["apply_donor_rotation"](primary, 1., inner, outer)
                np.testing.assert_array_equal(primary.get_velocity(), inherited)


def fitted_rotation(position, velocity, mass):
    """Independent least-squares fit, excluding the selected region's translation."""
    position = position - np.average(position, axis=0, weights=mass)
    velocity = velocity - np.average(velocity, axis=0, weights=mass)
    design = np.cross(np.eye(3)[:, None, :], position[None, :, :])
    design = design.transpose(1, 2, 0).reshape(-1, 3)
    weight = np.sqrt(np.repeat(mass / mass.sum(), 3))
    return np.linalg.lstsq(design * weight[:, None], velocity.ravel() * weight, rcond=None)[0]


class RotationPhysicsTests(unittest.TestCase):
    def setUp(self):
        self.namespace = dict(np=np, rsol=RSUN)
        execute_nodes([node for node in ast.parse(BUILDER.read_text()).body
                       if isinstance(node, ast.FunctionDef)], self.namespace, BUILDER)

    def test_wind_rigid_inner_and_constant_outer_specific_angular_momentum(self):
        radius = np.geomspace(.1, 100., 2001)
        weight = self.namespace["rotation_profile_weights"](radius, 1., 4.)
        np.testing.assert_array_equal(weight[radius <= 1.], 1.)
        np.testing.assert_allclose(weight[radius >= 4.] * radius[radius >= 4.]**2,
                                   4., rtol=1e-13)
        self.assertGreaterEqual(np.diff(radius**2 * weight).min(), -1e-12)
        self.assertLessEqual(np.diff(weight).max(), 1e-12)
        self.assertGreater(weight[-1], 0.)

    def test_wind_profile_has_smooth_endpoint_slopes(self):
        function = self.namespace["rotation_profile_weights"]
        step = 1e-5
        for boundary, derivative in ((1., 0.), (4., -.125)):
            sample = function(np.array([boundary-step, boundary, boundary+step]), 1., 4.)
            left_slope, right_slope = np.diff(sample) / step
            self.assertAlmostEqual(left_slope, derivative, delta=2e-6)
            self.assertAlmostEqual(right_slope, derivative, delta=2e-6)

    def test_synchronize_removes_tilted_spin_but_preserves_radial_pulsations(self):
        inner = np.r_[np.eye(3), -np.eye(3)] * .5
        position = np.r_[inner, inner * 6., inner * 16., [[0., 0., 0.]]]
        mass = np.r_[np.ones(6), np.full(6, .1), np.full(6, .01), 4.]
        weights = self.namespace["rotation_profile_weights"](np.linalg.norm(position, axis=1), 1., 2.)
        inherited_spin = np.array([.3, -.2, .7])
        radial_pulsations = position * np.repeat([.4, -.6, .3], [6, 6, 7])[:, None]
        bulk = np.array([11., -7., 3.])
        velocity = radial_pulsations + weights[:, None] * np.cross(inherited_spin, position) + bulk
        data = dict(pos=position, mass=mass, vel=velocity, type=np.r_[np.zeros(18), 1],
                    id=np.arange(19))
        primary = Component(data)
        with contextlib.redirect_stdout(io.StringIO()):
            self.namespace["apply_donor_rotation"](primary, 1.4, 1., 2., spin_mode="synchronize")
        target_rotation = weights[:, None] * np.cross([0., 0., 1.4], position)
        np.testing.assert_allclose(primary.get_velocity() - target_rotation - bulk,
                                   radial_pulsations, atol=1e-13)
        np.testing.assert_allclose(fitted_rotation(position[:6], primary.get_velocity()[:6], mass[:6]),
                                   [0., 0., 1.4], atol=1e-13)
        np.testing.assert_array_equal(primary.get_velocity()[-1], velocity[-1])

    def test_asymmetric_donor_no_bulk_kick_residuals_and_other_fields_preserved(self):
        position = np.array([[.5, 0., 0.], [0., .6, 0.], [0., 0., .7],
                             [-.4, -.3, -.1], [2., .8, .1], [-3., 1., .2]])
        mass = np.array([2., 3., 1., 4., .2, .3])
        position -= np.average(position, axis=0, weights=mass)
        velocity = np.array([[.1, .5, .7], [-.5, -.1, .8], [.7, .1, -.3],
                             [.3, -.6, -.4], [.1, .2, .5], [-.1, -.3, .1]])
        data = dict(pos=position, mass=mass, vel=velocity, id=np.arange(6),
                    type=np.array([0, 0, 0, 1, 0, 0]), u=np.arange(5.), rho=np.arange(5.) + 1.)
        primary = Component(data)
        inner = np.linalg.norm(position, axis=1) <= 1.
        old_spin = fitted_rotation(position[inner], velocity[inner], mass[inner])
        with contextlib.redirect_stdout(io.StringIO()):
            self.namespace["apply_donor_rotation"](primary, 2.1, 1., 3., spin_mode="synchronize")
        addition = primary.get_velocity() - velocity
        np.testing.assert_allclose(np.average(addition, axis=0, weights=mass), 0., atol=1e-14)
        np.testing.assert_allclose(fitted_rotation(position[inner], primary.get_velocity()[inner], mass[inner]),
                                   [0., 0., 2.1], atol=1e-13)
        def residual(vel, spin):
            return vel[inner] - np.average(vel[inner], axis=0, weights=mass[inner]) - np.cross(
                spin, position[inner] - np.average(position[inner], axis=0, weights=mass[inner]))
        np.testing.assert_allclose(residual(primary.get_velocity(), [0., 0., 2.1]),
                                   residual(velocity, old_spin), atol=1e-13)
        # Inner gas/core relative radial motion survives the change in mean spin.
        separation = position[0] - position[3]
        self.assertAlmostEqual(np.dot(primary.get_velocity()[0] - primary.get_velocity()[3], separation),
                               np.dot(velocity[0] - velocity[3], separation), delta=1e-13)
        for key in data.keys() - {"vel"}:
            np.testing.assert_array_equal(primary.data[key], data[key])
        synchronized = primary.get_velocity().copy()
        with contextlib.redirect_stdout(io.StringIO()):
            self.namespace["apply_donor_rotation"](primary, 2.1, 1., 3., spin_mode="synchronize")
        np.testing.assert_allclose(primary.get_velocity(), synchronized, atol=1e-13)

    def test_spin_fit_and_application_across_chunk_boundary(self):
        count = 100003
        axes = np.r_[np.eye(3), -np.eye(3)] * .5
        position = axes[np.arange(count) % 6]
        position[-3:] *= 8.
        mass = np.ones(count)
        position -= np.average(position, axis=0, weights=mass)
        old_spin = np.array([.3, -.2, .7])
        velocity = np.cross(old_spin, position) + .4 * position + np.array([3., -1., 7.])
        data = dict(pos=position, vel=velocity, mass=mass)
        primary = Component(data)
        weight = self.namespace["rotation_profile_weights"](np.linalg.norm(position, axis=1), 1., 2.)
        center = np.average(position, axis=0, weights=mass * weight)
        expected = velocity + weight[:, None] * np.cross(np.array([0., 0., 1.4]) - old_spin,
                                                         position - center)
        with contextlib.redirect_stdout(io.StringIO()):
            self.namespace["apply_donor_rotation"](primary, 1.4, 1., 2., spin_mode="synchronize")
        np.testing.assert_allclose(primary.get_velocity(), expected, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(np.average(primary.get_velocity() - velocity, axis=0, weights=mass),
                                   0., atol=1e-12)


@unittest.skipUnless(LIBRARY.is_file(), "sibling stellar_ics source is unavailable")
class OrbitalPhysicsTests(unittest.TestCase):
    def test_rotation_preserves_residuals_and_binary_momentum(self):
        tree = ast.parse(LIBRARY.read_text())
        system = next(node for node in tree.body
                      if isinstance(node, ast.ClassDef) and node.name == "MultipleSystem")
        namespace = dict(np=np, G=G, rsol=RSUN, msol=MSUN,
                         roche_distance=lambda q: .4)
        builder_tree = ast.parse(BUILDER.read_text())
        execute_nodes([node for node in builder_tree.body if isinstance(node, ast.FunctionDef)],
                      namespace, BUILDER)
        methods = {"add_component", "add_components_as_binary"}
        execute_nodes([node for node in system.body
                       if isinstance(node, ast.FunctionDef) and node.name in methods],
                      namespace, LIBRARY)
        Binary = type("Binary", (), {name: namespace[name] for name in methods})

        for factor in (0., .35, 1.):
            with self.subTest(factor=factor), contextlib.redirect_stdout(io.StringIO()):
                primary, secondary = Component(donor_data()), point_mass(1e33)
                inherited = primary.get_velocity().copy()
                binary = Binary()
                binary._components = []
                binary.add_components_as_binary(
                    primary, secondary, distance_fraction_rlof=1.7,
                    corotating_at_rlof=False, corotation_factor=0., e=0.)
                separation = binary.orbital_separation
                omega = np.array([0., 0., np.sqrt(
                    G * (primary.mass + secondary.mass) / separation**3)])
                namespace["apply_donor_rotation"](primary, factor * omega[2], 1e13, 2e13,
                                                   spin_mode="add")
                added_spin = factor * np.cross(omega, primary.get_position())
                orbital_velocity = np.cross(omega, primary.offset)
                np.testing.assert_allclose(
                    primary.get_velocity() - orbital_velocity - added_spin,
                    inherited, rtol=1e-12, atol=1e-8)
                primary_mean = np.average(primary.get_velocity(), axis=0,
                                          weights=primary.data["mass"])
                np.testing.assert_allclose(
                    secondary.get_velocity()[0] - primary_mean,
                    np.cross(omega, secondary.offset - primary.offset), atol=1e-8)
                np.testing.assert_allclose(
                    (primary.mass * primary_mean + secondary.mass * secondary.get_velocity()[0])
                    / (primary.mass + secondary.mass), 0., atol=1e-8)
                np.testing.assert_allclose(
                    secondary.get_velocity()[0], np.cross(omega, secondary.offset), atol=1e-8)


if __name__ == "__main__":
    unittest.main()
