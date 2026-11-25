"""Post-process single-star CGS snapshots through loadmodules; never write snapshots.

Run separately for each branch, including snapshot 291 as the reservoir baseline.
Reads gas fields lazily, one at a time; only keeps the selected radial region.
Gamma is kappa*F_radial/(c*G*M_enclosed/r**2), a spherical force diagnostic.
"""
import argparse
import csv
import inspect
import json
from pathlib import Path
import sys

import numpy as np

RSUN = 6.955e10
MSUN = 1.989e33
YEAR = 31536000.0
GRAVITY = 6.6738e-8
LIGHT_SPEED = 2.99792458e10
PARAMETERS = {
    "DustCondensationTemp": ("temperature", 1100.0),
    "BowenOpacitydT": ("width", 50.0),
    "BowenOpacityKappamax": ("kappa_max", 14.821),
    "BowenOpacityRhomin": ("rho_min", 1e-15),
}


def take_region(data, center, core_mass, reference_velocity, rmin, rmax,
                block_size=250000, allow_zero_scatter=False):
    """Consume a gas-only data mapping, including the reader's lazy mapping.

    Interior gas contributes to enclosed gravity even below the plotted region.
    Each full field is accessed once and then released from the mapping.
    """
    required = {"pos", "mass", "vel", "rho", "temp", "ka_r", "frad"}
    if not allow_zero_scatter:
        required.add("ka_s")
    missing = required.difference(data.keys())
    if missing:
        raise ValueError(f"Missing gas fields: {sorted(missing)}. Gamma cannot be inferred from absent fields.")

    positions = data["pos"]
    radius = np.empty(len(positions), dtype=np.float64)
    for first in range(0, len(positions), block_size):
        stop = min(first + block_size, len(positions))
        offset = positions[first:stop].astype(np.float64) - center
        radius[first:stop] = np.linalg.norm(offset, axis=1)
        del offset
    if not np.isfinite(radius).all():
        raise ValueError("Nonfinite gas positions.")
    selected = np.flatnonzero((radius >= rmin) & (radius < rmax))
    selected = selected[np.argsort(radius[selected], kind="stable")]
    distance = radius[selected]
    direction = positions[selected].astype(np.float64, copy=False) - center
    direction /= distance[:, None]
    del data["pos"], positions

    masses = data["mass"]
    if not np.isfinite(masses).all() or np.any(masses < 0):
        raise ValueError("Gas masses must be finite and nonnegative.")
    interior_mass = core_mass + masses[radius < rmin].sum(dtype=np.float64)
    cells = {"radius": distance, "mass": masses[selected].astype(np.float64, copy=False)}
    cells["enclosed_mass"] = interior_mass + np.cumsum(cells["mass"], dtype=np.float64)
    del data["mass"], masses, radius

    for field in ("vel", "frad", "rho", "temp", "ka_r", "ka_s"):
        if field == "ka_s" and field not in data:
            cells[field] = np.zeros(len(selected))
            continue
        full_field = data[field]  # Trigger lazy loading before deleting this key.
        values = full_field[selected].astype(np.float64, copy=False)
        del data[field], full_field
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite {field} values in the selected radial region.")
        if field == "vel":
            cells["vr"] = np.einsum("ij,ij->i", values - reference_velocity, direction)
        elif field == "frad":
            cells["flux_radial"] = np.einsum("ij,ij->i", values, direction)
        else:
            cells[field] = values
        del values
    if np.any(cells["rho"] <= 0) or np.any(cells["temp"] <= 0):
        raise ValueError("Density and temperature must be positive in the selected region.")
    return cells


def add_force_diagnostics(cells, law):
    """The trial opacity changes neither the radiation field nor gas state."""
    dense = cells["rho"] > law["rho_min"]
    exponent = np.clip((cells["temp"] - law["temperature"]) / law["width"], -700, 700)
    cells["kappa_trial"] = np.where(dense, law["kappa_max"] / (1 + np.exp(exponent)), 0.0)
    gravity = GRAVITY * cells["enclosed_mass"] / cells["radius"] ** 2
    force_per_opacity = cells["flux_radial"] / (LIGHT_SPEED * gravity)
    cells["gamma_base"] = cells["ka_r"] * force_per_opacity
    cells["gamma_added"] = cells["ka_s"] * force_per_opacity
    cells["gamma_total"] = cells["gamma_base"] + cells["gamma_added"]
    # Replace the stored added opacity for this counterfactual, never add it twice.
    cells["gamma_trial_total"] = cells["gamma_base"] + cells["kappa_trial"] * force_per_opacity


def summarize(cells, edges, law):
    """Mass-weighted diagnostics; inward/outward fluxes are positive magnitudes."""
    count = len(edges) - 1
    bins = np.searchsorted(edges, cells["radius"], side="right") - 1
    included = (bins >= 0) & (bins < count)
    bins = bins[included]
    mass = cells["mass"][included]

    def total(values):
        return np.bincount(bins, weights=values, minlength=count)

    def divide(numerator, denominator):
        return np.divide(numerator, denominator, out=np.full(count, np.nan), where=denominator > 0)

    shell_mass = total(mass)
    means = lambda values: divide(total(mass * values[included]), shell_mass)
    inward_cells = cells["vr"][included] < 0
    inward_mass = total(mass * inward_cells)
    outward = total(mass * np.maximum(cells["vr"][included], 0)) / np.diff(edges) * YEAR / MSUN
    inward = total(mass * np.maximum(-cells["vr"][included], 0)) / np.diff(edges) * YEAR / MSUN
    dense = cells["rho"] > law["rho_min"]
    columns = {
        "inner_Rsun": edges[:-1] / RSUN,
        "outer_Rsun": edges[1:] / RSUN,
        "radius_Rsun": np.sqrt(edges[:-1] * edges[1:]) / RSUN,
        "cell_count": np.bincount(bins, minlength=count),
        "gas_mass_Msun": shell_mass / MSUN,
        "density_volume_mean_g_cm3": divide(shell_mass, total(mass / cells["rho"][included])),
        "temperature_mass_mean_K": means(cells["temp"]),
        "radial_velocity_mass_mean_km_s": means(cells["vr"]) / 1e5,
        "outward_Msun_per_yr": outward,
        "inward_Msun_per_yr": inward,
        "net_Msun_per_yr": outward - inward,
        "inward_mass_fraction": divide(inward_mass, shell_mass),
        "density_eligible_mass_fraction": means(dense),
        "cool_and_dense_mass_fraction": means(dense & (cells["temp"] <= law["temperature"])),
        "stored_added_kappa_mass_mean_cm2_g": means(cells["ka_s"]),
        "trial_kappa_mass_mean_cm2_g": means(cells["kappa_trial"]),
        "halfmax_added_opacity_mass_fraction": means(cells["ka_s"] >= 0.5 * law["kappa_max"]),
        "gamma_base_mass_mean": means(cells["gamma_base"]),
        "gamma_added_mass_mean": means(cells["gamma_added"]),
        "gamma_total_mass_mean": means(cells["gamma_total"]),
        "gamma_trial_total_mass_mean": means(cells["gamma_trial_total"]),
        "gamma_total_gt1_mass_fraction": means(cells["gamma_total"] > 1),
        "fraction_of_inward_mass_gamma_total_gt1": divide(
            total(mass * inward_cells * (cells["gamma_total"][included] > 1)), inward_mass),
        "fraction_of_inward_mass_gamma_trial_total_gt1": divide(
            total(mass * inward_cells * (cells["gamma_trial_total"][included] > 1)), inward_mass),
    }
    return [{key: float(values[i]) for key, values in columns.items()} for i in range(count)]


def read_snapshot(args, number):
    from loadmodules import gadget_readsnap

    options = dict(hdf5=True, cosmological=False, applytransformationfacs=False)
    core = gadget_readsnap(number, str(args.source_dir), loadonlytype=[1],
                           loadonly=["id", "pos", "vel", "mass"], **options)
    if len(core.mass) != 1:
        raise ValueError("This diagnostic requires a single type-1 stellar core.")
    units = [core.get_header_attribute(key) for key in
             ("UnitLength_in_cm", "UnitMass_in_g", "UnitVelocity_in_cm_per_s")]
    if units != [1.0, 1.0, 1.0] or core.get_header_attribute("ComovingIntegrationOn") != 0 or core.hubbleparam != 1:
        raise ValueError(f"Requires noncosmological CGS snapshots with unit factors and HubbleParam = 1; got {units}.")
    if not core.config:
        raise ValueError("Missing /Config metadata: cannot establish fixed-core or Bowen status.")
    # Reject binaries: spherical single-star gravity is not their full force field.
    if np.any(np.asarray(core.header["NumPart_Total"])[2:] != 0):
        raise ValueError("Use a single-star snapshot containing only gas and the type-1 core.")
    fixed = "DMFIXED" in core.config
    reference_velocity = np.zeros(3) if fixed else core.vel[0].astype(np.float64)
    law, parameter_sources = {}, {}
    for parameter, (name, default) in PARAMETERS.items():
        law[name] = float(core.parameters.get(parameter, getattr(args, name)))
        parameter_sources[parameter] = "snapshot" if parameter in core.parameters else "CLI fallback"
    if law["width"] <= 0 or law["temperature"] <= 0 or law["kappa_max"] <= 0 or law["rho_min"] < 0:
        raise ValueError("Invalid diagnostic Bowen parameters.")
    gas = gadget_readsnap(number, str(args.source_dir), loadonlytype=[0], lazy_load=True,
                          loadonly=["pos", "vel", "mass", "rho", "temp", "ka_r", "ka_s", "frad"], **options)
    zero_scatter = "ka_s" not in gas.data and "BOWEN_OPACITY" not in core.config
    metadata = {
        "snapshot": number, "time_s": float(core.time), "time_yr": float(core.time) / YEAR,
        "source_dir": str(args.source_dir.resolve()), "core_id": int(core.id[0]),
        "python_interpreter": sys.executable, "reader_path": inspect.getfile(gadget_readsnap),
        "core_position_cm": core.pos[0].tolist(), "stored_core_velocity_cm_s": core.vel[0].tolist(),
        "reference_velocity_cm_s": reference_velocity.tolist(),
        "fixed_core": fixed, "bowen_enabled": "BOWEN_OPACITY" in core.config,
        "idort_neglect_vel": "IDORT_NEGLECT_VEL" in core.config,
        "missing_scatter_explicitly_zero": zero_scatter,
        "bowen_law": law, "parameter_sources": parameter_sources,
    }
    print(f"Snapshot {number}: t={metadata['time_yr']:.6f} yr, Bowen={metadata['bowen_enabled']}, fixed_core={fixed}")
    print(f"  Diagnostic opacity law: {law}; sources: {parameter_sources}")
    cells = take_region(gas.data, core.pos[0], float(core.mass[0]), reference_velocity,
                        args.rmin * RSUN, args.rmax * RSUN, args.block_size, zero_scatter)
    add_force_diagnostics(cells, law)
    return cells, metadata


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_profiles(rows, output):
    import matplotlib.pyplot as plt
    plt.switch_backend("Agg")
    panels = [
        ("gas_mass_Msun", "Gas mass per radial bin [solar masses]", "log"),
        ("density_volume_mean_g_cm3", "Mean density [g / cm3]", "log"),
        ("temperature_mass_mean_K", "Mass-weighted gas temperature [K]", "log"),
        ("radial_velocity_mass_mean_km_s", "Core-relative radial velocity [km / s]", "linear"),
        ("stored_added_kappa_mass_mean_cm2_g", "Stored added opacity [cm2 / g]", "symlog"),
        ("gamma_total_mass_mean", "Mass-weighted total Gamma (force proxy)", "symlog"),
        ("net_Msun_per_yr", "Net shell flux [solar masses / yr]", "symlog"),
        ("fraction_of_inward_mass_gamma_total_gt1", "Fraction of infalling mass with Gamma > 1", "linear"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(12, 14), sharex=True, constrained_layout=True)
    for ax, (field, label, scale) in zip(axes.flat, panels):
        for number in sorted({row["snapshot"] for row in rows}):
            subset = [row for row in rows if row["snapshot"] == number]
            ax.plot([row["radius_Rsun"] for row in subset], [row[field] for row in subset],
                    label=f"{number}: {subset[0]['time_yr']:.3f} yr")
        ax.set_xscale("log")
        if scale == "symlog":
            ax.set_yscale(scale, linthresh=1e-8 if field == "net_Msun_per_yr" else 0.01)
        else:
            ax.set_yscale(scale)
        ax.set_ylabel(label)
        ax.grid(alpha=0.2)
        if field == "gamma_total_mass_mean":
            ax.axhline(1, color="black", linewidth=0.7, linestyle="--")
        if field in ("net_Msun_per_yr", "radial_velocity_mass_mean_km_s"):
            ax.axhline(0, color="black", linewidth=0.7)
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Distance from core [solar radii]")
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--label", required=True, help="For example bowen or no_bowen")
    parser.add_argument("--snapshots", type=int, nargs="+", default=[291, 413, 479, 496])
    parser.add_argument("--output-dir", type=Path, default=Path("bowen_launch_analysis"))
    parser.add_argument("--rmin", type=float, default=300, help="Solar radii")
    parser.add_argument("--rmax", type=float, default=3000, help="Solar radii")
    parser.add_argument("--bins", type=int, default=48)
    parser.add_argument("--block-size", type=int, default=250000)
    parser.add_argument("--plot", action="store_true")
    for parameter, (name, default) in PARAMETERS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=float, default=default,
                            help=f"Fallback {parameter} if absent from snapshot; default {default:g}")
    args = parser.parse_args()
    if not 0 < args.rmin < args.rmax or args.bins < 1 or args.block_size < 1:
        parser.error("Require 0 < rmin < rmax, bins >= 1 and block-size >= 1.")
    if Path(args.label).name != args.label or args.label in (".", ".."):
        parser.error("label must be a filename label without a directory.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fine_edges = np.geomspace(args.rmin, args.rmax, args.bins + 1) * RSUN
    broad_edges = np.unique([args.rmin, *[r for r in [500, 700, 1000, 1500, 2000]
                                        if args.rmin < r < args.rmax], args.rmax]) * RSUN
    profiles, shells, metadata, baseline = [], [], [], {}
    for number in sorted(set(args.snapshots)):
        cells, info = read_snapshot(args, number)
        for name, edges, output in (("profiles", fine_edges, profiles), ("shells", broad_edges, shells)):
            rows = summarize(cells, edges, info["bowen_law"])
            for row in rows:
                key = (name, row["inner_Rsun"], row["outer_Rsun"])
                if number == 291:
                    baseline[key] = row["gas_mass_Msun"]
                row["delta_gas_mass_from_291_Msun"] = row["gas_mass_Msun"] - baseline.get(key, np.nan)
                row.update(label=args.label, snapshot=number, time_s=info["time_s"], time_yr=info["time_yr"])
            output.extend(rows)
            write_csv(args.output_dir / f"{args.label}_{name}.csv", output)
        metadata.append(info)
        (args.output_dir / f"{args.label}_metadata.json").write_text(json.dumps({
            "notes": [
                "CGS, single star; spherical gravity includes core and all interior gas.",
                "Gamma is a lab-frame kappa*F/c proxy, not the exact conservative IDORT momentum update or a pressure-inclusive acceleration.",
                "Trial opacity uses the recorded/fallback gas-temperature Bowen law on the unchanged snapshot radiation field; no feedback is predicted.",
                "ka_r is the stored effective base opacity; it may already contain table dust and opacity switching.",
                "Mass changes are Eulerian shell inventories, not tracking individual parcels or proving escape.",
                "Empty-bin means or fractions with no infalling mass are NaN; empty-bin mass and flux are zero.",
            ], "snapshots": metadata,
        }, indent=2), encoding="utf-8")
        del cells
    if args.plot:
        plot_profiles(profiles, args.output_dir / f"{args.label}_profiles.png")
    print(f"Saved {args.label}_profiles.csv, {args.label}_shells.csv and metadata in {args.output_dir}")


if __name__ == "__main__":
    main()
