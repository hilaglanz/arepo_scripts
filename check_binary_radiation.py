"""Read-only thermal/radiation diagnostics for a CGS binary IC and snapshot 0.

Loads selected gas fields through loadmodules; writes no snapshots or ICs.
Example: python check_binary_radiation.py binary_1.7.dat.ic.hdf5 output/snapshot_000.hdf5
IC-only comparison: python check_binary_radiation.py new.ic.hdf5 --compare-ic old.ic.hdf5
"""

import argparse
from pathlib import Path

from loadmodules import gadget_readsnapname, np


# Match IDORT's clight and rad_const (src/allvars.h).
LIGHT_SPEED = 2.99792458e10
RADIATION_CONSTANT = 7.5657e-15
WATCH_IDS = {25913384, 6910363, 2606999}  # Printed in the NYX failure log.


def summarize(label, values, ids, selected_ids):
    """Print the distribution and retain IDs of extrema, without sorting all cells."""
    finite = np.isfinite(values)
    bad_indices = np.flatnonzero(~finite)
    selected_ids.update(int(i) for i in ids[bad_indices[:3]])
    print(f"\n{label}: nonfinite={bad_indices.size}, "
          f"negative={np.count_nonzero(values < 0)}, zero={np.count_nonzero(values == 0)}")
    indices = np.flatnonzero(finite)
    if indices.size == 0:
        return
    finite_values = values[indices]
    print("  min / median / p99 / max:", np.percentile(finite_values, [0, 50, 99, 100]))
    lowest = indices[np.argmin(finite_values)]
    count = min(3, indices.size)
    highest = indices[np.argpartition(finite_values, -count)[-count:]]
    highest = highest[np.argsort(values[highest])[::-1]]
    selected_ids.update(int(i) for i in ids[np.append(highest, lowest)])
    print("  largest (ID, value):", [(int(ids[i]), float(values[i])) for i in highest])


def inspect_ic(path, selected_ids):
    print("\nIC:", Path(path).resolve())
    ic = gadget_readsnapname(path, loadonly=["id", "mass", "u"], loadonlytype=[0], nommap=True)
    print("Gas cells:", len(ic.id), "  Total gas mass [g]:", ic.mass.sum())
    summarize("IC u [erg/g]", ic.u, ic.id, selected_ids)
    energy = ic.mass * ic.u
    print("Sum(m*u) [erg]:", energy.sum())
    summarize("IC m*u per cell [erg]", energy, ic.id, selected_ids)
    print("This is the energy encoded in InternalEnergy; it may include restored EOS radiation.")
    return ic


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ic", help="The IC actually used by the run being checked")
    parser.add_argument("snapshot0", nargs="?", help="Optional snapshot 0 from that same run")
    parser.add_argument("--compare-ic", help="Earlier IC to summarize independently; IDs are not paired across builds")
    args = parser.parse_args()
    if args.compare_ic:
        previous = inspect_ic(args.compare_ic, set())
        del previous  # Release the previous IC before loading the new one.
        print("\nOld/new summaries are independent; filtering may change counts and particle IDs.")
    selected_ids = set(WATCH_IDS)
    ic = inspect_ic(args.ic, selected_ids)
    if args.snapshot0 is None:
        print("\nIC-only check assumes CGS units and physical gas masses, as in the NYX run.")
        print("No temperature or radiation field is inferred from IC InternalEnergy.")
        return

    print("\nSnapshot 0:", Path(args.snapshot0).resolve())
    del ic.data["mass"]  # The matched-cell table only needs the IC IDs and energies.
    snap = gadget_readsnapname(args.snapshot0,
                              loadonly=["id", "rho", "temp", "u", "erad", "frad"],
                              loadonlytype=[0], nommap=True)
    units = (snap.UnitLength_in_cm, snap.UnitMass_in_g, snap.UnitVelocity_in_cm_per_s)
    print("Snapshot units (length, mass, velocity):", units)
    if (units != (1.0, 1.0, 1.0) or snap.parameters.get("ComovingIntegrationOn") != 0
            or snap.hubbleparam != 1.0):
        raise ValueError("This diagnostic assumes the noncosmological CGS-unit=1 NYX run.")
    for name, data, required in (("IC", ic.data, {"id", "u"}),
                                 ("snapshot 0", snap.data, {"id", "rho", "temp", "u", "erad", "frad"})):
        missing = required.difference(data)
        if missing:
            raise ValueError(f"Missing {name} fields: {sorted(missing)}")

    print(f"Gas counts: IC={len(ic.id)}, snapshot0={len(snap.id)}; snapshot time={snap.time:.17g}")
    print("IC u may include restored radiation; do not assume it should equal snapshot u.")
    for field in ("rho", "temp", "u", "erad"):
        summarize(f"snapshot0 {field}", snap.data[field], snap.id, selected_ids)
    print("Snapshot0 cells with T >= 1e6 K:", np.count_nonzero(snap.temp >= 1e6))
    print("Snapshot0 cells with T >= 1e8 K:", np.count_nonzero(snap.temp >= 1e8))

    # Use the vector flux magnitude, not just the signed radial component.
    flux = np.hypot.reduce(snap.data.pop("frad"), axis=1)
    summarize("snapshot0 |F_rad| [erg/cm^2/s]", flux, snap.id, selected_ids)
    positive_energy = np.isfinite(snap.erad) & (snap.erad > 0)
    flux_factor = np.full(len(snap.id), np.nan)
    np.divide(flux, LIGHT_SPEED * snap.erad, out=flux_factor, where=positive_energy)
    print("Flux factor undefined where E_rad is nonpositive/nonfinite:",
          np.count_nonzero(~positive_energy))
    summarize("snapshot0 |F_rad|/(c E_rad)", flux_factor, snap.id, selected_ids)
    print("Finite flux factors > 1 + 1e-6:",
          np.count_nonzero(np.isfinite(flux_factor) & (flux_factor > 1.0 + 1e-6)))

    # Read optional opacity fields one at a time rather than keeping more full arrays.
    for field in ("ka_r", "ka_p", "ka_s"):
        opacity = gadget_readsnapname(args.snapshot0, loadonly=[field], loadonlytype=[0], nommap=True)
        if field in opacity.data:
            summarize(f"snapshot0 {field} [cm^2/g]", opacity.data[field], snap.id, selected_ids)
        else:
            print(f"snapshot0 {field}: not stored")
        del opacity

    # Only match the handful of diagnostic IDs; do not assume equal row ordering.
    targets = np.array(sorted(selected_ids), dtype=np.uint64)
    ic_rows = np.flatnonzero(np.isin(ic.id, targets))
    snap_rows = np.flatnonzero(np.isin(snap.id, targets))
    ic_energy = {int(ic.id[i]): float(ic.u[i]) for i in ic_rows}
    found = {int(snap.id[i]) for i in snap_rows}
    print("\nSelected IDs absent from snapshot0:", sorted(selected_ids - found))
    print("Selected IDs absent from IC:", sorted(selected_ids - ic_energy.keys()))
    print("\nID  u_IC  rho_0  T_0  u_0  E_rad_0  T_rad_0  |F_rad|_0  flux_factor")
    for i in snap_rows:
        particle_id = int(snap.id[i])
        radiation_temperature = ((snap.erad[i] / RADIATION_CONSTANT)**0.25
                                 if positive_energy[i] else np.nan)
        values = (ic_energy.get(particle_id, np.nan), snap.rho[i], snap.temp[i],
                  snap.u[i], snap.erad[i], radiation_temperature, flux[i], flux_factor[i])
        print(particle_id, *(f"{value:.6e}" for value in values))
    print("\nT_rad need not equal gas T outside LTE. InternalEnergy is not conserved total Energy.")
    print("A missing failure ID may have been created by refinement after snapshot0.")
    print("Filtering renumbers IC IDs: match IDs within each IC/snapshot pair, not across different IC builds.")


if __name__ == "__main__":
    main()
