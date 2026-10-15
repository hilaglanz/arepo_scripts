import argparse
import sys
import os

from loadmodules import *
from stellar_ics.tools import *
from stellar_ics.multiple_star import MultipleSystem, SnapshotComponent, PointMassComponent
from BinariesICs import *
from plot_multiples import get_obj_index
def initialized_new_data(npart, npart0):
    data = {}
    data['pos'] = np.zeros((npart, 3))
    data['vel'] = np.zeros((npart, 3))
    data['mass'] = np.zeros(npart)
    data['u'] = np.zeros(npart0)
    data['count'] = npart
    data['type'] = np.zeros(npart)

    return data

def copy_old_data(snapshot):
    data = initialized_new_data(snapshot.nparticles[0], snapshot.nparticles[0])
    data['pos'] = snapshot.data['pos'][:snapshot.nparticles[0]]
    data['vel']= snapshot.data['vel'][:snapshot.nparticles[0]]
    data['mass'] = snapshot.data['mass'][:snapshot.nparticles[0]]
    data['u'] = snapshot.data['u']
    if 'bfld' in snapshot.data.keys():
        data['bfld'] = snapshot.data['bfld']
    if 'pass' in snapshot.data.keys():
        data['pass'] = snapshot.data['pass']
    data['boxsize'] = snapshot.boxsize

    return data

def remove_bulk_velocity(data):
    """Remove COM velocity while preserving all relative velocities."""
    bulk = np.einsum("i,ij->j", data["mass"], data["vel"]) / data["mass"].sum()
    data["vel"] -= bulk

    return bulk

def rotation_profile_weights(radius, inner_radius, outer_radius, profile="wind"):
    """Angular velocity / inner angular velocity; radii in cm.

    The wind profile smoothly changes d(log Omega)/d(log r) from 0 to -2,
    with nondecreasing specific angular momentum along a radial streamline.
    Outside outer_radius, Omega / Omega_inner = inner_radius*outer_radius/r**2.
    The optional taper profile instead sets the added rotation to zero outside.
    """
    if (not np.isfinite(inner_radius) or not np.isfinite(outer_radius)
            or inner_radius < 0 or outer_radius <= inner_radius):
        raise ValueError("Corotation radii must satisfy 0 <= inner < outer and be finite")
    if profile == "taper":
        transition = np.clip((radius - inner_radius) / (outer_radius - inner_radius), 0.0, 1.0)
        return 1.0 - transition**2 * (3.0 - 2.0 * transition)
    if profile != "wind":
        raise ValueError("corotation_profile must be 'wind' or 'taper'")
    if inner_radius <= 0:
        raise ValueError("The wind rotation profile requires a positive inner radius")
    log_width = np.log(outer_radius / inner_radius)
    log_radius = np.log(np.maximum(radius, inner_radius) / inner_radius)
    transition = np.clip(log_radius / log_width, 0.0, 1.0)
    return np.exp(-2.0 * log_width * (transition**3 - 0.5 * transition**4)
                  - 2.0 * np.maximum(log_radius - log_width, 0.0))


def measure_inner_rotation(component, inner_radius):
    """Mass-weighted rigid-rotation fit within inner_radius, excluding translation.

    Return the best-fitting 3D angular velocity [rad/s], not an individual core
    velocity or an assumption about convection. Work in chunks for large donors.
    """
    position, velocity = component.get_position(), component.get_velocity()
    mass = component.data["mass"]
    total_mass = 0.0
    position_sum, velocity_sum, angular_momentum = np.zeros((3, 3))
    second_moment = np.zeros((3, 3))
    for start in range(0, len(mass), 100000):
        pos, vel = position[start:start+100000], velocity[start:start+100000]
        inner_mass = np.where(np.einsum("ij,ij->i", pos, pos) <= inner_radius**2,
                              mass[start:start+100000], 0.0)
        total_mass += inner_mass.sum()
        position_sum += np.einsum("i,ij->j", inner_mass, pos)
        velocity_sum += np.einsum("i,ij->j", inner_mass, vel)
        second_moment += np.einsum("i,ij,ik->jk", inner_mass, pos, pos)
        angular_momentum += np.einsum("i,ij->j", inner_mass, np.cross(pos, vel))
    if total_mass <= 0:
        raise ValueError("Full-spin region contains no donor mass")
    second_moment -= np.outer(position_sum, position_sum) / total_mass
    angular_momentum -= np.cross(position_sum, velocity_sum) / total_mass
    inertia = np.trace(second_moment) * np.eye(3) - second_moment
    try:
        return np.linalg.solve(inertia, angular_momentum)
    except np.linalg.LinAlgError as error:
        raise ValueError("Full-spin region must contain a spatially extended donor") from error


def apply_donor_rotation(component, angular_frequency, inner_radius, outer_radius,
                         profile="wind", spin_mode="synchronize"):
    """Set mean inner spin along +z or add spin, retaining velocity residuals.

    The same radial weight multiplies the change in fitted angular velocity;
    the old wind's unknown rotation history is retained in the residual field.
    Radial weights use distance from the donor COM. A weighted rotation center
    ensures zero added linear momentum. For an asymmetric donor this introduces
    an offset from that COM, so the concentric wind law is approximate; the offset
    is printed. All retained donor particles, including the core, share the field.
    Existing donor background is also included; this does not identify ambient gas.
    """
    if not np.isfinite(angular_frequency):
        raise ValueError("Corotation angular frequency must be finite")
    # Validate before fitting or changing velocities.
    rotation_profile_weights(np.zeros(1), inner_radius, outer_radius, profile)
    if spin_mode not in ("synchronize", "add"):
        raise ValueError("corotation_spin_mode must be 'synchronize' or 'add'")
    old_spin = (measure_inner_rotation(component, inner_radius)
                if spin_mode == "synchronize" else np.zeros(3))
    spin_change = np.array([0.0, 0.0, angular_frequency]) - old_spin
    position = component.get_position()
    velocity = component.get_velocity()
    mass = component.data["mass"]
    rotating_mass, position_sum = 0.0, np.zeros(3)
    for start in range(0, len(mass), 100000):
        pos = position[start:start+100000]
        radius = np.sqrt(np.einsum("ij,ij->i", pos, pos))
        weighted_mass = mass[start:start+100000] * rotation_profile_weights(
            radius, inner_radius, outer_radius, profile)
        rotating_mass += weighted_mass.sum()
        position_sum += np.einsum("i,ij->j", weighted_mass, pos)
    if rotating_mass <= 0:
        raise ValueError("Corotation region contains no donor mass")
    rotation_center = position_sum / rotating_mass
    maximum_speed_change = 0.0
    for start in range(0, len(mass), 100000):
        pos = position[start:start+100000]
        radius = np.sqrt(np.einsum("ij,ij->i", pos, pos))
        weight = rotation_profile_weights(radius, inner_radius, outer_radius, profile)
        change = weight[:, None] * np.cross(spin_change, pos - rotation_center)
        velocity[start:start+100000] += change
        maximum_speed_change = max(maximum_speed_change, np.sqrt(np.einsum("ij,ij->i", change, change)).max())
    component.set_velocity(velocity)
    print(f"Donor rotation: profile={profile}, mode={spin_mode}; inner/outer radii "
          f"{inner_radius / rsol:g}/{outer_radius / rsol:g} Rsun; "
          f"target/added omega_z={angular_frequency:.8e} rad/s")
    if spin_mode == "synchronize":
        print(f"Fitted inherited omega [rad/s]={old_spin}")
    else:
        print("Existing donor spin retained (add mode)")
    print(f"Applied omega change [rad/s]={spin_change}; "
          f"maximum speed change={maximum_speed_change / 1e5:.8g} km/s")
    print(f"Momentum-preserving rotation-center offset [cm]={rotation_center}; "
          f"offset / inner radius={np.linalg.norm(rotation_center) / max(inner_radius, 1.0):.8g}")
    if profile == "wind":
        print(f"Wind keeps rotation beyond the outer radius (v_phi approximately 1/r); "
              f"effective launch radius={np.sqrt(inner_radius * outer_radius) / rsol:g} Rsun")

def corotation_radius_from_density(component, inner_radius, density_cut):
    """First outward density crossing in 128 logarithmic shells beyond inner_radius.

    Radius is measured from the donor COM [cm]; density is volume averaged
    (sum(mass) / sum(mass/rho)) in g/cm^3. Interpolate the crossing in log-log
    space. Later dense shells do not extend the transition boundary. This sets
    a radial profile parameter, not a cell-by-cell spin mask or gas-removal cut.
    """
    if not np.isfinite(density_cut) or density_cut <= 0:
        raise ValueError("corotation_density_cut must be finite and positive")
    if not np.isfinite(inner_radius) or inner_radius <= 0:
        raise ValueError("Density-based corotation requires a positive finite inner radius")
    gas = component.data["type"] == 0
    position = component.get_position()
    radius = np.sqrt(np.einsum("ij,ij->i", position, position))[gas]
    mass = component.data["mass"][gas]
    density = component.data["rho"]  # Gas-only field, ordered like the gas particles.
    if np.any(~np.isfinite(density) | (density <= 0)):
        raise ValueError("Density-based corotation requires finite positive gas densities")
    if radius.size == 0 or radius.max() <= inner_radius:
        raise ValueError("No donor gas outside the full-spin radius to locate a density crossing")
    edges = np.geomspace(inner_radius, radius.max(), 129)
    shell_mass = np.histogram(radius, bins=edges, weights=mass)[0]
    shell_volume = np.histogram(radius, bins=edges, weights=mass / density)[0]
    occupied = shell_volume > 0
    shell_density = shell_mass[occupied] / shell_volume[occupied]
    shell_radius = np.sqrt(edges[:-1] * edges[1:])[occupied]
    crossings = np.flatnonzero(shell_density <= density_cut)
    if crossings.size == 0:
        raise ValueError("No outward density crossing found; use a higher spin cutoff or explicit outer radius")
    crossing = crossings[0]
    if crossing == 0:
        raise ValueError("Density already below spin cutoff near inner radius; lower cutoff or inner radius")
    outer_radius = np.exp(np.interp(np.log(density_cut),
                                   np.log(shell_density[crossing-1:crossing+1])[::-1],
                                   np.log(shell_radius[crossing-1:crossing+1])[::-1]))
    print(f"Corotation density boundary {density_cut:.8e} g/cm^3: "
          f"estimated outer radius {outer_radius / rsol:.8g} Rsun")
    return outer_radius

def enclosing_boxsize(components, minimum_size, padding):
    """Size of a COM-centered cube enclosing all particle centers."""
    half_extent = 0.0
    for component in components:
        pos = component.get_position()
        lower = pos.min(axis=0) + component.offset
        upper = pos.max(axis=0) + component.offset
        half_extent = max(
            half_extent,
            np.abs(lower).max(),
            np.abs(upper).max(),
        )
    return max(minimum_size, 2.0 * (half_extent + padding))

def filter_donor_background(component, density_cut=1e-20):
    """Filter a CGS SnapshotComponent in place, preserving all nongas particles."""
    if density_cut is None or density_cut == 0:
        return None

    data = component.data
    keep_gas = data["rho"] >= density_cut
    keep_particles = np.ones(len(data["type"]), dtype=bool)
    keep_particles[data["type"] == 0] = keep_gas
    removed_mass = data["mass"][~keep_particles].sum()
    gas_mass = data["mass"][data["type"] == 0].sum()
    for key, values in data.items():
        if not np.isscalar(values):
            mask = keep_gas if len(values) == keep_gas.size else keep_particles
            data[key] = values[mask]
    data.update(count=int(keep_particles.sum()), boxsize=component.boxsize)

    component._mass = component.get_mass()  # Refresh the constructor's cached mass.
    component.recenter()
    print(f"Removed {np.count_nonzero(~keep_gas)} donor gas cells: "
          f"{removed_mass / msol:.6e} Msun, "
          f"{removed_mass / gas_mass:.6e} of donor gas mass")
    return keep_gas


def energy_for_standard_arepo_ic(snapshot, keep_gas=None):
    required = {"EOS_OPAL", "IRT_STAR_PH_EOS_SWITCH"}
    unsupported = {"EOS_NO_RAD", "EOS_SUBTRACT_RAD"}

    if (not required.issubset(snapshot.config)
            or unsupported.intersection(snapshot.config)):
        raise ValueError(
            "This correction requires a donor evolved with switched OPAL EOS, "
            "without EOS_NO_RAD or EOS_SUBTRACT_RAD."
        )

    units = (
        snapshot.UnitLength_in_cm,
        snapshot.UnitMass_in_g,
        snapshot.UnitVelocity_in_cm_per_s,
    )
    if units != (1.0, 1.0, 1.0):
        raise ValueError("This correction currently supports CGS unit=1 snapshots only.")

    if (snapshot.parameters.get("ComovingIntegrationOn") != 0
            or snapshot.hubbleparam != 1.0):
        raise ValueError("This correction requires non-cosmological snapshot data.")
    selection = slice(None) if keep_gas is None else keep_gas
    density = np.asarray(snapshot.data["rho"][selection], dtype=np.float64)
    temperature = np.asarray(snapshot.data["temp"][selection], dtype=np.float64)
    photosphere_density = snapshot.parameters["IrtStarPhDensity"]

    retained_fraction = np.clip(
        (density / photosphere_density - 0.1) / 9.9, 0.0, 1.0
    )
    # Match the constants used in AREPO's opal_eos.h.
    radiation_constant = 4.0 * 5.67051e-5 / 2.9979245e10

    return (
        snapshot.data["u"][selection]
        + (1.0 - retained_fraction)
        * radiation_constant * temperature**4 / density
    )

class BinarySystemWithLogIDs(MultipleSystem):
    def write_ics(self, filename=None):
        if self.parameters.get("binary_log_ids", False):
            types = self.data["type"]
            ids = self.data["id"]
            core, sink = types == 1, types == 5

            if np.count_nonzero(core) != 1 or np.count_nonzero(sink) != 1:
                raise ValueError("Binary-log IDs require one core and one sink")

            reserved = (ids == 1000000000) | (ids == 1000000001)
            if np.any(reserved & ~(core | sink)):
                raise ValueError("A requested binary ID belongs to another particle")

            ids[core] = 1000000000
            ids[sink] = 1000000001

        return super().write_ics(filename=filename)
    
def AddPointMassToFile(snapshot_file, new_file_name, point_mass, separation, rlof_factor=1.0, giant_radius_rsol=None,
                       padding_cm=0, restore_eos_radiation=False,
                       donor_density_cut=1e-20, corotation_factor=0.0,
                       corotation_inner_radius_rsol=None, corotation_outer_radius_rsol=None,
                       corotation_density_cut=None, corotation_profile="wind",
                       corotation_spin_mode="synchronize"):
    """Build a circular binary, optionally adding aligned donor spin.

    A zero corotation_factor disables spin modification. Otherwise it sets the
    target mean inner spin / orbital frequency along +z in synchronize mode, or
    the added spin in add mode. The rigid mean spin is fitted and replaced while
    retaining residual convection/pulsation velocities. The default wind profile
    transitions to approximately conserved specific spin angular momentum outside
    the outer radius; taper instead removes the added rotation outside that radius.
    The inner radius defaults to the giant radius. Supply an outer radius or a
    density crossing to locate the end of the transition (both measured from the
    donor COM). This is an initial velocity prescription, not tidal equilibrium
    or a reconstruction of the old wind's launch history. Existing donor ambient
    gas cannot be identified automatically and receives the same radial profile.
    """
    if not np.isfinite(corotation_factor):
        raise ValueError("corotation_factor must be finite")
    if corotation_profile not in ("wind", "taper"):
        raise ValueError("corotation_profile must be 'wind' or 'taper'")
    if corotation_spin_mode not in ("synchronize", "add"):
        raise ValueError("corotation_spin_mode must be 'synchronize' or 'add'")
    if corotation_density_cut is not None and corotation_outer_radius_rsol is not None:
        raise ValueError("Choose corotation_density_cut or corotation_outer_radius_rsol, not both")
    if corotation_factor != 0.0 and corotation_outer_radius_rsol is None and corotation_density_cut is None:
        raise ValueError("Specify corotation_outer_radius_rsol or corotation_density_cut when spin is enabled")
    if corotation_density_cut is not None and (not np.isfinite(corotation_density_cut) or corotation_density_cut <= 0):
        raise ValueError("corotation_density_cut must be finite and positive")
    snapshot=gadget_readsnapname(snapshot_file)
    new_size = snapshot.boxsize

    giant = SnapshotComponent(data=snapshot.data, boxsize=snapshot.boxsize, radius=giant_radius_rsol)
    keep_gas = filter_donor_background(giant, donor_density_cut)
    # Set ambient properties from retained material before radiation restoration.
    grid_xnuc = giant.data["xnuc"][0].copy()
    grid_u = min(float(giant.data["u"].min()), 1e10)
    grid_rho = min(float(giant.data["rho"].min()), 1e-20)
    if donor_density_cut is not None:
        if grid_rho < donor_density_cut:
            print("rho cut is higher than min density, using rho_cut for grid rho")
        grid_rho = max(grid_rho, donor_density_cut)

    if restore_eos_radiation:
        giant.data["u"] = energy_for_standard_arepo_ic(snapshot, keep_gas)
    del keep_gas
    remove_bulk_velocity(giant.data)
    companion = PointMassComponent(mass=point_mass)
    companion.data['type'] = np.array([5])
    q = giant.mass / point_mass
    if giant_radius_rsol is None:
        giant_radius_rsol = giant.get_radius() / rsol
        print("calculated radius=", giant_radius_rsol)
    current_rlof = giant_radius_rsol  / roche_distance(q)
    print("Roche-filling separation= ", current_rlof, " Rsun")
    if separation is None:
        print("calculating separation from Roche Lobe")
        separation = current_rlof * rlof_factor
        print("separation= ", separation)
    else:
        print("using given separation of ", separation, "Rsun")
        rlof_factor = separation / current_rlof

    print("Roche factor = ", rlof_factor)
    rlof_factor *= (giant_radius_rsol * rsol / giant.get_radius())
    print("rlof_factor according to radius calculation = ", rlof_factor)
    binary = BinarySystemWithLogIDs(newsize=new_size, reset_dm_ids=True, binary_log_ids=True, ndir=32, 
                                    grid_xnuc=grid_xnuc, grid_rho=grid_rho,
                                    grid_u=grid_u)
    # Add only the common orbital velocity here; spin is applied separately.
    binary.add_components_as_binary(giant, companion, distance_fraction_rlof=rlof_factor,
                                    corotating_at_rlof=False, corotation_factor=0.0, e=0.0)
    if corotation_factor != 0.0:
        inner_radius = (giant_radius_rsol if corotation_inner_radius_rsol is None
                        else corotation_inner_radius_rsol) * rsol
        outer_radius = (corotation_radius_from_density(giant, inner_radius, corotation_density_cut)
                        if corotation_density_cut is not None else corotation_outer_radius_rsol * rsol)
        angular_frequency = corotation_factor * np.sqrt(
            G * (giant.mass + companion.mass) / binary.orbital_separation**3)
        apply_donor_rotation(giant, angular_frequency, inner_radius, outer_radius,
                              profile=corotation_profile, spin_mode=corotation_spin_mode)
    binary.newsize = enclosing_boxsize((giant, companion), snapshot.boxsize, padding_cm)
    print(f"Use BoxSize {binary.newsize:.17g} in the AREPO parameter file")
    binary.create_ics(model=new_file_name)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--giant_snapshot_file', type=str, help='path to snapshot file containing the giant star', default="")
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types, '
                                                                  'if there is a point mass companion one '
                                                                  'should also load type 1 or 5', default=[0,1])
    parser.add_argument('--orbital_separation', type=float,
                        help='initial binary separation in Rsun', default=None)
    parser.add_argument('--giant_radius', type=float,
                        help='initial giant radius in Rsun', default=None)
    parser.add_argument('--point_mass', type=float, help='new object mass in msun', default=1)
    parser.add_argument('--rlof_factor', type=float, help='if relative to RL, by what factor?', default=1)
    parser.add_argument('--corotation_factor', type=float, default=0.0,
                        help='Inner spin / orbital angular frequency about +z: '
                             '0 disables modification; 1 synchronizes the mean inner rotation. '
                             'Convection and pulsation residuals are retained. '
                             'Default: %(default)s.')
    parser.add_argument('--corotation_profile', choices=('wind', 'taper'), default='wind',
                        help='wind: smoothly join solid rotation to v_phi approximately 1/r; '
                             'taper: added rotation decreases to zero outside the outer radius. '
                             'Default: %(default)s.')
    parser.add_argument('--corotation_spin_mode', choices=('synchronize', 'add'), default='synchronize',
                        help='synchronize: replace the fitted inner rigid rotation with the requested spin; '
                             'add: retain existing spin and add the requested profile. Default: %(default)s.')
    parser.add_argument('--corotation_inner_radius', type=float, default=None,
                        help='Rigid mean rotation inside this radius [Rsun], measured from the donor COM; '
                             'default: giant radius.')
    spin_boundary = parser.add_mutually_exclusive_group()
    spin_boundary.add_argument('--corotation_outer_radius', type=float, default=None,
                        help='End of transition [Rsun], measured from the donor COM: '
                             'wind keeps approximately constant specific angular momentum beyond it; '
                             'taper reaches zero added spin. Must exceed inner radius.')
    spin_boundary.add_argument('--corotation_density_cut', type=float, default=None,
                        help='Find the outer transition radius at the first shell-volume-mean density crossing '
                             'below this threshold [g/cm^3], beyond the inner radius. '
                             'Does not remove gas. Supply this or corotation_outer_radius to enable spin.')
    parser.add_argument("--box_padding_rsun", type=float,
                        help="Clearance beyond current retained particle centers, in solar radii", default=0.0)
    parser.add_argument("--restore_eos_radiation", action=argparse.BooleanOptionalAction,
                        help="Prepare evolved switched-OPAL energies for normal AREPO IC startup (off by default).",
                        default=False)
    parser.add_argument("--donor_density_cut", type=float, default=0,
                        help="Remove donor gas below this density [g/cm^3]; also the minimum added "
                             "background density. Default: %(default)s; 0 disables filtering.")
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="tce.ic.dat")
    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    AddPointMassToFile(args.giant_snapshot_file, new_file_name=args.ic_file_name,
                           separation=args.orbital_separation, point_mass=args.point_mass * msol,
                       rlof_factor=args.rlof_factor, giant_radius_rsol=args.giant_radius,
                       padding_cm=args.box_padding_rsun * rsol, restore_eos_radiation=args.restore_eos_radiation,
                       donor_density_cut=args.donor_density_cut, corotation_factor=args.corotation_factor,
                       corotation_inner_radius_rsol=args.corotation_inner_radius,
                       corotation_outer_radius_rsol=args.corotation_outer_radius,
                       corotation_density_cut=args.corotation_density_cut,
                       corotation_profile=args.corotation_profile,
                       corotation_spin_mode=args.corotation_spin_mode)
