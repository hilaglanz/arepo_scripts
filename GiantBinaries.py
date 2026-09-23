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

def AddPointMassToFile(snapshot_file, new_file_name, point_mass, separation, rlof_factor=1.0, giant_radius_rsol=None, padding_cm=0):
    snapshot=gadget_readsnapname(snapshot_file)
    new_size = snapshot.boxsize
    
    giant = SnapshotComponent(data=snapshot.data, boxsize=snapshot.boxsize, radius=giant_radius_rsol)
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
    binary = MultipleSystem(newsize=new_size,
                            reset_dm_ids=True, ndir=32, grid_xnuc=snapshot.data['xnuc'][0],
                            grid_rho=min([snapshot.rho.min(), 1e-20]),
                            grid_u=min([snapshot.data['u'].min(), 1e10]))
    binary.add_components_as_binary(giant, companion, distance_fraction_rlof=rlof_factor, corotating_at_rlof=False, corotation_factor=0.0, e=0.0)
    binary.newsize = enclosing_boxsize((giant, companion), snapshot.boxsize, padding_cm)
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
    parser.add_argument("--box_padding_rsun", type=float, help="Clearance beyond current retained particle centers, in solar radii", default=0.0)
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
                       rlof_factor=args.rlof_factor, giant_radius_rsol=args.giant_radius, padding_cm=args.box_padding_rsun * rsol)
