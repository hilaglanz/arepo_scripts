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

def AddPointMassToFile(snapshot_file, new_file_name, point_mass, separation, rlof_factor=1.0):
    snapshot=gadget_readsnapname(snapshot_file)
    new_size = snapshot.boxsize
    if separation > new_size/100:
        new_size *= 100
    binary = MultipleSystem(newsize=new_size,
                            reset_dm_ids=True, ndir=32, grid_xnuc=snapshot.data['xnuc'][0],
                            grid_rho=min([snapshot.rho.min(), 1e-20]),
                            grid_u=min([snapshot.data['u'].min(), 1e10]))
    giant = SnapshotComponent.from_snapshot_name(snapshot_file)
    companion = PointMassComponent(mass=point_mass)
    companion.data['type'] = [5]
    q = giant.mass / point_mass
    current_rlof = (giant.get_radius() / rsol) / roche_distance(q)
    print("current Roche lobe size= ", current_rlof, " Rsun")
    if separation is None:
        print("calculating separation from Roche Lobe")
        separation = current_rlof * rlof_factor
    else:
        print("using given separation of ", separation, "Rsun")
        rlof_factor = separation / current_rlof

    print("Roche factor = ", rlof_factor)
    binary.add_components_as_binary(giant, companion, distance_fraction_rlof=rlof_factor)
    binary.create_ics(model=new_file_name)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--giant_snapshot_file', type=str, help='path to snapshot file containing the giant star', default="")
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types, '
                                                                  'if there is a point mass companion one '
                                                                  'should also load type 1 or 5', default=[0,1])
    parser.add_argument('--orbital_separation', type=float,
                        help='initial binary separation in Rsun', default=None)
    parser.add_argument('--point_mass', type=float, help='new object mass in msun', default=1)
    parser.add_argument('--rlof_factor', type=float, help='if relative to RL, by what factor?', default=2)
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="tce.ic.dat")
    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    AddPointMassToFile(args.giant_snapshot_file, new_file_name=args.ic_file_name,
                           separation=args.orbital_separation, point_mass=args.point_mass * msol, rlof_factor=args.rlof_factor)
