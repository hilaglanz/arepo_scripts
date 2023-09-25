import argparse
import sys
import os

from loadmodules import *
from stellar_ics.tools import *
from stellar_ics.multiple_star import MultipleSystem, SnapshotComponent, PointMassComponent
from BinariesICs import *


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

def AddPointMassToFile(snapshot_file, new_file_name, point_mass, separation,
                       initial_rg_radius=None, initial_inner_separation=None, velocity=None):
    snapshot=gadget_readsnapname(snapshot_file)
    new_size = snapshot.boxsize
    if separation > new_size/100:
        new_size *= 100
    triple = MultipleSystem(newsize=new_size,
                            reset_dm_ids=True, ndir=32, grid_xnuc=snapshot.data['xnuc'][0],
                            grid_rho=min([snapshot.rho.min(), 1e-20]),
                            grid_u=min([snapshot.data['u'].min(), 1e10]))
    inner_binary = SnapshotComponent.from_snapshot_name(snapshot_file)
    tertiary = PointMassComponent(mass=point_mass)
    rlof_factor = 1.0
    current_rlof = inner_binary.get_radius() * roche_distance(point_mass / inner_binary.mass)
    print("current Roche lobe size= ", current_rlof / rsol, " Rsun")
    if initial_rg_radius is not None and initial_inner_separation is not None:
            print("calculating minimum distance from stability criteria")
            minimum_a = 2.8 * initial_inner_separation * (1+point_mass/inner_binary.mass)**(2.0/5)
            print("minimal separation from stability criteria = ", minimum_a, " Rsun, using ",
                  1.1*minimum_a)
            separation = 1.1 * minimum_a * rsol
    else:
        print("using given separation of ", separation/rsol, "Rsun")

    rlof_factor = separation / current_rlof
    triple.add_components_as_binary(inner_binary, tertiary, distance_fraction_rlof=separation/rlof_factor)
    triple.create_ics(model=new_file_name)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--giant_snapshot_file', type=str, help='path to snapshot file containing the giant star', default="")
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types, '
                                                                  'if there is a point mass companion one '
                                                                  'should also load type 1', default=[0,1])
    parser.add_argument('--giant_initial_radius', type=float,
                        help='initial radius of the giant in Rsun before CE', default=None)
    parser.add_argument('--initial_inner_separation', type=float,
                        help='initial inner binary separation in Rsun before CE', default=None)
    parser.add_argument('--outer_mass', type=float, help='new object mass in msun', default=1)
    parser.add_argument('--outer_separation', type=float, help='initial separation between the binary objects in Rsun',
                        default=1000)
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="tce.ic.dat")
    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    AddPointMassToFile(args.giant_snapshot_file, new_file_name=args.ic_file_name,
                       separation=args.outer_separation * rsol, point_mass=args.outer_mass * msol,
                       initial_rg_radius=args.giant_initial_radius,
                       initial_inner_separation=args.initial_inner_separation)
