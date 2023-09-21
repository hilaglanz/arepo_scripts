import argparse
import sys
import os

from loadmodules import *
from BinariesICs import *


def initialized_new_data(npart):
    data = {}
    data['pos'] = np.zeros((npart, 3))
    data['vel'] = np.zeros((npart, 3))
    data['bfld'] = np.zeros((npart, 3))
    data['mass'] = np.zeros(npart)
    data['u'] = np.zeros(npart)
    data['pass'] = np.zeros((npart, 2))
    data['count'] = npart

    return data


def copy_old_data(snapshot):
    data = initialized_new_data(snapshot.npart + 1)
    data['pos'][:-1] = snapshot.data['pos']
    data['vel'][:-1] = snapshot.data['vel']
    data['mass'][:-1] = snapshot.data['mass']
    data['u'][:-1] = snapshot.data['u']
    if 'bfld' in snapshot.data.keys():
        data['bfld'][:-1] = snapshot.data['bfld']
    if 'pass' in snapshot.data.keys():
        data['pass'][:-1] = snapshot.data['pass']

    return data
def AddPointMassToFile(snapshot_file, new_file_name, loadtypes, point_mass, separation, velocity=None):
    loaded_snap = gadget_readsnap(snapshot_file, loadonlytype=loadtypes)
    data = copy_old_data(loaded_snap)
    data['pos'][-1] = loaded_snap.center + np.array([separation, 0.0, 0.0])
    if velocity is None:
        total_mass = loaded_snap.mass.sum() + point_mass
        velocity = np.array([0.0, (G * total_mass / separation).sqrt(), 0.0])
    data['vel'][-1] = velocity
    data['mass'][-1] = point_mass
    if (separation > 0.01 * loaded_snap.boxsize).any():
        gadget_add_grids(data, loaded_snap.boxsize * [10, 100], 32)
        data['boxsize'] = loaded_snap.boxsize*100
    print("added a point mass of ",point_mass," at ", data['pos'][-1])
    gadget_write_ics(new_file_name, data, format='hdf5', double=True)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--giant_snapshot_file', type=str, help='path to snapshot file containing the giant star', default="")
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types, '
                                                                  'if there is a point mass companion one '
                                                                  'should also load type 1', default=[0,1])
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
    AddPointMassToFile(args.giant_snapshot_file, new_file_name=args.ic_file_name, loadtypes=args.load_types,
                       separation=args.outer_separation * rsol, point_mass=args.point_mass * msol)
