import argparse
import sys
import os

from loadmodules import *
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
    data = initialized_new_data(snapshot.npart + 1, snapshot.nparticles[0])
    data['pos'][:-1] = snapshot.data['pos']
    data['vel'][:-1] = snapshot.data['vel']
    data['mass'][:-1] = snapshot.data['mass']
    data['type'][:-1] = snapshot.data['type']
    data['u'] = snapshot.data['u']
    if 'bfld' in snapshot.data.keys():
        data['bfld'] = snapshot.data['bfld']
    if 'pass' in snapshot.data.keys():
        data['pass'] = snapshot.data['pass']

    return data
def AddPointMassToFile(snapshot_file, new_file_name, loadtypes, point_mass, separation, velocity=None):
    loaded_snap = gadget_readsnapname(snapshot_file, loadonlytype=loadtypes)
    data = copy_old_data(loaded_snap)
    data['pos'][-1] = loaded_snap.center + np.array([separation, 0.0, 0.0])
    if velocity is None:
        total_mass = loaded_snap.mass.sum() + point_mass
        velocity = np.array([0.0, (G * total_mass / separation)**0.5, 0.0])
    data['vel'][-1] = velocity
    data['mass'][-1] = point_mass
    data['type'][-1] = 1     #dm particle
    gadget_add_grid(data, loaded_snap.boxsize * 10, 32)
    if (separation > 0.01 * loaded_snap.boxsize).any():
        print("expanding the box")
        gadget_add_grid(data, loaded_snap.boxsize * 100, 32)
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
                       separation=args.outer_separation * rsol, point_mass=args.outer_mass * msol)
