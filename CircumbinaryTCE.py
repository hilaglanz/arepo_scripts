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


def ReplaceInnerBinaryWithPointMass(snapshot_file, new_file_name, obj1_id, obj2_id, remove_to_radius=None):
    snapshot = gadget_readsnapname(snapshot_file)
    obj1_ind = get_obj_index(snapshot, obj1_id)
    obj2_ind = get_obj_index(snapshot, obj2_id)
    num_gas_to_remove = 0
    inds_to_remove = np.array([])
    gas_inds_to_remove = np.array([])
    new_pos = ((snapshot.pos[obj1_ind] * snapshot.mass[obj1_ind]).sum() +
               (snapshot.pos[obj2_ind] * snapshot.mass[obj2_ind]).sum()) / \
              (snapshot.mass[obj1_ind] + snapshot.mass[obj2_ind])

    new_vel = ((snapshot.vel[obj1_ind] * snapshot.mass[obj1_ind]).sum() +
               (snapshot.vel[obj2_ind] * snapshot.mass[obj2_ind]).sum()) / \
              (snapshot.mass[obj1_ind] + snapshot.mass[obj2_ind])
    inds_to_remove = np.array([obj1_ind, obj2_ind])
    if remove_to_radius is not None:
        gas_inds_to_remove = np.where((snapshot.type == 0) &
                                 (((snapshot.pos - new_pos)**2).sum(axis=1)**0.5 < remove_to_radius))[0]
        num_gas_to_remove = len(inds_to_remove)
        if num_gas_to_remove > 0:
            inds_to_remove = np.concatenate((inds_to_remove, gas_inds_to_remove))

    new_data = initialized_new_data(snapshot.npart - 1 - num_gas_to_remove, snapshot.npart - 3 - num_gas_to_remove)
    new_data['boxsize'] = snapshot.boxsize
    new_data['count'] = snapshot.npart - 1 - num_gas_to_remove

    for value in snapshot.data:
        print("updating value: ",value)
        if value in ["boxsize","count"]:
            continue
        if value in ["u","temp","B","rho","cmce",'grar', 's', 'u', 'bfld', 'divb', 'dvba','pres', 'grap',
                     'csnd', 'temp', 'tstp','grav', 'vol'] or 'xnuc' in value:
            if num_gas_to_remove > 0:
                new_data[value] = np.delete(snapshot.data[value], gas_inds_to_remove, axis=0)
            else:
                new_data[value] = snapshot.data[value]
        else:
            if value not in new_data.keys():
                if len(snapshot.data[value].shape) > 1:
                    new_data[value] = np.zeros((new_data['count'],3))
                else:
                    new_data[value] = np.zeros(new_data['count'])
            new_data[value][:-1] = np.delete(snapshot.data[value], inds_to_remove, axis=0)

    new_data['pos'][-1] = new_pos
    new_data['vel'][-1] = new_vel
    new_data['mass'][-1] = snapshot.mass[obj1_ind] + snapshot.mass[obj2_ind]
    new_data['type'][-1] = snapshot.type[obj1_ind]
    new_data['id'][-1] = obj1_id

    if 'soft' in new_data.keys():
        new_soft = snapshot.soft[obj1_ind] + (((snapshot.pos[obj1_ind] - snapshot.pos[obj2_ind])**2).sum()**0.5) / 2
        new_data['soft'][-1] = new_soft
    print("removed ", num_gas_to_remove, " gas cells")
    
    gadget_write_ics(new_file_name, new_data, format='hdf5', double=True)




def AddPointMassToFile(snapshot_file, new_file_name, point_mass, separation, initial_inner_separation=None,
                       velocity=None):
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
    q = inner_binary.mass / point_mass
    current_rlof = (inner_binary.get_radius() / rsol) / roche_distance(q)
    print("current Roche lobe size= ", current_rlof, " Rsun")
    if initial_inner_separation is not None:
        print("calculating minimum distance from stability criteria")
        minimum_a = 2.8 * initial_inner_separation * (1 + q)**(2.0/5)
        print("minimal separation from stability criteria = ", minimum_a, " Rsun, using ",
              1.03*minimum_a)
        separation = 1.03 * minimum_a
    else:
        print("using given separation of ", separation, "Rsun")

    rlof_factor = separation / current_rlof
    triple.add_components_as_binary(inner_binary, tertiary, distance_fraction_rlof=rlof_factor)
    triple.create_ics(model=new_file_name)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--giant_snapshot_file', type=str, help='path to snapshot file containing the giant star', default="")
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types, '
                                                                  'if there is a point mass companion one '
                                                                  'should also load type 1', default=[0,1])
    parser.add_argument('--initial_inner_separation', type=float,
                        help='initial inner binary separation in Rsun before CE', default=None)
    parser.add_argument('--outer_mass', type=float, help='new object mass in msun', default=1)
    parser.add_argument('--outer_separation', type=float, help='initial separation between the binary objects in Rsun',
                        default=1000)
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="tce.ic.dat")
    parser.add_argument('--replace_id1', type=int, help='id of the first obj to remove and replace', default=None)
    parser.add_argument('--replace_id2', type=int, help='id of the second obj to remove and replace', default=None)
    parser.add_argument('--remove_to_radius', type=float, help='radius around the com of the merging point mass to '
                                                               'remove all gas cells', default=None)
    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    if args.replace_id1 is not None and args.replace_id2 is not None:
        ReplaceInnerBinaryWithPointMass(snapshot_file=args.giant_snapshot_file, new_file_name=args.ic_file_name,
                                        obj1_id=args.replace_id1, obj2_id=args.replace_id2,
                                        remove_to_radius=args.remove_to_radius)
        print("replaced 2 point masses with 1")
    else:
        AddPointMassToFile(args.giant_snapshot_file, new_file_name=args.ic_file_name,
                           separation=args.outer_separation, point_mass=args.outer_mass * msol,
                           initial_inner_separation=args.initial_inner_separation)
