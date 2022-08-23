import os.path

import numpy as np
import argparse
from loadmodules import *


class BinariesICs:
    def __init__(self, species_file="species55.txt"):
        print("initializing general binary info")
        self.init_species(species_file)
        self.init_object_properties()
        self.init_orbital_parameters()

    def init_species(self, species_file):
        self.sp = loaders.load_species(species_file)
        self.iC12 = self.sp['names'].index('c12')
        self.iO16 = self.sp['names'].index('o16')
        print("species file set")

    def init_object_properties(self):
        self.m1 = self.snapshot1.mass[self.i1].sum()
        self.m2 = self.snapshot2.mass[self.i2].sum()
        self.total_mass = self.m1 + self.m2
        self.r1 = self.snapshot1.r()[self.i1].max()
        self.r2 = self.snapshot2.r()[self.i2].max()
        self.npart1 = size(self.i1)
        self.npart2 = size(self.i2)
        self.npart = self.npart1 + self.npart2
        print("objects properties set")

    def init_orbital_parameters(self):
        self.pos1 = (self.snapshot1.mass[self.i1] * self.self.snapshot1.pos[self.i1]).sum(axis=0) / self.m1
        self.pos2 = (self.self.snapshot2.mass[self.i2] * self.self.snapshot2.pos[self.i2]).sum(axis=0) / self.m2
        self.loaded_separation = self.pos2 - self.pos1
        self.center_of_mass = (self.m1 * self.pos1 + self.m2 * self.pos2) / self.total_mass

        self.v1 = (self.snapshot1.mass[self.i1] * self.snapshot1.v[self.i1]).sum(axis=0) / self.m1
        self.v2 = (self.snapshot2.mass[self.i2] * self.snapshot2.v[self.i2]).sum(axis=0) / self.m2
        self.velocity_difference = self.v2 - self.v1

        self.angular_momentum = np.cross(self.loaded_separation, self.velocity_difference)

        self.eccentricity = np.cross(self.velocity_difference, self.angular_momentum) / G * (self.m1 + self.m2) - \
                            self.loaded_separation / np.dot(self.loaded_separation, self.loaded_separation) ** 0.5

        self.semimajor = (self.angular_momentum ** 2) / (G * (self.m1 + self.m2) * (1 - self.eccentricity ** 2))

        print("current orbital parameters: a=", self.semimajor, " e=", self.eccentricity, " r=", self.loaded_separation,
              " l=", self.angular_momentum, " v=", self.velocity_difference)

    def initialized_new_data(self):
        self.data = {}
        self.data['pos'] = np.zeros((self.npart, 3))
        self.data['vel'] = np.zeros((self.npart, 3))
        self.data['bfld'] = np.zeros((self.npart, 3))
        self.data['mass'] = np.zeros(self.npart)
        self.data['u'] = np.zeros(self.npart)
        self.data['xnuc'] = np.zeros((self.npart, self.sp['count']))
        self.data['pass'] = np.zeros((self.npart, 2))
        self.data['count'] = self.npart

    def copy_old_data(self):
        self.data['mass'][:self.npart1] = self.snapshot1.mass[self.i1]
        self.data['u'][:self.npart1] = self.snapshot1.data['u'][self.i1]
        self.data['xnuc'][:self.npart1, :] = self.snapshot1.data['xnuc'][self.i1, :]
        self.data['pass'][:self.npart1, 0] = 1.0

        self.data['mass'][self.npart1:] = self.snapshot2.mass[self.i2]
        self.data['u'][self.npart1:] = self.snapshot2.data['u'][self.i2]
        self.data['xnuc'][self.npart1:, :] = self.snapshot2.data['xnuc'][self.i2, :]
        self.data['pass'][self.npart1:, 1] = 1.0

        self.data['pos'][:npart1, :] = self.snapshot1.pos[i1, :]
        self.data['vel'][:npart1, :] = self.snapshot1.vel[i1, :]

        self.data['pos'][npart1:, :] = self.snapshot2.pos[i2, :]
        self.data['vel'][npart1:, :] = self.snapshot2.vel[i2, :]

        if 'bfld' in self.snapshot1.data.keys():
            self.data['bfld'][:self.npart1, :] = self.snapshot1.data['bfld'][self.i1, :]
            self.data['bfld'][self.npart1:, :] = self.snapshot2.data['bfld'][self.i2, :]

        self.data['boxsize'] = 1e10

    def find_new_borders(self):
        return 4.0 * np.sqrt(self.data['pos'][:, 0] ** 2, self.data['pos'][:, 1] ** 2,
                             self.data['pos'][:, 2] ** 2).max()

    def add_grids_and_save_ic(self, ic_file_name):
        boxsize = self.data['boxsize']
        print("using inner boxsize= ", boxsize)
        xnuc = np.zeros(self.sp['count'])
        xnuc[self.iC12] = 0.5
        xnuc[self.iO16] = 0.5
        gadget_add_grids(self.data, [boxsize, 10 * boxsize, 100 * boxsize], 32, xnuc=xnuc)
        gadget_write_ics(ic_file_name, self.data, double=True, format="hdf5")
        print("ic file saved to ", ic_file_name)

    def create_ic_merger_keplerian(self,ic_file_name="bin.dat.ic", T=120, rhocut=1, relative_to_RL=True, factor=2):
        a0 = self.calculate_RL()
        w0 = sqrt(G * m / a0 ** 3)
        T0 = 2.0 * pi / w0
        print('Egg:', T0, a0, w0)

        w = 2.0 * pi / T
        d = (G * m / (w * w)) ** (1. / 3.)
        print('Using:', T, d, w)

        self.initialized_new_data()
        self.copy_old_data()

        self.relative_x = d
        self.relative_y = 0
        self.relative_vx = 0
        self.relative_vy = 0
        self.create_new_position_array()
        self.place_objects_at_new_pos()

        print("changing relative velocities according to new positions")
        self.relative_vx = w * self.data['pos'][:,1]
        self.relative_vy = -w * self.data['pos'][:,0]
        self.create_new_velocity_array()
        self.change_objects_velocity()

        self.add_magnetic_field()

        self.add_grids_and_save_ic(ic_file_name)

    def create_ic_at_apastron(self, semimajor, eccentricity=0):
        raise Exception("not implemented")

    def create_ic_collision(self, impact_parameter, velocity=16e8, ic_file_name="bin.dat.ic"):
        self.relative_y = impact_parameter
        self.relative_x = self.calculate_RL()

        self.relative_vx = velocity
        self.relative_vy = 0.0

        self.create_new_position_velocity_arrays()
        self.place_objects_at_new_pos()
        self.change_objects_velocity()
        self.add_grids_and_save_ic(ic_file_name)

    def create_ic_for_next_interaction(self, ic_file_name="bin.dat.ic", relative_to_RL=True, factor=2, dist=None):
        self.initialized_new_data()
        self.copy_old_data()
        if relative_to_RL:
            dist = factor * self.calculate_RL() # factor*RL
        elif dist is None:
            print("cannot find desired position")
            return
        print("starting new system at a separation= ", dist)

        self.find_keplerian_orbit_at_distance(dist)
        self.place_objects_at_new_pos()
        self.change_objects_velocity()
        print("objects new positions and velocities set adding grid and saving")
        self.add_grids_and_save_ic(ic_file_name)

    def GetRelativeVelocityAtAngel(self, orbitalPhase):
        coefficient = (G * self.total_mass / (self.semimajor * (1 - self.eccentricity ** 2))).sqrt()
        vTangelntail = coefficient * (self.eccentricity + math.cos(orbitalPhase))
        vRadial = coefficient * math.sin(orbitalPhase)
        return [vRadial, vTangelntail]

    def find_keplerian_orbit_at_distance(self, dist):
        orbitalPhaseGoal = -1.0 * np.arccos(self.semimajor * (1 - self.eccentricity ** 2) /
                                            (self.eccentricity * dist) - 1.0 / self.eccentricity)
        print("using orbital phase=", orbitalPhaseGoal)
        self.relative_x = dist * math.sin((pi - orbitalPhaseGoal) % (2 * pi))
        self.relative_y = -dist * np.cos(orbitalPhaseGoal)
        print("relative_x= ", self.relative_x, " relative_y= ",self.relative_y)

        [self.relative_vy, self.relative_vx] = self.GetRelativeVelocityAtAngel(orbitalPhaseGoal)

        self.create_new_position_velocity_arrays()

    def create_new_position_velocity_arrays(self):
        self.create_new_position_array()
        self.create_new_velocity_array()

    def create_new_position_array(self):
        self.new_x1 = - self.m2 * relative_x / self.total_mass
        self.new_x2 = self.m1 * relative_x / self.total_mass

        self.new_y1 = - self.m2 * relative_y / self.total_mass
        self.new_y2 = self.m1 * relative_y / self.total_mass

        self.new_pos1 = np.array([self.new_x1, self.new_y1, self.pos1[2]])
        self.new_pos2 = np.array([self.new_x2, self.new_y2, self.pos2[2]])

    def create_new_velocity_array(self):
        self.new_vx1 = - self.m2 * relative_vx / self.total_mass
        self.new_vx2 = self.m1 * relative_vx / self.total_mass

        self.new_vy1 = - self.m2 * relative_vy / self.total_mass
        self.new_vy2 = self.m1 * relative_vy / self.total_mass

        self.new_v1 = np.array([self.new_vx1, self.new_vy1, self.v1[2]])
        self.new_v2 = np.array([self.new_vx2, self.new_vy2, self.v2[2]])

    def place_objects_at_new_pos(self):
        self.data['pos'][:npart1, :] = self.data['pos'][:self.npart1, :] + self.new_pos1 - self.pos1
        self.data['pos'][npart1:, :] = self.data['pos'][self.npart1:, :] + self.new_pos2 - self.pos2

        self.data['boxsize'] = max(self.find_new_borders(), 1e10)
        self.data['pos'] += 0.5 * self.data['boxsize']

    def change_objects_velocity(self):
        self.data['vel'][:npart1, :] = self.data['vel'][:self.npart1, :] + self.new_v1 - self.v1
        self.data['vel'][npart1:, :] = self.data['vel'][self.npart1:, :] + self.new_v1 - self.v2

    def add_magnetic_field(self):
        mm = np.array([0., 0., 1e3 * 1e9 ** 3 / 2.])  # 1e3 G at 1e9 cm
        c1 = (self.data['pos'] * self.data['mass'][:,None] * self.data['pass'][:,0][:,None]).sum( axis=0 ) / self.m1
        r1 = np.maximum(np.sqrt(((self.data['pos'] - c1[None, :]) ** 2).sum(axis=1)), 3e7).astype('float64')
        i, = np.where(r1 < 1e10)
        rad1 = self.data['pos'][i, :] - c1[None, :]
        self.data['bfld'][i, :] = 3. * rad1 * (mm[None, :] * rad1).sum(axis=1)[:, None] / (r1[i] ** 5)[:, None] - \
                                  mm[None,:] / (r1[i] ** 3)[:, None]

        c2 = (self.data['pos'] * self.data['mass'][:, None] * self.data['pass'][:, 1][:, None]).sum(axis=0) / self.m2
        r2 = np.maximum(np.sqrt(((data['pos'] - c2[None, :]) ** 2).sum(axis=1)), 3e7).astype('float64')
        i, = np.where(r2 < 1e10)
        rad2 = data['pos'][i, :] - c2[None, :]
        self.data['bfld'][i, :] += 3. * rad2 * (mm[None, :] * rad2).sum(axis=1)[:, None] / (r2[i] ** 5)[:, None] - \
                              mm[None, :] / (r2[i] ** 3)[:, None]

    def calculate_RL(self):
        q = self.m2 / self.m1
        RL2 = self.r2 *  (0.6 * q ** (2. / 3.) + log(1 + q ** (1. / 3.))) / (0.49 * q ** (2. / 3.))

        q = self.m1 / self.m2
        RL1 = self.r1 *  (0.6 * q ** (2. / 3.) + log(1 + q ** (1. / 3.))) / (0.49 * q ** (2. / 3.))

        return max(RL1, RL2)

class BinariesLoader(BinariesICs):
    def __init__(self, snapshot_file, conditional_axis=0, rhocut=1, species_file="species55.txt"):
        print("initializeing binaries from a single snapshot")
        self.binary_snapshot = gadget_readsnapname(snapshot_file, hdf5=True, loadonlytype=[0])
        if conditional_axis is not None:
            if rhocut is not None:
                self.get_objects_position_density(conditional_axis,rhocut)
            else:
                self.get_objects_position(conditional_axis)
        else:
            print("using passive scalars to distinguish")
            if rhocut is not None:
                self.get_objects_density(rhocut)
            else:
                self.get_objects()
        self.snapshot1 = self.binary_snapshot
        self.snapshot2 = self.binary_snapshot

        super().__init__(species_file)

    def get_objects_position_density(self, conditional_axis, rhocut):
        self.i1 = np.where((self.binary_snapshot.pos[:,conditional_axis] - self.binary_snapshot.center[conditional_axis] > 0) &
                        (self.binary_snapshot.rho > rhocut))
        self.i2 = np.where((self.binary_snapshot.pos[:,conditional_axis] - self.binary_snapshot.center[conditional_axis] < 0) &
                        (self.binary_snapshot.rho > rhocut))

    def get_objects_position(self, conditional_axis):
        self.i1 = np.where((self.binary_snapshot.pos[:,conditional_axis] - self.binary_snapshot.center[conditional_axis] > 0))
        self.i2 = np.where((self.binary_snapshot.pos[:,conditional_axis] - self.binary_snapshot.center[conditional_axis] < 0))

    def get_objects_density(self, rhocut):
        self.i1 = np.where((self.binary_snapshot.pass00 == 1) & (self.binary_snapshot.rho > rhocut))
        self.i2 = np.where((self.binary_snapshot.pass01 == 1) & (self.binary_snapshot.rho > rhocut))

    def get_objects(self):
        self.i1 = np.where((self.binary_snapshot.pass00 == 1))
        self.i2 = np.where((self.binary_snapshot.pass01 == 1))


class SeparateBinariesLoader(BinariesICs):
    def __init__(self, snapshot_file1, snapshot_file2, rhocut=1, species_file="species55.txt"):
        self.snapshot1 = gadget_readsnapname(snapshot_file1, hdf5=True, loadonlytype=[0])
        self.snapshot2 = gadget_readsnapname(snapshot_file2, hdf5=True, loadonlytype=[0])
        self.get_objects_density(rhocut)
        super().__init__(species_file)

    def get_objects_density(self, rhocut):
        self.i1 = np.where((self.snapshot.rho > rhocut))
        self.i2 = np.where((self.snapshot.rho > rhocut))


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--snapshot_file', type=str,  help='path to snapshot file', default="")
    parser.add_argument('--snapshot_file2', type=str,  help='path to secind snapshot file', default="")
    parser.add_argument('--species_file', type=str,  help='path to species file', default="species55.txt")
    parser.add_argument('--conditional_axis', type=int,  help='axis of motion', default=0)
    parser.add_argument('--rhocut', type=float,  help='lower cutoff of density', default=1)
    parser.add_argument('--period', type=float,  help='orbital period in seconds for mergers', default=0)
    parser.add_argument('--impact_parameter_rhocut', type=float,
                        help='calculate impact parameter according to this density cutoff', default=0)
    parser.add_argument('--relative_velocity', type=float, help='', default=16e8)
    parser.add_argument('--relative_to_RL', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='is the distance should be relative to RL size?',
                        default=True)
    parser.add_argument('--RL_factor', type=float, help='if relative to RL, by what factor?', default=2)
    parser.add_argument('--find_next_interaction', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='is the distance should be relative to RL size?',
                        default=False)
    parser.add_argument('--separation', type=float,  help='initial separation between the binary objects', default=None)
    parser.add_argument('--ic_file_name', type=str,  help='path to save the ic file', default= "bin.ic.dat")
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    if args.snapshot_file2 is None:
        binary = BinariesLoader(args.snapshot_file,conditional_axis=args.conditional_axis,
                                rhocut=args.rhocut, species_file=args.species_file)
    else:
        binary = SeparateBinariesLoader(args.snapshot_file, args.snapshot_file2, args.rhocut, args.species_file)

    if args.period > 0:
        binary.create_ic_merger_keplerian(ic_file_name=args.ic_file_name, T= args.period, rhocut= args.rhocut,
                                          relative_to_RL=args.relative_to_RL, factor=args.RL_factor)
    elif args.impact_parameter_rhocut > 0:
        i1 = np.where(binary.snapshot1.rho > args.impact_parameter_rhocut)
        i2 = np.where(binary.snapshot2.rho > args.impact_parameter_rhocut)
        b = binary.snapshot1.r()[i1].max() + binary.snapshot2.r()[i2].max()
        binary.create_ic_collision(b, args.relative_velocity, args.ic_file_name)
    elif args.find_next_interaction:
        binary.create_ic_for_next_interaction(args.ic_file_name, args.relative_to_RL, args.RL_factor, args.separation)









