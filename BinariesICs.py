import os.path

import numpy as np
import argparse
from loadmodules import *


class SingleObject:
    def __init__(self, snapshot, rhocut=1e0, relevant_ind=[]):
        self.pos = None
        self.v = None
        self.r = None
        self.snapshot = snapshot
        if len(relevant_ind) > 0:
            self.i = relevant_ind
        else:
            self.get_objects_density(rhocut)
        self.m = self.snapshot.mass[self.i].astype('float64').sum()
        self.max_distance = self.snapshot.r()[self.i].max()  # maximum distance from the center of the box
        self.npart = size(self.i)
        self.init_orbital_parameters()

    def get_objects_density(self, rhocut=1.0):
        self.i, = np.where((self.snapshot.rho > rhocut))

    def init_orbital_parameters(self):
        print("setting orbital parameters")
        self.pos = (self.snapshot.pos[self.i, :].astype('f8') * self.snapshot.mass[self.i][:, None]).sum(
            axis=0) / self.m
        print("pos: ", self.pos)
        self.r = self.snapshot.r(center=self.pos)[self.i].max()  # radius of particles
        self.v = (self.snapshot.vel[self.i, :].astype('f8') * self.snapshot.mass[self.i][:, None]).sum(axis=0) / self.m


class BinariesICs:
    snapshot1 = None
    snapshot2 = None
    i1 = []
    i2 = []

    def __init__(self, species_file="species55.txt", rhocut=1.0):
        self.relative_z = None
        self.new_z1 = None
        self.new_z2 = None
        self.new_v1 = None
        self.new_v2 = None
        self.new_vx1 = None
        self.new_vx2 = None
        self.new_vy1 = None
        self.new_vy2 = None
        self.new_x1 = None
        self.new_x2 = None
        self.new_y1 = None
        self.new_y2 = None
        self.new_pos1 = None
        self.new_pos2 = None
        self.relative_x = None
        self.relative_y = None
        self.relative_vx = None
        self.relative_vy = None
        self.npart = None
        self.npart2 = None
        self.npart1 = None
        self.total_mass = None
        self.m2 = None
        self.m1 = None
        self.semimajor = None
        self.eccentricity = None
        self.eccentricity_vector = None
        self.angular_momentum = None
        self.velocity_difference = None
        self.v2 = None
        self.v1 = None
        self.center_of_mass = None
        self.loaded_separation = None
        self.r2 = None
        self.r1 = None
        self.pos2 = None
        self.pos1 = None
        if self.snapshot1 is None or self.snapshot2 is None:
            print("snapshot weren't initialized exiting")
            return

        self.obj1 = SingleObject(self.snapshot1, rhocut, self.i1)
        self.obj2 = SingleObject(self.snapshot2, rhocut, self.i2)

        print("initializing general binary info")
        if species_file is not None:
            self.init_species(species_file)
        self.init_object_properties()
        self.init_orbital_parameters()
        self.initialized_new_data()
        self.copy_old_data()

    def init_species(self, species_file):
        self.sp = loaders.load_species(species_file)
        self.iC12 = self.sp['names'].index('c12')
        self.iO16 = self.sp['names'].index('o16')
        print("species file set")

    def init_object_properties(self):
        self.m1 = self.obj1.m
        self.m2 = self.obj2.m
        self.total_mass = self.m1 + self.m2
        print("masses: ", self.m1, " , ", self.m2)
        self.npart1 = self.obj1.npart
        self.npart2 = self.obj2.npart
        self.npart = self.npart1 + self.npart2
        print("objects properties set")

    def init_orbital_parameters(self):
        print("setting orbital parameters")
        self.pos1 = self.obj1.pos
        self.pos2 = self.obj2.pos
        print("pos1: ", self.pos1, "pos2: ", self.pos2)
        self.r1 = self.obj1.r  # radius of particles1
        self.r2 = self.obj2.r  # radius of particles2
        self.loaded_separation = self.pos2 - self.pos1
        print("loaded separation: ", self.loaded_separation)
        self.center_of_mass = (self.m1 * self.pos1 + self.m2 * self.pos2) / self.total_mass
        print("binary center of mass: ", self.center_of_mass)
        self.v1 = self.obj1.v
        self.v2 = self.obj2.v
        self.velocity_difference = self.v2 - self.v1
        print("relative velocity: ", self.velocity_difference)
        self.angular_momentum = np.cross(self.loaded_separation, self.velocity_difference)
        print("angular momentum: ", self.angular_momentum)
        print("G= ", G)
        self.eccentricity_vector = np.cross(self.velocity_difference, self.angular_momentum) / (
                    G * (self.m1 + self.m2)) - \
                                   self.loaded_separation / np.sqrt((self.loaded_separation ** 2).sum())
        print("eccentricity vector: ", self.eccentricity_vector)
        self.eccentricity = np.sqrt((self.eccentricity_vector ** 2).sum())
        print("eccentricity: ", self.eccentricity)
        self.semimajor = ((self.angular_momentum ** 2).sum()) / (G * (self.m1 + self.m2) * (1 - self.eccentricity ** 2))

        print("current orbital parameters: a=", self.semimajor, " r=", self.loaded_separation)

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

    def copy_old_data_to_objects(self, i_begin, i_end, obj, passive_scalar):
        self.data['mass'][i_begin: i_end] = obj.snapshot.mass[obj.i]
        self.data['u'][i_begin: i_end] = obj.snapshot.data['u'][obj.i]
        self.data['xnuc'][i_begin: i_end, :] = obj.snapshot.data['xnuc'][obj.i, :]
        self.data['pass'][i_begin: i_end, passive_scalar] = 1.0
        self.data['pos'][i_begin: i_end, :] = obj.snapshot.pos[obj.i, :]
        self.data['vel'][i_begin: i_end, :] = obj.snapshot.vel[obj.i, :]
        if 'bfld' in obj.snapshot.data.keys():
            print("using bfld from snapshot")
            self.data['bfld'][i_begin: i_end, :] = obj.snapshot.data['bfld'][obj.i, :]

    def copy_old_data(self):
        self.copy_old_data_to_objects(0, self.npart1, self.obj1, 0)
        self.copy_old_data_to_objects(self.npart1, None, self.obj2, 1)

        self.data['boxsize'] = max(self.snapshot1.boxsize, self.snapshot2.boxsize)

    def find_new_borders(self):
        return 1.1 * np.sqrt(self.data['pos'][:, 0] ** 2 + self.data['pos'][:, 1] ** 2 +
                             self.data['pos'][:, 2] ** 2).max()

    def add_grids_and_save_ic(self, ic_file_name):
        xnuc = np.zeros(self.sp['count'])
        xnuc[self.iC12] = 0.5
        xnuc[self.iO16] = 0.5
        print("current boxsize= ", self.data["boxsize"])
        self.data['boxsize'] = max(self.data['boxsize'], self.find_new_borders())
        self.data['pos'] += 0.5 * self.data['boxsize']
        boxsize = self.data['boxsize']
        print("using inner boxsize= ", boxsize)

        gadget_add_grids(self.data, [boxsize, 10 * boxsize, 100 * boxsize], 32, xnuc=xnuc)
        gadget_write_ics(ic_file_name, self.data, double=True, format="hdf5")
        print("ic file saved to ", ic_file_name)

    def create_ic_merger_keplerian(self, ic_file_name="bin.dat.ic", initial_period=120, rhocut=1, relative_to_RL=True,
                                   factor=2):
        a0 = self.calculate_RL()
        w0 = sqrt(G * self.total_mass / a0 ** 3)
        T0 = 2.0 * pi / w0
        print('orbital parameters at Roche Lobe filling (T0,a0,w0): ', T0, a0, w0)

        if relative_to_RL:
            initial_distance = a0 * factor
            orbital_vel = sqrt(G * self.total_mass / initial_distance ** 3)
            initial_period = 2.0 * pi / orbital_vel
        else:
            orbital_vel = 2.0 * pi / initial_period
            initial_distance = (G * self.total_mass / (orbital_vel * orbital_vel)) ** (1. / 3.)
        print('Using orbital parameters (T,a,w): ', initial_period, initial_distance, orbital_vel)

        '''self.initialized_new_data()
        self.copy_old_data()'''  # TODO:remove

        self.relative_x = initial_distance
        self.relative_y = 0
        self.relative_z = 0
        self.relative_vx = 0
        self.relative_vy = 0
        self.create_new_position_array()
        self.place_objects_at_new_pos()

        print("changing relative velocities according to new positions")
        self.relative_vx = orbital_vel * self.data['pos'][:, 1]
        self.relative_vy = -orbital_vel * self.data['pos'][:, 0]

        self.create_new_velocity_array()
        self.change_objects_velocity()

        self.data['vel'][:, 0] += self.relative_vx
        self.data['vel'][:, 1] -= self.relative_vy

        #TODO: should update new_v1, new_v2?
        self.add_magnetic_field()

        self.add_grids_and_save_ic(ic_file_name)

    def create_ic_at_apastron(self, semimajor, eccentricity=0):
        raise Exception("not implemented")

    def create_ic_collision(self, impact_parameter, ic_file_name="bin.dat.ic", velocity=None, separation=None,
                            relative_to_RL=1):
        self.relative_y = impact_parameter
        if separation is None:
            if relative_to_RL is None:
                relative_to_RL = 1
            self.relative_x = relative_to_RL * self.calculate_RL()
        else:
            self.relative_x = separation
        print("relative positions= ", self.relative_x, self.relative_y)

        velocity = self.calculate_escape_velocity(velocity, self.relative_x)
        print("using relative velocity= ", velocity)
        self.relative_vx = velocity
        self.relative_vy = 0.0

        self.create_new_position_velocity_arrays()
        self.place_objects_at_new_pos()
        self.change_objects_velocity()
        self.add_magnetic_field(object_R=1e9)
        self.add_grids_and_save_ic(ic_file_name)

    def calculate_escape_velocity(self, velocity, distance):
        if velocity is None:
            print("computing escape velocity at separation ", distance)
            velocity = ((2 * G * (self.total_mass) / distance) ** 0.5)
        return velocity

    def create_ic_for_next_interaction(self, ic_file_name="bin.dat.ic", relative_to_RL=True, factor=2, dist=None):
        '''self.initialized_new_data()
                self.copy_old_data()'''  # TODO:remove
        if relative_to_RL:
            dist = factor * self.calculate_RL()  # factor*RL
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
        print("relative_x= ", self.relative_x, " relative_y= ", self.relative_y)

        [self.relative_vy, self.relative_vx] = self.GetRelativeVelocityAtAngel(orbitalPhaseGoal)

        self.create_new_position_velocity_arrays()

    def create_new_position_velocity_arrays(self):
        self.create_new_position_array()
        self.create_new_velocity_array()

    def create_new_position_array(self):
        self.new_x1 = - self.m2 * self.relative_x / self.total_mass
        self.new_x2 = self.m1 * self.relative_x / self.total_mass

        self.new_y1 = - self.m2 * self.relative_y / self.total_mass
        self.new_y2 = self.m1 * self.relative_y / self.total_mass

        self.new_z1 = - self.m2 * self.relative_z / self.total_mass
        self.new_z2 = self.m1 * self.relative_z / self.total_mass

        self.new_pos1 = np.array([self.new_x1, self.new_y1, self.new_z1])
        self.new_pos2 = np.array([self.new_x2, self.new_y2, self.new_z2])

    def create_new_velocity_array(self):
        if len(self.relative_vx) > 1:
            print("received a vector for relative velocity in x, using just 0 for relative velocity")
            self.new_v1 = np.array([0, 0, 0])
            self.new_v2 = np.array([0, 0, 0])
            return

        print("using relative velocity ", self.relative_vx, self.relative_vy, " for new relative velocity")
        self.new_vx1 = self.m2 * self.relative_vx / self.total_mass
        self.new_vx2 = -self.m1 * self.relative_vx / self.total_mass

        self.new_vy1 = self.m2 * self.relative_vy / self.total_mass
        self.new_vy2 = -self.m1 * self.relative_vy / self.total_mass

        self.new_v1 = np.array([self.new_vx1, self.new_vy1, 0])
        self.new_v2 = np.array([self.new_vx2, self.new_vy2, 0])

    def change_com_vector(self, value, i_begin, i_end, new_vector, old_vector):
        return self.data[value][i_begin:i_end, :] + new_vector[None, :] - old_vector[None, :] #TODO: check this for mergers

    def place_objects_at_new_pos(self):
        print("changing positions from ", self.pos1, self.pos2, " to ", self.new_pos1, self.new_pos2)
        self.data['pos'][:self.npart1, :] = self.change_com_vector('pos', 0, self.npart1, self.new_pos1, self.pos1)
        self.data['pos'][self.npart1:, :] = self.change_com_vector('pos', self.npart1, None, self.new_pos2, self.pos2)

    def change_objects_velocity(self):
        self.data['vel'][:self.npart1, :] = self.change_com_vector('vel', 0, self.npart1, self.new_v1, self.v1)
        self.data['vel'][self.npart1:, :] = self.change_com_vector('vel', self.npart1, None, self.new_v2, self.v2)

    def calculate_magnetic_field_around_pos(self, pos, seed_B=1e3, object_R=1e9, replace=True):
        mm = np.array([0., 0., seed_B * object_R ** 3 / 2.])  # 1e3 G at 1e9 cm
        relative_separation = np.maximum(np.sqrt(((self.data['pos'] - pos[None, :]) ** 2).sum(axis=1)), 3e7).astype('float64')
        i_box, = np.where(relative_separation < self.data['boxsize'])
        relative_pos = self.data['pos'][i_box, :] - pos[None, :]
        bfld= 3. * relative_pos * (mm[None, :] * relative_pos).sum(axis=1)[:, None] / \
               (relative_separation[i_box] ** 5)[:, None] - mm[None, :] / (relative_separation[i_box] ** 3)[:, None]

        if replace:
            self.data['bfld'][i_box,:] = bfld
        else:
            self.data['bfld'][i_box,:] += bfld

    def add_magnetic_field(self, seed_B=1000, object_R=None):
        c1 = (self.data['pos'] * self.data['mass'][:, None] * self.data['pass'][:, 0][:, None]).sum(axis=0) / self.m1
        print("c1 in magnetic field calculations= ", c1, " new_pos1= ", self.new_pos1)
        R1 = self.obj1.r if object_R is None else object_R
        R2 = self.obj2.r if object_R is None else object_R
        self.calculate_magnetic_field_around_pos(self.new_pos1, seed_B=seed_B, object_R=R1)
        print("added magnetic dipole inside object 1 with initial ", seed_B, " G at R_surface=", R1, " cm")
        self.calculate_magnetic_field_around_pos(self.new_pos2, seed_B=seed_B, object_R=R2, replace=False)
        print("added magnetic dipole inside object 2 with initial ", seed_B, " G at R_surface=", R2, " cm")

    def calculate_RL(self):
        q = self.m2 / self.m1
        RL2 = self.r2 * (0.6 * q ** (2. / 3.) + log(1 + q ** (1. / 3.))) / (0.49 * q ** (2. / 3.))

        q = self.m1 / self.m2
        RL1 = self.r1 * (0.6 * q ** (2. / 3.) + log(1 + q ** (1. / 3.))) / (0.49 * q ** (2. / 3.))

        return max(RL1, RL2)


class BinariesLoader(BinariesICs):
    def __init__(self, snapshot_file, conditional_axis=0, rhocut=1, species_file="species55.txt", load_types=[0]):
        print("initializeing binaries from a single snapshot")
        self.binary_snapshot = gadget_readsnapname(snapshot_file, hdf5=True, loadonlytype=load_types)
        if conditional_axis is not None:
            if rhocut is not None:
                self.get_objects_position_density(conditional_axis, rhocut)
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
        self.i1, = np.where(
            (self.binary_snapshot.pos[:, conditional_axis] - self.binary_snapshot.center[conditional_axis] > 0) &
            (self.binary_snapshot.rho > rhocut))
        self.i2, = np.where(
            (self.binary_snapshot.pos[:, conditional_axis] - self.binary_snapshot.center[conditional_axis] < 0) &
            (self.binary_snapshot.rho > rhocut))

    def get_objects_position(self, conditional_axis):
        self.i1, = np.where(
            (self.binary_snapshot.pos[:, conditional_axis] - self.binary_snapshot.center[conditional_axis] > 0))
        self.i2, = np.where(
            (self.binary_snapshot.pos[:, conditional_axis] - self.binary_snapshot.center[conditional_axis] < 0))

    def get_objects_density(self, rhocut):
        self.i1, = np.where((self.binary_snapshot.pass00 == 1) & (self.binary_snapshot.rho > rhocut))
        self.i2, = np.where((self.binary_snapshot.pass01 == 1) & (self.binary_snapshot.rho > rhocut))

    def get_objects(self):
        self.i1, = np.where((self.binary_snapshot.pass00 == 1))
        self.i2, = np.where((self.binary_snapshot.pass01 == 1))


class SeparateBinariesLoader(BinariesICs):
    def __init__(self, snapshot_file1, snapshot_file2, rhocut=1, species_file="species55.txt", load_types=[0]):
        self.snapshot1 = gadget_readsnapname(snapshot_file1, hdf5=True, loadonlytype=load_types)
        self.snapshot2 = gadget_readsnapname(snapshot_file2, hdf5=True, loadonlytype=load_types)
        super().__init__(species_file, rhocut=rhocut)

    def get_objects_density(self, rhocut):
        self.i1, = np.where((self.snapshot1.rho > rhocut))
        self.i2, = np.where((self.snapshot2.rho > rhocut))


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--snapshot_file', type=str, help='path to snapshot file', default="")
    parser.add_argument('--snapshot_file2', type=str, help='path to secind snapshot file', default="")
    parser.add_argument('--species_file', type=str, help='path to species file', default="species55.txt")
    parser.add_argument('--conditional_axis', type=int, help='axis of motion', default=0)
    parser.add_argument('--rhocut', type=float, help='lower cutoff of density', default=1)
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types', default=[0])
    parser.add_argument('--period', type=float, help='orbital period in seconds for mergers', default=0)
    parser.add_argument('--impact_parameter_rhocut', type=float,
                        help='calculate impact parameter according to this density cutoff', default=0)
    parser.add_argument('--impact_parameter', type=float,
                        help='use this initial impact parameter at the initial separation specified by either '
                             'separation or relative_to_RL, used only if not using impact_parameter_rho_cut', default=0)
    parser.add_argument('--relative_velocity', type=float, help='', default=None)
    parser.add_argument('--relative_to_RL', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='is the distance should be relative to RL size?',
                        default=True)
    parser.add_argument('--RL_factor', type=float, help='if relative to RL, by what factor?', default=2)
    parser.add_argument('--find_next_interaction', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='is the distance should be relative to RL size?',
                        default=False)
    parser.add_argument('--separation', type=float, help='initial separation between the binary objects', default=None)
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="bin.ic.dat")
    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    if args.snapshot_file2 is None:
        binary = BinariesLoader(args.snapshot_file, conditional_axis=args.conditional_axis,
                                rhocut=args.rhocut, species_file=args.species_file, load_types=args.load_types)
    else:
        binary = SeparateBinariesLoader(args.snapshot_file, args.snapshot_file2, args.rhocut, args.species_file,
                                        load_types=args.load_types)

    if args.period > 0:
        binary.create_ic_merger_keplerian(ic_file_name=args.ic_file_name, initial_period=args.period,
                                          rhocut=args.rhocut,
                                          relative_to_RL=args.relative_to_RL, factor=args.RL_factor)
    elif args.impact_parameter_rhocut > 0:
        i1, = np.where(binary.snapshot1.rho > args.impact_parameter_rhocut)
        i2, = np.where(binary.snapshot2.rho > args.impact_parameter_rhocut)
        b = binary.snapshot1.r()[i1].max() + binary.snapshot2.r()[i2].max()
        print("b = ", b)
        binary.create_ic_collision(b, args.ic_file_name, args.relative_velocity, args.separation,
                                   args.impact_parameter_rhocut)

    elif args.impact_parameter_rhocut > 0:
        b = args.impact_parameter
        print("b = ", b)
        binary.create_ic_collision(b, args.ic_file_name, args.relative_velocity, args.separation,
                                   args.impact_parameter_rhocut)

    elif args.find_next_interaction:
        binary.create_ic_for_next_interaction(args.ic_file_name, args.relative_to_RL, args.RL_factor, args.separation)
