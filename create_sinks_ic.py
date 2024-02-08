import os, sys
import argparse
import numpy as np
from loadmodules import *
from stellar_ics import background_grid
from stellar_ics.tools import get_grid_u


def initialize_dictionary_with_point_masses(point_mass, npart, boxsize):
    pointStar = {}
    pointStar['count'] = npart
    pointStar['pos'] = np.zeros( (npart,3) )
    pointStar['pos'] += 0.5 * boxsize
    pointStar['mass'] = np.array([point_mass] * npart)
    pointStar['vel'] = np.zeros( (npart,3) )
    pointStar['boxsize'] = boxsize

    return pointStar

def get_finest_grid_size_and_resolution(accretion_radius, sink_radius, surroundings = 10):
    finest_grid_size = accretion_radius
    highest_resolution = round(finest_grid_size / (sink_radius * pi / surroundings))
    while highest_resolution > 200:
        finest_grid_size /= 2.0
        print("decreasing finest grid size to ", finest_grid_size)
        highest_resolution = round(finest_grid_size / (sink_radius * pi / surroundings))

    return finest_grid_size, highest_resolution
def get_smoothed_sub_grid_sizes(boxsize, finest_grid_size):
    num_sub_grids = np.log2(boxsize / finest_grid_size)
    sub_grid_indices = np.array(range(1, int(num_sub_grids) + 1))
    sub_grid_sizes = finest_grid_size * 2. ** sub_grid_indices
    if len(sub_grid_sizes) > 1 and boxsize > 1.4 * sub_grid_sizes[-1]:
        sub_grid_sizes = np.append(sub_grid_sizes, boxsize)
    else:
        sub_grid_sizes[-1] = boxsize

    return sub_grid_sizes

def create_hard_sphere_boundary(mass, radius, background_data, point_mass_id=0, factor_u=10**-12):
    position = background_data['pos'][point_mass_id]
    sphere_cells = np.where(np.sqrt(((background_data['pos']-position)**2).sum(axis=1)) < radius)
    sphere_cells = np.delete(sphere_cells, [point_mass_id])
    if mass == 0:
        background_data['mass'][sphere_cells] *= factor_u
    else:
        background_data['mass'][sphere_cells] = 3.0 * mass / (4 * pi * radius ** 3)
    background_data['u'][sphere_cells] *= factor_u
    background_data['vel'][sphere_cells,:] *= factor_u
    background_data['bflg'][sphere_cells] = 1
    background_data['id'][sphere_cells] += 10**9

    return background_data

def create_a_radial_gowing_mesh(inner_sphere_radius, outer_sphere_radius, smallest_cell_radius, growth_factor=1.1,
                                box_center=[0,0,0], resolution  = 100):
    growing_radius = outer_sphere_radius - inner_sphere_radius
    largest_cell_radius = ((growth_factor - 1)*growing_radius + smallest_cell_radius) / growth_factor
    print("largest cell radius with a growth factor of", growth_factor, "is ", largest_cell_radius)
    print("box_center= ", box_center)
    last_index = ceil(log(largest_cell_radius / smallest_cell_radius) / log(growth_factor))
    pos_array = []
    cell_radius = smallest_cell_radius / growth_factor
    current_distance = inner_sphere_radius - smallest_cell_radius
    x0 = box_center[0]
    y0 = box_center[1]
    z0 = box_center[2]
    for i in range(last_index):
        current_distance += cell_radius
        x = x0
        y = y0
        z = z0
        for phi in np.linspace(0, 2*pi, resolution): #2*arcsin(cell_radius/(2*current_distance))):
            for psi in np.linspace(0, pi, resolution/2): #2*arcsin(cell_radius/(2*current_distance))):
                x = x0 + current_distance*cos(phi)*sin(psi)
                y = y0 + current_distance*sin(phi)*sin(psi)
                z = z0 + current_distance*cos(psi)
                pos_array.append(np.array([x,y,z]))
        cell_radius *= growth_factor
    pos_array = np.array(pos_array) + box_center
    print(pos_array)
    print(pos_array[:,0].max(), pos_array[:,1].max(), pos_array[:,2].max())
    print(pos_array[:,0].min(), pos_array[:,1].min(), pos_array[:,2].min())

    return cell_radius, current_distance + y0, pos_array

def create_ic_with_sink(ic_path, boxsize=32, G=6.672*10**-8, mach=1.4, cs=1, rho=1, gamma=5.0/3, Ra=1, Rs=0.02, res=100,
                        binary=False, semimajor = 2.5, supersonic_perscription=True, surroundings=10, hard_sphere=False,
                        use_wind_ids_for_region=None):
    vel = mach*cs
    accretion_radius = Ra
    last_sink_i = 0
    orbital_vel = 0
    if supersonic_perscription:
        sink_mass = accretion_radius*(vel**2)/(2*G)
    else:
        sink_mass = accretion_radius*(vel**2 + cs**2)/(2*G)
    if binary:
        last_sink_i = 1
        orbital_vel = 0.5 * (2 * G * sink_mass / semimajor) ** 0.5
    num_sinks = last_sink_i + 1
    print("using cs= ", cs, "v_inf= ", vel, "mach= ", mach, "rho_inf= ", rho, "Ra= ", accretion_radius, "G= ", G)

    pointStar = initialize_dictionary_with_point_masses(sink_mass, num_sinks, boxsize)
    if not binary:
        #background = initialize_dictionary_with_point_masses(rho,1,Rs*0.1)
        maximum_cell_radius, sphere_size, pos = create_a_radial_gowing_mesh(Rs, minimum(100 * Rs, boxsize),
                                                                            Rs / surroundings, resolution=surroundings)
        data={}
        data["pos"] = pos
        data["boxsize"] = sphere_size
        pointStar['mass'] = np.array([rho] * pos.shape[0])
        pointStar['vel'] = np.zeros( (pos.shape[0],3) )
        #gadget_add_grid(background, Rs * 0.5, res=min([res, highest_resolution])) # no need for so many cells well inside the sink
        '''
        bgSphere = background_grid.BackgroundGridAroundSphere(background, boxsize=Rs*1.5, ndir=ceil(0.8*res),
                                                              newsize=Rs*2.0, grid_rho=rho,
                                                               grid_u=(cs**2)/(gamma*(gamma-1)))
        '''

        bgSphere = background_grid.BackgroundGridAroundSphere(data, boxsize=Rs * 1.5, ndir=ceil(0.8 * res),
                                                              newsize=sphere_size, grid_rho=rho,
                                                              grid_u=(cs ** 2) / (gamma * (gamma - 1)))
        background = bgSphere.add_grid()
        print("sphere_size= ", sphere_size, "background boxsize= ", background["boxsize"])
        print("added background grid of size ", background["boxsize"], " around the sphere of size ", 1.5*Rs)
        print("minimum volume= ", (Rs/surroundings)**3.0)
        #gadget_add_grid(background, background['boxsize'], res) #filling the inner sphere
        background['pos'] += boxsize/2.0 - 0.5 * background['boxsize']
        for key in pointStar.keys():
            if key == 'count':
                pointStar[key] += background[key]
                continue
            if key == 'boxsize':
                continue
            else:
                pointStar[key] = np.append(pointStar[key], background[key], axis=0)
        print(pointStar["pos"])
        #finest_grid_size, highest_resolution = get_finest_grid_size_and_resolution(accretion_radius, Rs, surroundings)
        finest_grid_size = sphere_size * 1.01
        highest_resolution = maximum_cell_radius
        #gadget_add_grid(pointStar, Rs * 0.8, res=ceil(mean([res, highest_resolution])))  # no need for so many cells well inside the sink
        gadget_add_grid(pointStar, finest_grid_size, res=ceil(highest_resolution*0.8)) # should have many close to its surface
    else:
        finest_grid_size, highest_resolution = get_finest_grid_size_and_resolution(accretion_radius, Rs, surroundings)
        gadget_add_grid(pointStar, finest_grid_size, res=ceil(highest_resolution)) # should have many close to its surface

    print("added inner grid with size of ", finest_grid_size / accretion_radius, "Ra")
    print("minimum vol =", (finest_grid_size ** 3) / highest_resolution ** 3)
    print(pointStar["pos"])

    sub_grid_sizes = get_smoothed_sub_grid_sizes(boxsize, finest_grid_size)
    gadget_add_grids(pointStar, sub_grid_sizes, res=res)
    print(sub_grid_sizes)
    print("added {0} sub-grids around".format(len(sub_grid_sizes)))
    print("max vol =", boxsize**3 / res ** 3.0)
    pointStar['type']=np.zeros(pointStar['count'])
    if num_sinks > 0:
        pointStar['type'][:num_sinks] = [5] * num_sinks
    pointStar['mass'][num_sinks:] = rho #3e-2 with read mass as density will give same densities to all subgrids cells
    pointStar['vel'][num_sinks:,0] = vel
    pointStar['u'] = np.zeros(pointStar['count'])
    pointStar['u'][:] = (cs**2)/(gamma*(gamma-1))
    print("u: ", (cs**2)/(gamma*(gamma-1)))

    pointStar['bflg'] = np.zeros(pointStar['count'])
    pointStar['bfld'] = np.zeros((pointStar['count'],3))
    if use_wind_ids_for_region is not None:
        pointStar['id'] = np.array([i + 1 for i in range(pointStar['count'])])
        pointStar['id'][np.where(pointStar['pos'][:,0] < use_wind_ids_for_region)] += 10**9
        pointStar['id'][np.where(pointStar['pos'][:, 0] > (pointStar['boxsize'] - use_wind_ids_for_region))] += 10 ** 9

    if hard_sphere:
        pointStar = create_hard_sphere_boundary(0, 0.8*Rs, pointStar, 0,1)

    print(pointStar.keys())
    if num_sinks > 0:
        if binary:
            pointStar['pos'][0, 0] += semimajor / 2.0
            pointStar['pos'][1, 0] -= semimajor / 2.0
            pointStar['vel'][0, 1] += orbital_vel
            pointStar['vel'][1, 1] -= orbital_vel

        print(pointStar['vel'][:,0])
        for key in pointStar.keys():
            print (key)
            if key != 'count' and key!= 'boxsize' and pointStar[key].shape[0]>1:
                ax = 0
            else:
                ax = None
            pointStar[key] = np.flip(pointStar[key], axis=ax)

    print(pointStar['vel'][:, 0])
    print(pointStar['mass'])

    gadget_write_ics(ic_path, pointStar, double=True, format="hdf5")
    print("ic file saved to ", ic_path)


def InitParser():
    parser = argparse.ArgumentParser(description='This file creates a new ic file for simulating '
                                                 'accretion on a sink particle moving inside a wind tunnel')
    parser.add_argument('--ic_path', type=str,  help='path to output file including file name',
                        default="tasupersonic.dat.ic")
    parser.add_argument('--boxsize', type=float,  help='size of the total cubical box', default=32)
    parser.add_argument('--cs', type=float,  help='sound speed of the medium in internal units', default=1)
    parser.add_argument('--mach', type=float,  help='mach number of the medium in internal units', default=1.4)
    parser.add_argument('--rho', type=float,  help='density of the medium in internal units', default=1.0)
    parser.add_argument('--Ra', type=float,  help='accretion radius in internal units', default=1.0)
    parser.add_argument('--Rs', type=float,  help='sink particle radius in internal units', default=0.02)
    parser.add_argument('--G', type=float,  help='gravitational constant internal units', default=6.672e-8)
    parser.add_argument('--gamma', type=float,  help='adiabatic index', default=5.0/3.0)
    parser.add_argument('--res', type=int,  help='subgrids resolution', default=100)
    parser.add_argument('--sink_surroundings', type=int, help='how many neighbours surround the sink on its plane', default=10)
    parser.add_argument('--binary', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='do we have a binary accreting?',
                        default=False)
    parser.add_argument('--binary_separation', type=float,  help='initial separation between the binary objects', default=None)
    parser.add_argument('--hard_sphere', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='do we have a hard sphere instead of a point mass?',
                        default=False)
    parser.add_argument('--use_wind_ids_region', type=float,  help='wind_large_ids_region', default=None)
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    create_ic_with_sink(ic_path=args.ic_path, boxsize=args.boxsize, G=args.G, mach=args.mach, cs=args.cs, rho=args.rho,
                        gamma=args.gamma, Ra=args.Ra, Rs=args.Rs, res=args.res, binary=args.binary,
                        semimajor=args.binary_separation, surroundings=args.sink_surroundings,
                        hard_sphere=args.hard_sphere, use_wind_ids_for_region=args.use_wind_ids_region)


