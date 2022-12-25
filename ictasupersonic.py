import os, sys
import argparse
import numpy as np
from loadmodules import *

def initialize_dictionary_with_point_masses(point_mass, npart, boxsize):
    pointStar = {}
    pointStar['count'] = npart
    pointStar['pos'] = np.zeros( (npart,3) )
    pointStar['pos'] += 0.5 * boxsize
    pointStar['mass'] = np.copy(point_mass)
    pointStar['vel'] = np.zeros( (npart,3) )
    pointStar['boxsize'] = boxsize

    return pointStar

def get_finest_grid_size_and_resolution(accretion_radius, sink_radius):
    finest_grid_size = accretion_radius
    highest_resolution = round(10.0 * finest_grid_size / sink_radius)
    while highest_resolution > 200:
        finest_grid_size /= 2.0
        print("decreasing finest grid size to ", finest_grid_size)
        highest_resolution = round(10.0 * finest_grid_size / sink_radius)

    return finest_grid_size, highest_resolution
def get_smoothed_sub_grid_sizes(boxsize, finest_grid_size):
    num_sub_grids = np.log2(boxsize / finest_grid_size)
    sub_grid_indices = np.array(range(1, int(num_sub_grids) + 1))
    sub_grid_sizes = finest_grid_size * 2. ** sub_grid_indices
    if boxsize > 1.4 * sub_grid_sizes[-1]:
        sub_grid_sizes = np.append(sub_grid_sizes, boxsize)
    else:
        sub_grid_sizes[-1] = boxsize

    return sub_grid_sizes

##Pass array to this function that contains data for all sinks(!)
def create_ic_with_sink_from_data(ic_path, sink_data, boxsize=32, G=6.672*10**-8, mach=1.4, cs=1, rho=1, gamma=5.0/3, Rs=0.02, res=100, supersonic_perscription=True):
    vel = mach*cs
    last_sink_i = 0
    orbital_vel = 0

    num_sinks = len(sink_data)
    ##First column corresponds to the masses
    sink_mass = sink_data[:,0]
    ##Accretion radius for the whole collection of sink particles
    if supersonic_perscription:
        accretion_radius=G*np.sum(sink_mass)/(vel**2)
    else: 
        accretion_radius=G*np.sum(sink_mass)/(vel**2+cs**2)
    ##Possibly make sure Rs and sink radius are always related...
    print("using cs= ", cs, "v_inf= ", vel, "mach= ", mach, "rho_inf= ", rho, "Ra= ", accretion_radius, "G= ", G)
    pointStar = initialize_dictionary_with_point_masses(sink_mass, num_sinks, boxsize)

    finest_grid_size, highest_resolution = get_finest_grid_size_and_resolution(accretion_radius, Rs)
    gadget_add_grid(pointStar, finest_grid_size, res=highest_resolution)
    print("added inner grid with size of ", finest_grid_size/accretion_radius, "Ra")
    print("minimum vol =", (finest_grid_size**3)/highest_resolution**3)

    sub_grid_sizes = get_smoothed_sub_grid_sizes(boxsize, finest_grid_size)
    gadget_add_grids(pointStar, sub_grid_sizes, res=res)
    print(sub_grid_sizes)
    print("added {0} sub-grids around".format(len(sub_grid_sizes)))
    print("max vol =", (boxsize**3.0 - (sub_grid_sizes[-2])**3)/res**3)
    
    pointStar['type']=np.zeros(pointStar['count'])
    pointStar['type'][:num_sinks] = [5] * num_sinks
    pointStar['pos'][:num_sinks]+=np.copy(sink_data[:,1:4])
    pointStar['vel'][:num_sinks]+=np.copy(sink_data[:,4:])

    pointStar['mass'][num_sinks:] = rho #3e-2 with read mass as density will give same densities to all subgrids cells
    pointStar['vel'][num_sinks:,0] = vel
    pointStar['u'] = np.zeros(pointStar['count'])
    pointStar['u'][:] = (cs**2)/(gamma*(gamma-1))
    print("u: ", (cs**2)/(gamma*(gamma-1)))
    print(pointStar.keys())


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

def create_ic_with_sink(ic_path, boxsize=32, G=6.672*10**-8, mach=1.4, cs=1, rho=1, gamma=5.0/3, Ra=1, Rs=0.02, res=100,
                        binary=False, semimajor = 2.5, supersonic_perscription=True):
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

    pointStar = initialize_dictionary_with_point_masses(np.array([sink_mass] * num_sinks), num_sinks, boxsize)

    finest_grid_size, highest_resolution = get_finest_grid_size_and_resolution(accretion_radius, Rs)
    gadget_add_grid(pointStar, finest_grid_size, res=highest_resolution)
    print("added inner grid with size of ", finest_grid_size/accretion_radius, "Ra")
    print("minimum vol =", (finest_grid_size**3)/highest_resolution**3)

    sub_grid_sizes = get_smoothed_sub_grid_sizes(boxsize, finest_grid_size)
    gadget_add_grids(pointStar, sub_grid_sizes, res=res)
    print(sub_grid_sizes)
    print("added {0} sub-grids around".format(len(sub_grid_sizes)))
    print("max vol =", (boxsize**3.0 - (sub_grid_sizes[-2])**3)/res**3)
    
    pointStar['type']=np.zeros(pointStar['count'])
    pointStar['type'][:num_sinks] = [5] * num_sinks

    pointStar['mass'][num_sinks:] = rho #3e-2 with read mass as density will give same densities to all subgrids cells
    pointStar['vel'][num_sinks:,0] = vel
    pointStar['u'] = np.zeros(pointStar['count'])
    pointStar['u'][:] = (cs**2)/(gamma*(gamma-1))
    print("u: ", (cs**2)/(gamma*(gamma-1)))
    print(pointStar.keys())
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
    parser.add_argument('--binary', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='do we have a binary accreting?',
                        default=False)
    parser.add_argument('--binary_separation', type=float,  help='initial separation between the binary objects', default=None)
    parser.add_argument('--data_file', type=str, help='if specified read the sink data from a file...')
    parser.add_argument('--super', type=bool, help='controlling the prescription for the accretion radius...')
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    print(args.data_file)

    if args.data_file:
        sink_data=np.genfromtxt(args.data_file)
        create_ic_with_sink_from_data(args.ic_path, sink_data, boxsize=args.boxsize, G=args.G, mach=args.mach, cs=args.cs, rho=args.rho,
            gamma=args.gamma, Rs=args.Rs, res=args.res, supersonic_perscription=args.super)
    else:
        create_ic_with_sink(ic_path=args.ic_path, boxsize=args.boxsize, G=args.G, mach=args.mach, cs=args.cs, rho=args.rho,
                        gamma=args.gamma, Ra=args.Ra, Rs=args.Rs, res=args.res, binary=args.binary,
                        semimajor=args.binary_separation, supersonic_perscription=args.super)


