import os, sys
import argparse
import numpy as np
from loadmodules import *


def create_ic_with_sink(ic_path, boxsize=32, G=6.672*10**-8, mach=1.4, cs=1, rho=1, gamma=5.0/3, Ra=1, Rs=0.02, res=100,
                        binary=False, semimajor = 2.5):
    vel = mach*cs
    accretion_radius = Ra
    last_sink_i = 0
    orbital_vel = 0
    sink_mass = 1.0*(vel**2)/(2*G)
    if binary:
        last_sink_i = 1
        orbital_vel = 0.5 * (2 * G * sink_mass / semimajor) ** 0.5
    num_sinks = last_sink_i + 1
    print("using cs= ", cs, "v_inf= ", vel, "mach= ", mach, "rho_inf= ", rho, "Ra= ", accretion_radius, "G= ", G)

    pointStar = {}
    pointStar['count'] = num_sinks
    pointStar['pos'] = np.zeros( (num_sinks,3) )
    pointStar['pos'] += 0.5 * boxsize
    pointStar['mass'] = np.array([sink_mass] * num_sinks)
    pointStar['vel'] = np.zeros( (num_sinks,3) )
    pointStar['boxsize'] = boxsize

    gadget_add_grid(pointStar, accretion_radius, res=round(accretion_radius/Rs))
    print("added inner grid with size of Ra= ", accretion_radius)
    print("minimum vol =", (accretion_radius**3)/round(accretion_radius/Rs)**3)
    zones=np.log(boxsize/accretion_radius)/np.log(2)
    grids=2.**np.array(range(1,int(zones)+1))*accretion_radius
    if boxsize>1.4*grids[-1]:
        grids=np.append(grids, boxsize)
    else:
        grids[-1]=boxsize
    gadget_add_grids(pointStar, grids, res=res)
    print(grids)
    print("added {0} sub-grids around".format(len(grids)))
    print("max vol =", (boxsize**3.0 - (grids[-2])**3)/res**3)
    
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
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    create_ic_with_sink(ic_path=args.ic_path, boxsize=args.boxsize, G=args.G, mach=args.mach, cs=args.cs, rho=args.rho,
                        gamma=args.gamma, Ra=args.Ra, Rs=args.Rs, res=args.res, binary=args.binary,
                        semimajor=args.binary_separation)


