import os, sys
import argparse
from loadmodules import *


def create_ic_with_sink(ic_path, boxsize=32, G=6.672*10**-8, mach=1.4, cs=1, rho=1, gamma=5.0/3, Ra=1, Rs=0.02, res=100):
    vel = mach*cs
    accretion_radius = Ra
    print("using cs= ", cs, "v_inf= ", vel, "mach= ", mach, "rho_inf= ", rho, "Ra= ", accretion_radius, "G= ", G)

    pointStar = {}
    pointStar['count'] = 1
    pointStar['pos'] = np.zeros( (1,3) )
    pointStar['pos'] += 0.5 * boxsize
    pointStar['mass'] = np.array([1.0*(vel**2)/(2*G)])
    pointStar['vel'] = np.zeros( (1,3) )
    pointStar['boxsize'] = boxsize

    gadget_add_grid(pointStar, accretion_radius, res=round(accretion_radius/Rs))
    print("added inner grid with size of Ra= ", accretion_radius)
    print("minimum vol =", (accretion_radius**3)/round(accretion_radius/Rs)**3)
    gadget_add_grids(pointStar, [2*accretion_radius, 4*accretion_radius, 8*accretion_radius, 16*accretion_radius,
                                 boxsize], res=res)
    print("added 5 sub-grids around")
    print("max vol =" , (boxsize**3.0 - (16.0*accretion_radius)**3)/res**3)
    print("mid vol =" , (16**3.0 - (8.0*accretion_radius)**3)/res**3)

    pointStar['type']=np.zeros(pointStar['count'])
    pointStar['type'][0] = 5
    pointStar['mass'][1:] = rho #3e-2 with read mass as density will give same densities to all subgrids cells
    pointStar['vel'][1:,0] = vel
    pointStar['u'] = np.zeros(pointStar['count'])
    pointStar['u'][:] = (cs**2)/(gamma*(gamma-1))
    print("u: ", (cs**2)/(gamma*(gamma-1)))

    pointStar['bflg']=np.zeros(pointStar['count'], dtype=np.uint32)
    pointStar['bflg'][0]=1
    print(pointStar.keys())

    print(pointStar['vel'][:,0])
    for key in pointStar.keys():
        print (key)
        if key != 'count' and key!= 'boxsize' and pointStar[key].shape[0]>1:
            ax = 0
        else:
            ax=None
        pointStar[key]= np.flip(pointStar[key], axis=ax)
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
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    create_ic_with_sink(ic_path=args.ic_path, boxsize=args.boxsize, G=args.G, mach=args.mach, cs=args.cs, rho=args.rho,
                        gamma=args.gamma, Ra=args.Ra, Rs=args.Rs, res=args.res)


