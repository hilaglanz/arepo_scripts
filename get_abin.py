import glob
import h5py
import numpy as np
import rebound
import cgs_const as cgs
import argparse
from loadmodules import *

parser = argparse.ArgumentParser(description='Get orbital elements from output')
parser.add_argument('--out_path', type=str,  help='path to output files',
                        default="output")
parser.add_argument('--type', type=int,  help='Particle type for binary', default=5)
args=parser.parse_args()
out_path=args.out_path
ptype=args.type


pos_all1=[]
pos_all2=[]
abins=[]
ebins=[]
sep=[]
snaps=glob.glob(out_path+'/snapshot*hdf5')
print(snaps)
ts=[]



for ii in range(len(snaps)):
    sim=rebound.Simulation()
    sim.G=cgs.G
    snap=gadget_readsnap(ii, out_path)
    with h5py.File(out_path+'/snapshot_{0:03d}.hdf5'.format(ii)) as ff:
        pos12=ff['PartType'+str(ptype)]['Coordinates'][...]
        vel12=ff['PartType'+str(ptype)]['Velocities'][...]
        mass12=ff['PartType'+str(ptype)]['Masses'][...]
        ids=ff['PartType'+str(ptype)]['ParticleIDs'][...]
        id1=(np.where(ids==ids[0])[0][0])
        id2=(np.where(ids==ids[-1])[0][0])
        pos_all1.append(pos12[id1])
        pos_all2.append(pos12[id2])
    print(id1, id2)
    print(pos12)
    print(vel12)
    sim.add(m=mass12[id1], x=pos12[id1][0], y=pos12[id1][1], z=pos12[id1][2],\
           vx=vel12[id1][0], vy=vel12[id1][1], vz=vel12[id1][2])
    sim.add(m=mass12[id2], x=pos12[id2][0], y=pos12[id2][1], z=pos12[id2][2],\
           vx=vel12[id2][0], vy=vel12[id2][1], vz=vel12[id2][2])
    orbs=sim.calculate_orbits(primary=sim.particles[0])
    abins.append(orbs[0].a)
    ebins.append(orbs[0].e)
    sep.append(np.linalg.norm(np.array(pos12[id1])-np.array(pos12[id2])))
    ts.append(snap.time)
    
pos_all1=np.array(pos_all1)
pos_all2=np.array(pos_all2)
abins=np.array(abins)
ebins=np.array(ebins)

np.savetxt('abin',np.transpose([ts, abins]))
np.savetxt('ebin',np.transpose([ts, ebins]))
np.savetxt('sep_bin',np.transpose([ts, sep]))
