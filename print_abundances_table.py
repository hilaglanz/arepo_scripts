import numpy as np
from loadmodules import *
import sys,os
import argparse
from make_plots import *

def get_abundances_table(snap, low_dens=1e2):
    i, = np.where(snap.rho > low_dens)
    vel_bound = (snap.vel[i, :] * snap.mass.astype('f8')[i, None]).sum(axis=0) / snap.mass[i].astype('f8').sum()
    bound, = np.where((snap.pot + 0.5 * ((snap.vel - vel_bound[None, :]) ** 2.).sum(axis=1) + snap.u <= 0.))
    unbound, = np.where((snap.pot + 0.5 * ((snap.vel - vel_bound[None, :]) ** 2.).sum(axis=1) + snap.u >= 0.))
    bound_abundances = np.zeros(55)
    unbound_abundances = np.zeros(55)
    element_names = np.array([
        '${}^{}$n',
        '${}^{}$H',
        '${}^{4}$He',
        '${}^{11}$B',
        '${}^{12}$C',
        '${}^{13}$C',
        '${}^{13}$N',
        '${}^{14}$N',
        '${}^{15}$N',
        '${}^{15}$O',
        '${}^{16}$O',
        '${}^{17}$O',
        '${}^{18}$F',
        '${}^{19}$Ne',
        '${}^{20}$Ne',
        '${}^{21}$Ne',
        '${}^{22}$Ne',
        '${}^{23}$Ne',
        '${}^{23}$Mg',
        '${}^{24}$Mg',
        '${}^{25}$Mg',
        '${}^{26}$Mg',
        '${}^{25}$Al',
        '${}^{26}$Al',
        '${}^{27}$Al',
        '${}^{28}$Si',
        '${}^{29}$Si',
        '${}^{30}$Si',
        '${}^{29}$P',
        '${}^{30}$P',
        '${}^{31}$P',
        '${}^{31}$S',
        '${}^{32}$S',
        '${}^{33}$S',
        '${}^{33}$Cl',
        '${}^{34}$Cl',
        '${}^{353}$Cl',
        '${}^{36}$Ar',
        '${}^{37}$Ar',
        '${}^{38}$Ar',
        '${}^{39}$Ar',
        '${}^{39}$K',
        '${}^{40}$Ca',
        '${}^{43}$Sc',
        '${}^{44}$Ti',
        '${}^{47}$V',
        '${}^{48}$Cr',
        '${}^{51}$Mn',
        '${}^{52}$Fe',
        '${}^{56}$Fe',
        '${}^{55}$Co',
        '${}^{56}$Ni',
        '${}^{58}$Ni',
        '${}^{59}$Ni'])
    sum_table = []
    for element in range(55):
        bound_abundances[element] = (snap.xnuc[bound,element]*snap.mass[bound]).sum()/msol
        unbound_abundances[element] = (snap.xnuc[unbound,element]*snap.mass[unbound]).sum()/msol

    for element in range(55):
        sum_table.append(element_names[element] + r' & \num{' + str(bound_abundances[element]) + r'}' + r' & \num{' +
                         str(unbound_abundances[element]) + r'}\\')

    return sum_table

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--source_dir', type=str,  help='path to snapshot files directory', default= sys.argv[0])
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "plots")
    parser.add_argument('--snap_num', type=int,  help='number of the snapshot to analyze', default=0)
    parser.add_argument('--low_dens', type=float, help='cutoff density', default=1e2)

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    s = gadget_readsnap(args.snap_num,args.source_dir,loadonlytype=[0])
    sum_table = get_abundances_table(s,low_dens=args.low_dens)
    print("Writing elements")
    with open(args.saving_dir+'/abundances'+'_'+str(args.snap_num)+'.txt') as out_file:
        out_file.writelines(sum_table)
    print("Done.")
