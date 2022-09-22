import os
import glob
import argparse
import numpy as np
import pylab
from loadmodules import *
from BinariesICs import *

def set_new_fig_properties():
    fig = figure(figsize=(36, 20))
    rcParams.update({'font.size': 40, 'font.family': 'Serif'})
    rcParams['text.usetex'] = True
    rcParams['lines.linewidth'] = 3.0

def compute_cumulative_mass(snapshot, center):
    rsort = snapshot.r(center=center).argsort()

    mcum = np.zeros(snapshot.npart)
    mcum[0] = snapshot.mass[rsort[0]]
    for i in range(1, snapshot.npart):
        mcum[rsort[i]] = mcum[rsort[i - 1]] + snapshot.mass[rsort[i]]
    snapshot.data['cum_mass'] = mcum / msol
    return

def compute_gas_to_magnetic_pressure(snapshot):
    if "B" not in snapshot.data.keys():
        compute_value(snapshot, "B")
    snapshot.data["betta"] = 8.0 * pi * snapshot.pres/snapshot.data["B"]**2
    print("average betta = ", 8.0 * pi * snapshot.pres.mean()/(snapshot.data["B"]**2).mean())

    return "betta"

def compute_magnetic_stress_parameter(snapshot, center): #Shakura Sunyaev parameter
    if "B_R" not in snapshot.data.keys():
        phi = np.arctan(snapshot.pos[:,1]/snapshot.pos[:,0])
        r = np.sqrt((snapshot.pos[:,0] - center[0]) ** 2 + (snapshot.pos[:,1] - center[1]) ** 2)
        B_R = snapshot.bfld[:,0] * np.cos(phi) + snapshot.bfld[:,1] * np.sin(phi)
        B_phi = -snapshot.bfld[:,0] * np.sin(phi) + snapshot.bfld[:,1] * np.cos(phi)
        snapshot.data["B_R"] = B_R
        snapshot.data["B_phi"] = B_phi
        snapshot.data["B_z"] = snapshot.bfld[:,2]
    snapshot.data["alpha_m"] = -snapshot.data["B_R"] * snapshot.data["B_phi"]/(4 * pi * snapshot.pres)
    print("average alpha_m = ", -snapshot.data["B_R"].mean() * snapshot.data["B_phi"].mean()/ (4 * pi * snapshot.pres.mean()))

    return "alpha_m"

def compute_value(s, testing_value, center=None):
    if testing_value == "bfld" or testing_value == "B":
        print("adding magnetic field size")
        s.data["B"] = np.sqrt((s.data['bfld'] * s.data['bfld']).sum(axis=1))
        testing_value = "B"

    if testing_value == "betta":
        return compute_gas_to_magnetic_pressure(s)

    if center is None:
        center = s.center

    if testing_value == "alpha_m":
        return compute_magnetic_stress_parameter(s, center)

    if testing_value == "cum_mass":
        print("adding cummulative nass")
        compute_cumulative_mass(s, center)

    return testing_value

def plot_profile_test(output_dir,snapshot_name,plotting_dir,testing_value="rho",snapshot_number_array=[0,8,10],
                      center=False, log=True,new_fig=True, around_objects=False, around_density_peak=False, motion_axis= 0, object_num=0):
    if not os.path.exists(plotting_dir):
        os.mkdir(plotting_dir)

    if new_fig:
        set_new_fig_properties()

    line_colors = pylab.rcParams['axes.prop_cycle'].by_key()['color']
    labels = []
    suffix = ""
    for index, snapshot_number in enumerate(snapshot_number_array):
        if around_objects:
            if object_num != 0:
                snapshot_file = "%s/%s%03d" % (output_dir, snapshot_name, snapshot_number)
                binary = BinariesLoader(snapshot_file, conditional_axis=motion_axis)
                s = binary.binary_snapshot
                if object_num == 1:
                    i = binary.i1
                    center = binary.pos1
                    suffix = "1"
                    print("doing object 1")
                else:
                    i = binary.i2
                    center = binary.pos2
                    suffix = "2"
                    print("doing object 2")
            else:
                print("doing single object")
                s = gadget_readsnap(snapshot_number, output_dir, snapshot_name)
                if around_density_peak:
                    print("calculating around density peak")
                    center = s.pos[np.where(s.rho == s.rho.max())][-1, :]
                else:
                    center = s.centerofmass()
                print("around center: ", center)
                i = np.where(s.rho > 100)
            testing_value = compute_value(s,testing_value, center)
            nshells = 200
            dr = 0
            p = calcGrid.calcRadialProfile(s.data['pos'].astype('float64')[i],
                                           s.data[testing_value].astype('float64')[i], 2, nshells, dr,
                                           center[0],
                                           center[1], center[2])
            print("plotting")
            print("color= ", line_colors[index])
            if log:
                pylab.semilogy(p[1, :], p[0, :], color=line_colors[index])
            else:
                pylab.plot(p[1, :], p[0, :], color=line_colors[index])
        else:
            s = gadget_readsnap(snapshot_number, output_dir, snapshot_name)
            testing_value = compute_value(s, testing_value, center)
            s.plot_radprof(testing_value, log=log, color=line_colors[index], center=center)
        print("used color: ", line_colors[index], s.time)
        labels.append("snap " + str(snapshot_number) + "," + str(round(s.time, 2)) + " [s]")
    if len(snapshot_number_array) > 1:
        legend(labels)
    filename = plotting_dir + "/" + testing_value + "_profile_" + suffix + "_".join([str(snap_num) for snap_num
                                                                                     in snapshot_number_array]) + ".png"
    print("saving to: ", filename)
    savefig(filename)
    print("saved fig")

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--output_dir', type=str, help='diretory containing the snapshots', default="output")
    parser.add_argument('--snapshot_name', type=str, help='snapshots initilas', default="snapshot_")
    parser.add_argument('--plotting_dir', type=str, help='diretory to plot to', default="plots")
    parser.add_argument('--value', type=str, help='value to plot and compare profile', default="rho")
    parser.add_argument('--snapshot_nums', type=int, nargs='+', help='array of snapshot numbers to compare',
                        default=[0 ,8, 10])
    parser.add_argument('--logplot', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='logplot',
                        default=True)
    parser.add_argument('--around_objects', type=lambda x: (str(x).lower() in ['true', '1', 'yes']), 
                        help='should plot around each of the object in the binary system',
                        default=False)
    parser.add_argument('--motion_axis', type=int, help='axis of motion when plotting around objects', default=None)
    parser.add_argument('--take_single_object', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot around one  of the object',
                        default=False)
    parser.add_argument('--around_density_peak', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should calculate the center according to the density peak?',
                        default=False)

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    if args.around_objects and not args.take_single_object:
        plot_profile_test(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                          testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                          around_objects=args.around_objects, motion_axis=args.motion_axis,
                          around_density_peak=args.around_density_peak, object_num=1)
        plot_profile_test(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                          testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                          around_objects=args.around_objects, motion_axis=args.motion_axis,
                          around_density_peak=args.around_density_peak, object_num=2,new_fig=True)
    else:
        plot_profile_test(output_dir= args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                          testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                          around_objects=args.around_objects, motion_axis=args.motion_axis,
                          around_density_peak=args.around_density_peak)
