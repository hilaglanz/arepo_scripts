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

def plot_profile_test(output_dir,snapshot_name,plotting_dir,testing_value="rho",snapshot_number_array=[0,8,10],
                      center=False, log=True,new_fig=True, around_objects=False, motion_axis= 0, object_num=0):
    if not os.path.exists(plotting_dir):
        os.mkdir(plotting_dir)

    if new_fig:
        set_new_fig_properties()

    line_colors = pylab.rcParams['axes.prop_cycle'].by_key()['color']
    labels = []
    suffix = ""
    for index, snapshot_number in enumerate(snapshot_number_array):
        if around_objects:
            snapshot_file = "%s/%s%03d" % (output_dir, snapshot_name, snapshot_number)
            if object_num != 0:
                binary = BinariesLoader(snapshot_file, conditional_axis=motion_axis)
                s = binary.binary_snapshot
                if object_num == 1:
                    i = binary.i1
                    print("doing object 1")
                else:
                    i = binary.i2
                    print("doing object 2")
            else:
                s = gadget_readsnap(snapshot_number, output_dir, snapshot_name)
                center = s.centerofmass()
                i = np.where(s.rho > 10)
                print("doing single object")
            if testing_value == "bfld" or testing_value == "B":
                print("adding magnetic field size")
                s.data["B"] = np.sqrt((s.data['bfld'] * s.data['bfld']).sum(axis=1))
                testing_value = "B"
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
            if testing_value == "bfld" or testing_value == "B":
                s.data["B"] = np.sqrt((s.data['bfld'] * s.data['bfld']).sum(axis=1))
                testing_value = "B"
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

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    if args.around_objects:
        if not args.take_single_object:
            plot_profile_test(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                              testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                              around_objects=args.around_objects, motion_axis=args.motion_axis, object_num=1)
            plot_profile_test(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                              testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                              around_objects=args.around_objects, motion_axis=args.motion_axis,object_num=2,new_fig=True)
    else:
        plot_profile_test(output_dir= args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                          testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                          around_objects=args.around_objects, motion_axis= args.motion_axis)
