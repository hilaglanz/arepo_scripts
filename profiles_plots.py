import os
import glob
import argparse
import numpy as np
from loadmodules import *

def set_new_fig_properties():
    fig = figure(figsize=(36, 20))
    rcParams.update({'font.size': 40, 'font.family': 'Serif'})
    rcParams['text.usetex'] = True
    rcParams['lines.linewidth'] = 3.0

def plot_profile_test(output_dir,snapshot_name,plotting_dir,testing_value="rho",snapshot_number_array=[0,8,10],log=True,new_fig=True):
    if new_fig:
        set_new_fig_properties()
    evenly_spaced_interval = np.linspace(0, 1, len(snapshot_number_array))
    line_colors = [cm.rainbow(x) for x in evenly_spaced_interval]
    labels=[]
    for index, snapshot_number in enumerate(snapshot_number_array):
        s = gadget_readsnap(snapshot_number,output_dir,snapshot_name)
        s.plot_radprof(testing_value, log=log,color=line_colors[index])
        labels.append("snap " + str(snapshot_number) + "," + str(s.time)+ " [s]")
    if len(snapshot_number_array) > 1:
        legend(labels)
    filename = plotting_dir + "/" + value + "_profile_" + "_".join([str(snap_num) for snap_num in snapshot_number_array]) + ".png"
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

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    plot_profile_test(output_dir= args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                      testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot)