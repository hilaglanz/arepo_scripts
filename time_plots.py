import os
import glob
import argparse
import numpy as np
import pylab
from loadmodules import *
from BinariesICs import *
from profiles_plots import *
from make_plots import *

def calculate_particle_value_diff_rate(snapshot, particle_index, value, old_val, old_time):
    return calculate_particle_value_diff(snapshot, particle_index, value, old_val) / (snapshot.time - old_time)

def calculate_particle_value_diff(snapshot, particle_index, value, old_val):
    return calculate_particle_value(snapshot, particle_index, value) - old_val

def calculate_particle_value(snapshot, particle_index, value):
    return snapshot.data[value][particle_index]

def calculate_mean_value(snapshot, value):
    if value in snapshot.data.keys():
        if len(snapshot.data[value].shape) > 1:
            return np.sqrt((snapshot.data[value] * snapshot.data[value]).sum(axis=1)).mean()

        return snapshot.data[value].mean()
def calculate_max_value(snapshot, value):
    if value in snapshot.data.keys():
        if len(snapshot.data[value].shape) > 1:
            return np.sqrt((snapshot.data[value] * snapshot.data[value]).sum(axis=1)).max()

        return snapshot.data[value].max()

def calculate_value(snapshot, value, sink_value=False, sink_id=0):
    if value not in snapshot.data.keys():
        if "drag" in value:
            snapshot.data[value] = snapshot.data["acce"] * snapshot.mass[:,None]

    if sink_value:
        sink_idks = np.where(snapshot.type == 5)
        sink_idk = sink_idks[0][sink_id]
        print("doing for sink particle")
        if len(snapshot.data[value].shape) > 1:
            curr_val = np.sqrt((snapshot.data[value][sink_idk,:]**2).sum())
        else:
            curr_val = abs(snapshot.data[value][sink_idk])

        return curr_val

    values = (snapshot.data[value] * snapshot.mass[:, None]).sum(axis=0) / snapshot.mass.sum()
    return np.sqrt((values ** 2).sum(axis=1))

def calculate_value_over_time(snapshots_number_list, snapshot_dir="output", value="mass",
                    mean=False, max=False, sink_value=False, sink_id=0):
    value_over_time = []
    times = []
    value_to_calc = value
    for snapshot_num in snapshots_number_list:
        snapshot = gadget_readsnap(snapshot_num, snapshot_dir)
        times.append(snapshot.time)
        if "dot" in value:
            print("plotting time different")
            value_to_calc = value.split("dot")[0]
        if "diff" in value:
            print("plotting difference of value")
            value_to_calc = value.split("diff")[0]
        if mean:
            value_over_time.append(calculate_mean_value(snapshot, value_to_calc))
        else:
            if max:
                value_over_time.append(calculate_max_value(snapshot, value_to_calc))
            else:
                value_over_time.append(calculate_value(snapshot, value_to_calc, sink_value, sink_id))
    print("added ", value_to_calc, " to the time evolution")
    if value_to_calc != value:
        value_diff_over_time = []
        prev_value = value_over_time[0]
        prev_time = times[0]
        for i in range(len(value_over_time)):
            if "diff" in value:
                value_diff_over_time.append(value_over_time[i] - prev_value)
            else:
                if times[i] == prev_time:
                    value_diff_over_time.append(0)
                else:
                    value_diff_over_time.append((value_over_time[i] - prev_value) / (times[i] - prev_time))
            prev_value = value_over_time[i]
            prev_time = times[i]

        return value_diff_over_time, times

    else:
        return value_over_time, times

def make_time_plots(snapshots_number_list, snapshot_dir="output", plotting_dir="times_plots", value="mass", log=False,
                    mean=False, max=False,sink_value=False, sink_id=0):

    value_over_time, times = calculate_value_over_time(snapshots_number_list,
                                                       snapshot_dir, value, mean, max,
                                                       sink_value, sink_id)
    set_new_fig_properties()
    if log:
        pylab.semilogy(times, value_over_time)
    else:
        pylab.plot(times, value_over_time)
    suptitle(value + " time evolution", fontsize='x-large')
    rcParams.update({'font.size': 40, 'font.family': 'Serif'})
    rcParams['text.usetex'] = True
    filename = plotting_dir + "/" + value + "_over_time_" + \
               str(snapshots_number_list[0]) + "_to_" + str(snapshots_number_list[-1]) + ".png"
    print("saving to: ", filename)
    savefig(filename)
    print("saved fig")
    txt_file_name = filename.replace("png", "txt")
    print("saving txt file to: ", txt_file_name)
    with open(txt_file_name, "w") as opened_file:
        times_and_values = [str(times[i]) + "," + str(value_over_time[i]) + "\r\n" for i in range(len(value_over_time))]
        opened_file.writelines(times_and_values)

    print("saved txt")






def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--output_dir', type=str, help='diretory containing the snapshots', default="output")
    parser.add_argument('--snapshot_name', type=str, help='snapshots initilas', default="snapshot_")
    parser.add_argument('--beginStep', type=int,  help='first step', default=0)
    parser.add_argument('--lastStep', type=int,  help='last step', default=-1)
    parser.add_argument('--skipStep', type=int, help='number of steps to skip', default=1)
    parser.add_argument('--plotting_dir', type=str, help='diretory to plot to', default="plots")
    parser.add_argument('--value', type=str, help='value to plot and compare profile', default="rho")
    parser.add_argument('--logplot', type=lambda x: (str(x).lower() in ['true', '1', 'yes']), help='logplot',
                        default=True)
    parser.add_argument('--mean', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='calculate mean value', default=False)
    parser.add_argument('--max', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='calculate max value (only possible if not using mean value)', default=False)
    parser.add_argument('--sink_value', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='calculate value for sink', default=False)
    parser.add_argument('--sink_id', type=int, help='sink particle to plot for', default=0)

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    snapshot_number_list = get_snapshot_number_list(args.output_dir, snapshotName=args.snapshot_name,
                                                    firstSnap=args.beginStep, lastSnap=args.lastStep,
                                                    skipSteps= args.skipStep)
    if not os.path.exists(args.plotting_dir):
        os.mkdir(args.plotting_dir)

    make_time_plots(snapshot_number_list, snapshot_dir=args.output_dir, plotting_dir=args.plotting_dir, value=args.value,
                    log=args.logplot, mean=args.mean, max=args.max, sink_value=args.sink_value, sink_id=args.sink_id )
