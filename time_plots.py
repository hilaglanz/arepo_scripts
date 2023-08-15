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

def get_relevant_plotting_parameters_around_object(around_density_peak, motion_axis, object_num, output_dir,
                                                   snapshot_name, snapshot_number):
    suffix = ""
    if object_num != 0:
        s, i, center, suffix = get_relevant_plotting_parameters_for_binary(motion_axis, object_num, output_dir,
                                                                           snapshot_name, snapshot_number)
    else:
        s, i, center = get_relevant_plotting_parameters_for_single(around_density_peak, output_dir,
                                                                   snapshot_name, snapshot_number)
    return s, i, center, suffix
def get_relevant_plotting_parameters_for_binary(motion_axis, object_num, output_dir, snapshot_name, snapshot_number):
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
    return s, i, center, suffix
def get_relevant_plotting_parameters_for_single(around_density_peak, output_dir, snapshot_name, snapshot_number):
    print("doing single object")
    s = gadget_readsnap(snapshot_number, output_dir, snapshot_name)
    if around_density_peak:
        print("calculating around density peak")
        center = s.pos[np.where(s.rho == s.rho.max())][-1, :]
    else:
        center = s.centerofmass()
    print("around center: ", center)
    i = np.where(s.rho > 100)
    return s, i, center

def calculate_mean_value(snapshot, value, ind=[]):
    if len(ind) == 0:
        ind = snapshot.data['mass'] != 0
    if value in snapshot.data.keys():
        if len(snapshot.data[value].shape) > 1:
            return np.sqrt((snapshot.data[value][ind] * snapshot.data[value][ind]).sum(axis=1)).mean()

        return snapshot.data[value][ind].mean()
def calculate_max_value(snapshot, value, ind=[]):
    if len(ind) == 0:
        ind = snapshot.data['mass'] != 0
    if value in snapshot.data.keys():
        if len(snapshot.data[value].shape) > 1:
            return np.sqrt((snapshot.data[value][ind] * snapshot.data[value][ind]).sum(axis=1)).max()

        return snapshot.data[value][ind].max()

def calculate_value(snapshot, value, sink_value=False, sink_id=0, ind=[]):
    if len(ind) == 0:
        ind = snapshot.data['mass'] != 0
    if value not in snapshot.data.keys():
        if "drag" in value:
            snapshot.data[value] = snapshot.data["acce"][ind] * snapshot.mass[ind,None]

    if sink_value:
        sink_idks = np.where(snapshot.type == 5)
        sink_idk = sink_idks[0][sink_id]
        print("doing for sink particle")
        if len(snapshot.data[value].shape) > 1:
            curr_val = np.sqrt((snapshot.data[value][sink_idk,:]**2).sum())
        else:
            curr_val = abs(snapshot.data[value][sink_idk])

        return curr_val

    if len(snapshot.data[value].shape) > 1:
        values = (snapshot.data[value][ind] * snapshot.mass[ind, None]).sum(axis=0) / snapshot.mass[ind].sum()
        return np.sqrt((values ** 2).sum(axis=1))
    else:
        values = (snapshot.data[value][ind] * snapshot.mass[ind]).sum(axis=0) / snapshot.mass[ind].sum()
        return np.sqrt((values ** 2).sum())


def calculate_value_over_time(snapshots_number_list, snapshot_dir="output", value="mass",
                    mean=False, max=False, sink_value=False, sink_id=0, around_objects=False, around_density_peak=False,
                              object_num=0, motion_axis=0, along_axis_line=False, relative_to_motion=None):
    value_over_time = []
    times = []
    value_to_calc = value
    for snapshot_num in snapshots_number_list:
        if around_objects:
            snapshot, cell_indices, center, suffix = get_relevant_plotting_parameters_around_object(around_density_peak,
                                                                                             motion_axis,
                                                                                             object_num, snapshot_dir,
                                                                                             "snapshot_",
                                                                                             snapshot_num)
            #TODO: currently center is ignored- and is plotted around the center of the box
            center = snapshot.center

        else:
            snapshot = gadget_readsnap(snapshot_num, snapshot_dir)
            cell_indices = snapshot.data['mass'] != 0
            suffix = ""
            center = snapshot.center
        times.append(snapshot.time)
        if "dot" in value:
            print("plotting time different")
            value_to_calc = value.split("dot")[0]
        if "diff" in value:
            print("plotting difference of value")
            value_to_calc = value.split("diff")[0]
        if along_axis_line:
            relevant_cells = np.where((absolute(s.pos[:, (motion_axis + 1) % 3] - center[(motion_axis + 1) % 3]) < (
                        2 * snapshot.data["vol"] ** (1.0 / 3))) &
                                      (absolute(snapshot.pos[:, (motion_axis + 2) % 3] - center[(motion_axis + 2) % 3]) < (
                                                  2 * snapshot.data["vol"] ** (1.0 / 3))))
            cell_indices = np.intersect1d(cell_indices, relevant_cells)

        if relative_to_motion is not None:
            relevant_vector = calculate_value(snapshot=snapshot, value="vel", ind=cell_indices)
            snapshot = calculate_value_relative_to_vector(snapshot, value_to_calc, relevant_vector, cell_indices)
            if relative_to_motion == 0:
                value_to_calc += "_v"
            else:
                value_to_calc += "_u"
        if mean:
            value_over_time.append(calculate_mean_value(snapshot, value_to_calc, ind=cell_indices))
        else:
            if max:
                value_over_time.append(calculate_max_value(snapshot, value_to_calc, ind=cell_indices))
            else:
                value_over_time.append(calculate_value(snapshot, value_to_calc, sink_value, sink_id, ind=cell_indices))
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
                    mean=False, max=False,sink_value=False, sink_id=0, around_objects=False, motion_axis=0, object_num=0,
                    along_axis_line=False, relative_to_motion=None):

    value_over_time, times = calculate_value_over_time(snapshots_number_list,
                                                       snapshot_dir, value, mean, max,
                                                       sink_value, sink_id, around_objects=around_objects,
                                                       object_num=object_num, motion_axis=motion_axis,
                                                       along_axis_line=along_axis_line,
                                                       relative_to_motion=relative_to_motion)
    set_new_fig_properties()
    if log:
        pylab.semilogy(times, value_over_time)
    else:
        pylab.plot(times, value_over_time)
    suptitle(value + " time evolution", fontsize='x-large')
    rcParams.update({'font.size': 40, 'font.family': 'Serif'})
    rcParams['text.usetex'] = True
    if object_num != 0:
        suffix = str(object_num)
    else:
        suffix = ""
    if along_axis_line:
        suffix += "_along_" + chr(ord('x') + motion_axis)
    filename = plotting_dir + "/" + value + suffix + "_over_time_" + \
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
    parser.add_argument('--around_objects', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot around each of the object in the binary system',
                        default=False)
    parser.add_argument('--motion_axis', type=int, help='axis of motion when plotting around objects', default=0)
    parser.add_argument('--along_axis', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='considering only values along the line of motion', default=None)
    parser.add_argument('--relative_to_motion', type=int,
                        help='default is none, if 0- value is along the motion axis 1- perpendicular to the motion axis',
                        default=None)
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
    snapshot_number_list = get_snapshot_number_list(args.output_dir, snapshotName=args.snapshot_name,
                                                    firstSnap=args.beginStep, lastSnap=args.lastStep,
                                                    skipSteps= args.skipStep)
    if not os.path.exists(args.plotting_dir):
        os.mkdir(args.plotting_dir)

    if args.around_objects and not args.take_single_object:
        make_time_plots(snapshot_number_list, snapshot_dir=args.output_dir, plotting_dir=args.plotting_dir, value=args.value,
                        log=args.logplot, mean=args.mean, max=args.max, sink_value=args.sink_value, sink_id=args.sink_id,
                        around_objects=args.around_objects, motion_axis=args.motion_axis, object_num=1,
                        relative_to_motion=args.relative_to_motion)
        make_time_plots(snapshot_number_list, snapshot_dir=args.output_dir, plotting_dir=args.plotting_dir, value=args.value,
                        log=args.logplot, mean=args.mean, max=args.max, sink_value=args.sink_value, sink_id=args.sink_id,
                        around_objects=args.around_objects, motion_axis=args.motion_axis, object_num=2,
                        relative_to_motion=args.relative_to_motion)
    else:
        make_time_plots(snapshot_number_list, snapshot_dir=args.output_dir, plotting_dir=args.plotting_dir, value=args.value,
                        log=args.logplot, mean=args.mean, max=args.max, sink_value=args.sink_value, sink_id=args.sink_id,
                        relative_to_motion=args.relative_to_motion)
