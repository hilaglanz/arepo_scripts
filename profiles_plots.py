import os
import glob
import argparse
import numpy as np
import pylab
from loadmodules import *
from BinariesICs import *
from plot_multiples import *


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

def compute_cumulative_unbounded_mass(snapshot, center):
    if "unbounded_mass" not in snapshot.data:
        compute_unbounded_mass(snapshot)

    rsort = snapshot.r(center=center).argsort()

    mcum = np.zeros(snapshot.npart)
    mcum[0] = snapshot.data["unbounded_mass"][rsort[0]]
    for i in range(1, snapshot.npart):
        mcum[rsort[i]] = mcum[rsort[i - 1]] + snapshot.data["unbounded_mass"][rsort[i]]
    snapshot.data['cum_unbounded_mass'] = mcum / msol
    return

def compute_gas_to_magnetic_pressure(snapshot):
    if "B" not in snapshot.data.keys():
        compute_value(snapshot, "B")
    snapshot.data["betta"] = 8.0 * pi * snapshot.pres / snapshot.data["B"] ** 2
    print("average betta = ", 8.0 * pi * snapshot.pres.mean() / (snapshot.data["B"] ** 2).mean())

    return "betta"


def compute_magnetic_stress_parameter(snapshot, center):  # Shakura Sunyaev parameter
    if "B_R" not in snapshot.data.keys():
        phi = np.arctan(snapshot.pos[:, 1] / snapshot.pos[:, 0])
        r = np.sqrt((snapshot.pos[:, 0] - center[None, 0]) ** 2 + (snapshot.pos[:, 1] - center[None, 1]) ** 2)
        print(len(r))
        print("max distance = ", r.max())
        B_R = snapshot.bfld[:, 0] * np.cos(phi) + snapshot.bfld[:, 1] * np.sin(phi)
        B_phi = -snapshot.bfld[:, 0] * np.sin(phi) + snapshot.bfld[:, 1] * np.cos(phi)
        snapshot.data["B_R"] = B_R
        snapshot.data["B_phi"] = B_phi
        snapshot.data["B_z"] = snapshot.bfld[:, 2]
    snapshot.data["alpha_m"] = -snapshot.data["B_R"] * snapshot.data["B_phi"] / (4 * pi * snapshot.pres)
    print("average alpha_m = ",
          -snapshot.data["B_R"].mean() * snapshot.data["B_phi"].mean() / (4 * pi * snapshot.pres.mean()))

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

    if testing_value == "mach":
        print("adding mach")
        s.computeMach()

    if testing_value == "cs" or testing_value == "csound" or testing_value == "sound":
        if "sound" not in s.data.keys():
            print("adding mach and cs")
            s.computeMach()

    if testing_value == "entr":
        if "GAMMA" in s.config.keys():
            s.data["entr"] = s.data["pres"] / s.data["rho"] ** s.config["GAMMA"]

    if testing_value == "unbounded_mass":
        compute_unbounded_mass(s)

    if testing_value == "cum_unbounded_mass":
        print("adding cummulative unbounded nass")
        compute_cumulative_unbounded_mass(s, center)

    if testing_value not in s.data.keys():
        s.computeValueGas(testing_value)
        
    return testing_value


def compute_unbounded_mass(s):
    indgas = s.data['type'] == 0
    s.data['e'] = s.data['pot'][indgas] + 0.5 * ((s.data['vel'][indgas] ** 2).sum(axis=1))
    s.data["unbounded_mass"] = s.mass
    for i in np.where(s.data["e"] < 0)[0]:
        s.data["unbounded_mass"][i] = 0.0
    for i in np.where(s.type != 0):
        s.data["unbounded_mass"][i] = 0.0


def plot_profiles(output_dir, snapshot_name, plotting_dir, testing_value="rho", snapshot_number_array=[0, 8, 10],
                  center=False, log=True, new_fig=True, around_objects=False, around_density_peak=False,
                  line_profile=False, motion_axis=0, object_num=0, output_txt_files=False, relative_to_sink=False,
                  max_distance=None, central_id=None):
    if not os.path.exists(plotting_dir):
        os.mkdir(plotting_dir)

    if new_fig:
        set_new_fig_properties()

    line_colors = pylab.rcParams['axes.prop_cycle'].by_key()['color']
    labels = []
    suffix = ""
    for index, snapshot_number in enumerate(snapshot_number_array):
        if line_profile:
            p, s, suffix, testing_value = get_line_profile_for_snapshot(around_density_peak, around_objects, center,
                                                                          motion_axis, object_num, output_dir,
                                                                          snapshot_name, snapshot_number, testing_value,
                                                                        relative_to_sink, max_distance=max_distance,
                                                                        central_id=central_id)
            suffix += "_line_" + str(motion_axis)
        else:
            p, s, suffix, testing_value = get_radial_profile_for_snapshot(around_density_peak, around_objects, center,
                                                                      motion_axis, object_num, output_dir,
                                                                      snapshot_name, snapshot_number, testing_value,
                                                                          relative_to_sink, max_distance=max_distance,
                                                                          central_id=central_id)

        if output_txt_files:
            write_txt_file(p, plotting_dir, snapshot_number, s.time, suffix, testing_value)

        plot_to_figure(index, line_colors, log, p, s)
        labels.append("snap " + str(snapshot_number) + "," + str(round(s.time, 2)) + " [s]")
    if len(snapshot_number_array) > 1:
        legend(labels)
    filename = plotting_dir + "/" + testing_value + "_profile_" + suffix + "_".join([str(snap_num) for snap_num
                                                                                     in snapshot_number_array]) + ".png"
    print("saving to: ", filename)
    savefig(filename)
    print("saved fig")


def get_radial_profile_for_snapshot(around_density_peak, around_objects, center, motion_axis, object_num, output_dir,
                                    snapshot_name, snapshot_number, testing_value, relative_to_sink=False,
                                    max_distance=None, central_id=None):
    suffix = ""
    if around_objects:
        s, cell_indices, center, suffix = get_relevant_plotting_parameters_around_object(around_density_peak,
                                                                                         motion_axis,
                                                                                         object_num, output_dir,
                                                                                         snapshot_name,
                                                                                         snapshot_number)
        testing_value = compute_value(s, testing_value, center)

    else:
        s = gadget_readsnap(snapshot_number, output_dir, snapshot_name, loadonlytype=[0,1])
        testing_value = compute_value(s, testing_value, center)
        cell_indices = np.where(s.data['type'] == 0)

    center = get_center_array(s, center)
    if central_id is not None:
        center = s.pos[get_obj_index(s, central_id)][0]

    if max_distance is not None:
        relevant_cells = get_cells_inside_radius(s, center, max_distance)
        cell_indices = np.intersect1d(cell_indices, relevant_cells)

    if relative_to_sink:
        relevant_cells, sink_radius = get_cells_out_of_sink(s)
        cell_indices = np.intersect1d(cell_indices, relevant_cells)
        distance_right, val_right = plot_one_side(s, cell_indices, center, motion_axis, testing_value, right=True,
                                sink_size=s.parameters["SinkFormationRadius"])
        distance_left, val_left = plot_one_side(s, cell_indices, center, motion_axis, testing_value, right=False,
                               sink_size=s.parameters["SinkFormationRadius"])
        distance_left *= -1.0
        distances = np.concatenate((distance_left, distance_right))
        values = np.concatenate((val_left, val_right))
        sorted_ind = np.argsort(distances)
        p = np.row_stack((values[sorted_ind], distances[sorted_ind]))
    else:
        shells = 200
        if "cum_" in testing_value:
            shells = ceil(s.nparticlesall[0]/100)
        p = plot_snapshot_cells_around_center(cell_indices, center, s, testing_value,shells=shells)

    return p, s, suffix, testing_value


def get_center_array(s, center):
    if type(center) == list:
        center = pylab.array(center)
    elif type(center) != np.ndarray:
        center = s.center
    return center


def plot_one_side(s, cell_indices, center, motion_axis, testing_value, right=True, default_mode=-1,
                  sink_size=0):
    if right:
        half_cells = get_right_cells(s, cell_indices, center, motion_axis)
    else:
        half_cells = get_left_cells(s, cell_indices, center, motion_axis)
    #print("Duplicating ", len(half_cells), " positions and values")
    #pos_half = miror_positions(s, half_cells, center, motion_axis)
    #values_half, mode = mirror_values(s, half_cells, testing_value)
    distances = ((s.pos[half_cells] - center)**2).sum(axis=1)**0.5
    nshells = 200
    if sink_size > 0:
        nshells = max([nshells, ceil(distances.max()/(0.1*sink_size))])

    return get_averaged_data(distances, s.data[testing_value][half_cells], distances.max(), nshells)

    #if default_mode != -1:
    #    mode = default_mode
    #p = plot_profile_arrays_around_center(pos_half.astype('float64'), values_half, center, mode)

    #return p



def get_left_cells(s, cell_indices, center, motion_axis):
    left_cells = np.where(s.pos[:, motion_axis] < center[motion_axis])
    left_cells = np.intersect1d(left_cells, cell_indices)
    return left_cells


def get_right_cells(s, cell_indices, center, motion_axis):
    right_cells = np.where(s.pos[:, motion_axis] > center[motion_axis])
    right_cells = np.intersect1d(right_cells, cell_indices)
    return right_cells


def get_plotting_mode_and_array(s, cells_idx, testing_value):
    if testing_value == 'rho':
        values = s.data['mass'].astype('float64')[cells_idx]
        mode = 1
    else:
        values = s.data[testing_value].astype('float64')[cells_idx]
        mode = 2

    return values, mode

def miror_positions(s, cells, center, motion_axis):
    mirror_axis = np.array([1, 1, 1])
    mirror_axis[motion_axis] = -1
    positions = s.pos[cells] - center[motion_axis]
    positions = np.concatenate((positions[::-1] * mirror_axis, positions))
    positions += center[motion_axis]

    return positions

def mirror_values(s, cells, testing_value, right_to_left=True):
    half_values, mode = get_plotting_mode_and_array(s, cells, testing_value)
    half_values = np.concatenate((half_values[::-1], half_values))

    return half_values, mode

def plot_profile_arrays_around_center(distances, values, center, mode, shells=None):
    if shells is None:
        nshells = 200
    else:
        nshells = shells
    dr = 0
    if type(center) == list:
        center = pylab.array(center)
    p = calcGrid.calcRadialProfile(distances,
                                   values, mode, nshells, dr,
                                   center[0],
                                   center[1], center[2])
    return p

def plot_snapshot_cells_around_center(cell_indices, center, s, testing_value, shells=None):
    values, mode = get_plotting_mode_and_array(s, cell_indices, testing_value)
    if shells is None:
        shells = 200
    return plot_profile_arrays_around_center(s.data['pos'].astype('float64')[cell_indices],
                                   values, center, mode, shells=shells)


def get_line_profile_for_snapshot(around_density_peak, around_objects, center, motion_axis, object_num, output_dir,
                                    snapshot_name, snapshot_number, testing_value, relative_to_sink=False,
                                  max_distance=None, central_id=None):
    suffix = ""
    reference_r = 0
    if around_objects:
        s, cell_indices, center, suffix = get_relevant_plotting_parameters_around_object(around_density_peak,
                                                                                         motion_axis,
                                                                                         object_num, output_dir,
                                                                                         snapshot_name,
                                                                                         snapshot_number)
        testing_value = compute_value(s, testing_value, center)
    else:
        s = gadget_readsnap(snapshot_number, output_dir, snapshot_name, loadonlytype=[0])
        testing_value = compute_value(s, testing_value, center)
        center = get_center_array(s, center)
        cell_indices = np.where(s.data['mass'] != 0)
        if relative_to_sink:
            cell_indices, reference_r = get_cells_out_of_sink(s)

    if central_id is not None:
        center = s.pos[get_obj_index(s, central_id)][0]

    relevant_cells = np.where(
        (absolute(s.pos[:,(motion_axis + 1) % 3] - center[(motion_axis + 1) % 3]) < reference_r + 2 * s.data["vol"] ** (1.0 / 3)) &
        (absolute(s.pos[:,(motion_axis + 2) % 3] - center[(motion_axis + 2) % 3]) < reference_r + 2 * s.data["vol"] ** (1.0 / 3)))

    cell_indices = np.intersect1d(cell_indices, relevant_cells)
    if max_distance is not None:
        relevant_cells = get_cells_inside_radius(s, center, max_distance)
        cell_indices = np.intersect1d(cell_indices, relevant_cells)
    smoothed_pos_left, smoothed_val_left = \
        get_averaged_for_half_line(s, cell_indices, center, motion_axis, testing_value, right=False,
                                   max_size_shell=0.1*reference_r)
    smoothed_pos_left *= -1.0
    smoothed_pos_right, smoothed_val_right = \
        get_averaged_for_half_line(s, cell_indices, center, motion_axis, testing_value, right=True,
                                   max_size_shell=0.1*reference_r)

    smoothed_pos = np.concatenate((smoothed_pos_left, smoothed_pos_right))
    smoothed_val = np.concatenate((smoothed_val_left, smoothed_val_right))
    sorted_ind = np.argsort(smoothed_pos)
    p = np.row_stack((smoothed_val[sorted_ind], smoothed_pos[sorted_ind]))
    print(p.shape)

    return p, s, suffix, testing_value


def get_cells_out_of_sink(s):
    distance_from_sink = np.sqrt(((s.pos - s.center) ** 2).sum(axis=1))
    reference_r = s.parameters["SinkFormationRadius"]
    print("sink radius= ", s.parameters["SinkFormationRadius"])
    print("largest distance from sink: ", distance_from_sink.max())
    cell_indices = np.where((s.type ==0) & (distance_from_sink > (s.parameters["SinkFormationRadius"] + s.vol ** (1.0 / 3)))
                            & (s.mass != 0))
    print("largest x-distance from sink: ", (s.pos[cell_indices, 0] - s.center[0]).max())
    print("minimum x-distance from sink: ", (s.pos[cell_indices, 0] - s.center[0]).min())
    return cell_indices, reference_r

def get_cells_inside_radius(s, center, max_r):
    return np.where(((s.pos - center)**2).sum(axis=1)**0.5 < max_r)
def get_averaged_for_half_line(s, cell_indices, center, motion_axis, testing_value, right=True,
                               max_size_shell=0):
    if right:
        half_cells = get_right_cells(s, cell_indices, center, motion_axis)
    else:
        half_cells = get_left_cells(s, cell_indices, center, motion_axis)
    distances = absolute(s.data["pos"][half_cells, motion_axis] - center[motion_axis])
    values = s.data[testing_value][half_cells]
    if max_size_shell > 0:
        nshells = max([ceil(distances.max() / max_size_shell), 200])
    else:
        nshells = 200
    return get_averaged_data(distances, values, distances.max(), nshells=nshells)


def get_averaged_data(distances, values, size, nshells=200):
    dr = size / nshells
    print ("using dr= ", dr)
    smoothed_val = np.zeros(nshells)
    count_shells = np.zeros(nshells)
    smoothed_pos = np.zeros(nshells)
    for i in range(len(distances)):
        distance = distances[i]
        shell = floor(distance / dr)
        if shell < nshells:
            smoothed_val[shell] += values[i]
            count_shells[shell] += 1
    for shell in range(nshells):
        if count_shells[shell] > 0:
            smoothed_val[shell] /= count_shells[shell]
        smoothed_pos[shell] = (shell + 0.5) * dr
    relevant_i = np.where(smoothed_val != 0)

    return smoothed_pos[relevant_i], smoothed_val[relevant_i]


def plot_to_figure(index, line_colors, log, p, s):
    print("plotting")
    print("color= ", line_colors[index])
    if log:
        pylab.semilogy(p[1, :], p[0, :], color=line_colors[index])
    else:
        pylab.plot(p[1, :], p[0, :], color=line_colors[index])
    print("used color: ", line_colors[index], s.time)


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


def get_relevant_plotting_parameters_for_single(around_density_peak, output_dir, snapshot_name, snapshot_number, rho_cut=100):
    print("doing single object")
    s = gadget_readsnap(snapshot_number, output_dir, snapshot_name)
    if around_density_peak:
        print("calculating around density peak")
        center = s.pos[np.where(s.rho == s.rho.max())][-1, :]
    else:
        center = s.centerofmass()
    print("around center: ", center)
    i = np.where(s.rho > rho_cut)
    return s, i, center


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


def write_txt_file(p, plotting_dir, snapshot_number, snapshot_time, suffix, testing_value):
    txt_value_filename = plotting_dir + "/" + testing_value + "_profile_" + suffix + "_" + str(snapshot_number) \
                         + "_" + str(snapshot_time) + ".txt"
    print("creating txt file for ", testing_value, " with the name: ", txt_value_filename)
    with open(txt_value_filename, "w") as txt_file:
        txt_file.write("radius," + testing_value + "\r\n")
        for i in range(len(p[0])):
            txt_file.write(str(p[1, i]) + "," + str(p[0, i]) + "\r\n")
    print("txt file saved")


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--output_dir', type=str, help='diretory containing the snapshots', default="output")
    parser.add_argument('--snapshot_name', type=str, help='snapshots initilas', default="snapshot_")
    parser.add_argument('--plotting_dir', type=str, help='diretory to plot to', default="plots")
    parser.add_argument('--value', type=str, help='value to plot and compare profile', default="rho")
    parser.add_argument('--snapshot_nums', type=int, nargs='+', help='array of snapshot numbers to compare',
                        default=[0, 8, 10])
    parser.add_argument('--logplot', type=lambda x: (str(x).lower() in ['true', '1', 'yes']), help='logplot',
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
    parser.add_argument('--line_profile', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot only the value along the motion axis?',
                        default=False)
    parser.add_argument('--relative_to_sink', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot relative to sink position and radius, aassuming it is the last cell in the arrays',
                        default=False)
    parser.add_argument('--output_txt_files', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should also make txt files with plotting values?',
                        default=False)
    parser.add_argument('--max_distance', type=float,  help='radius around the object of interest to calculate for',
                        default=None)
    parser.add_argument('--central_id', type=int, help='id of the cell to centeralize around', default=None)

    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    if args.around_objects and not args.take_single_object:
        plot_profiles(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                      testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                      around_objects=args.around_objects, motion_axis=args.motion_axis,
                      around_density_peak=args.around_density_peak, line_profile=args.line_profile, object_num=1,
                      output_txt_files=args.output_txt_files, relative_to_sink=args.relative_to_sink,
                      max_distance=args.max_distance, central_id=args.central_id)
        plot_profiles(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                      testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                      around_objects=args.around_objects, motion_axis=args.motion_axis,
                      around_density_peak=args.around_density_peak, line_profile=args.line_profile, object_num=2,
                      output_txt_files=args.output_txt_files, new_fig=True,relative_to_sink=args.relative_to_sink,
                      max_distance=args.max_distance, central_id=args.central_id)
    else:
        plot_profiles(output_dir=args.output_dir, snapshot_name=args.snapshot_name, plotting_dir=args.plotting_dir,
                      testing_value=args.value, snapshot_number_array=args.snapshot_nums, log=args.logplot,
                      around_objects=args.around_objects, motion_axis=args.motion_axis,
                      around_density_peak=args.around_density_peak,  line_profile=args.line_profile,
                      output_txt_files=args.output_txt_files, relative_to_sink=args.relative_to_sink,
                      max_distance=args.max_distance, central_id=args.central_id)
