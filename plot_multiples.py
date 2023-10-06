from loadmodules import *
from time_plots import *

def get_obj_index(snapshot, obj_id):
    obj_index = np.where(snapshot.id == obj_id)

    return obj_index

def get_separation(snapshot, obj_id, center, take_inner_mass=False):
    obj_index = get_obj_index(snapshot, obj_id)
    if take_inner_mass:
        inner_cells = np.where(((snapshot.pos-center)**2).sum(axis=1) < ((snapshot.pos[obj_index]-center)**2).sum())
        inner_mass = snapshot.mass[inner_cells].sum()
        inner_pos = (snapshot.pos[inner_cells] * snapshot.mass[inner_cells][:,None]).sum(axis=0)/inner_mass
        center = inner_pos

    return ((snapshot.pos[obj_index] - center)**2).sum()**0.5 / rsol

def get_velocity(snapshot, obj_id, center, center_obj_id, take_inner_mass=False):
    obj_index = get_obj_index(snapshot, obj_id)
    if take_inner_mass:
        inner_cells = np.where(
            ((snapshot.pos - center) ** 2).sum(axis=1) < ((snapshot.pos[obj_index] - center) ** 2).sum())
        inner_mass = snapshot.mass[inner_cells].sum()
        inner_vel = (snapshot.vel[inner_cells] * snapshot.mass[inner_cells][:, None]).sum(axis=0) / inner_mass
    else:
        inner_vel = snapshot.vel[np.where(snapshot.id == center_obj_id)]

    return ((snapshot.vel[obj_index] - inner_vel) ** 2).sum() ** 0.5

def get_surrounding_rho(snapshot, obj_id, size):
    obj_index = get_obj_index(snapshot, obj_id)
    surrounding_cells = np.where((snapshot.type == 0 ) &
                                 (((snapshot.pos - snapshot.pos[obj_index]) ** 2).sum(axis=1)**0.5 < size))

    return snapshot.rho[surrounding_cells].mean()

def get_surrounding_value(snapshot, obj_id, size, value):
    if value == "rho":
        return get_surrounding_rho(snapshot, obj_id, size)

    if value == "unbounded_mass_frac":
        indgas = snapshot.data['type'] == 0
        snapshot.data["unbounded_mass"] = snapshot.data['pot'][indgas] + 0.5 * ((snapshot.data['vel'][indgas] ** 2).sum(axis=1))

    if value not in snapshot.data.keys():
        if snapshot.computeValueGas(value) != 0:
            print("could not add ", value, " for gas properties")
            return

    obj_index = get_obj_index(snapshot, obj_id)
    surrounding_cells = np.where((snapshot.type == 0) &
                                 (((snapshot.pos - snapshot.pos[obj_index]) ** 2).sum(axis=1)**0.5 < size))
    if value == "unbounded_mass_frac":
        return snapshot.data["unbounded_mass"][surrounding_cells].sum()

    com_value = (snapshot.data[value][surrounding_cells] * snapshot.mass[surrounding_cells][:, None]).sum(axis=0) / \
                snapshot.mass[surrounding_cells].sum()
    if len(com_value) > 1:
        return (com_value**2).sum()**0.5

    return com_value

def calculate_separation_vector(snapshot, obj_id, center, center_obj_id, take_inner_mass=False):
    obj_index = get_obj_index(snapshot, obj_id)
    if take_inner_mass:
        inner_cells = np.where(
            ((snapshot.pos - center) ** 2).sum(axis=1) < ((snapshot.pos[obj_index] - center) ** 2).sum())
        center = (snapshot.pos[inner_cells]**2).sum(axis=0)**0.5
    else:
        center = snapshot.pos[np.where(snapshot.id == center_obj_id)][0]

    return snapshot.pos[obj_index][0] - center

def get_drag(snapshot, obj_id, center, center_obj_id, take_inner_mass=False):
    obj_index = get_obj_index(snapshot, obj_id)
    separation = calculate_separation_vector(snapshot, obj_id, center, center_obj_id, take_inner_mass)
    snapshot = calculate_value_relative_to_vector(snapshot, 'acce', separation)

    return snapshot.data["acce_u"][obj_index]  * snapshot.mass[obj_index]

def get_enclosed_mass(snapshot, obj_id, center, center_obj_id=None):
    obj_index = get_obj_index(snapshot, obj_id)
    if center_obj_id is not None:
        central_index = get_obj_index(snapshot, center_obj_id)
        center = snapshot.pos[central_index]
    inner_cells = np.where(
        ((snapshot.pos - center) ** 2).sum(axis=1) < ((snapshot.pos[obj_index] - center) ** 2).sum())

    return snapshot.mass[inner_cells].sum()


def get_out_mass(snapshot, obj_id, center, center_obj_id=None):
    obj_index = get_obj_index(snapshot, obj_id)
    if center_obj_id is not None:
        central_index = get_obj_index(snapshot, center_obj_id)
        center = snapshot.pos[central_index]
    inner_cells = np.where(
        ((snapshot.pos - center) ** 2).sum(axis=1) > ((snapshot.pos[obj_index] - center) ** 2).sum())

    return snapshot.mass[inner_cells].sum()


def plot_value_range(snapshot_list, snapshot_dir, plotting_dir, value, core_id=1e9+1, secondary_id=1e9,
                     tertiary_id=1e9+2, take_inner_mass=True, surrounding_radius=10*rsol, around_object_id=1e9 + 2):
    times = []
    values = []
    ylab = value
    suffix = ""

    if "drag" in value or "enclosed" in value or "out_mass" in value:
        suffix += "_id" + str(around_object_id - 1e9)
    if "surrounding" in value:
        suffix += "_id" + str(around_object_id - 1e9) + "_" + str(surrounding_radius/rsol) + "rsol"
    elif take_inner_mass:
        suffix += "_inner_cells"

    if "separation" in value:
        ylab += " [" + r'$R_\odot$' + "]"
    elif "vel" in value:
        ylab += " [cm/s]"
    if "mass" in value:
        ylab += " [" + r'$M_\odot$' + "]"

    unbounded_mass_prev = 0
    time_prev = 0
    for snapshot_num in snapshot_list:
        snapshot = gadget_readsnap(snapshot_num, snapshot_dir)
        times.append(snapshot.time / day)
        if value == "inner_separation":
            values.append(get_separation(snapshot, secondary_id, snapshot.pos[get_obj_index(snapshot, core_id)],
                                         take_inner_mass=take_inner_mass))
        elif value == "outer_separation":
            values.append(get_separation(snapshot, tertiary_id, snapshot.pos[get_obj_index(snapshot, core_id)],
                                         take_inner_mass=take_inner_mass))
        elif value == "inner_velocity":
            values.append(get_velocity(snapshot, secondary_id, snapshot.pos[get_obj_index(snapshot, core_id)], core_id,
                                         take_inner_mass=take_inner_mass))
        elif value == "outer_velocity":
            values.append(get_velocity(snapshot, tertiary_id, snapshot.pos[get_obj_index(snapshot, core_id)], core_id,
                                         take_inner_mass=take_inner_mass))
        elif "surrounding" in value:
            values.append(get_surrounding_value(snapshot, around_object_id, surrounding_radius, value.split('_')[-1]))

        elif value == "drag":
            values.append(get_drag(snapshot, around_object_id, center=snapshot.pos[get_obj_index(snapshot, core_id)],
                                  center_obj_id=core_id, take_inner_mass=take_inner_mass))
        elif "unbounded_mass" in value and not "surrounding" in value:
            unbounded_mass = snapshot.computeUnboundMass()
            if "dot" in value:
                if time_prev == 0:
                    time_prev = float(snapshot.parameters['TimeBegin'])
                values.append((unbounded_mass - unbounded_mass_prev) / (snapshot.time - time_prev) / msol)
                time_prev = snapshot.time
                unbounded_mass_prev = unbounded_mass
            else:
                values.append(unbounded_mass / msol)
        elif "enclosed_mass" in value:
            values.append(get_enclosed_mass(snapshot, around_object_id, center=snapshot.pos[get_obj_index(snapshot, core_id)],
                                  center_obj_id=core_id) / msol)
        elif "out_mass" in value:
            values.append(get_out_mass(snapshot, around_object_id, center=snapshot.pos[get_obj_index(snapshot, core_id)],
                                  center_obj_id=core_id) / msol)




    plot_vs_time(value, values, times, False)
    filename = get_times_filename(snapshot_list, plotting_dir, value, suffix)
    xlabel("Time [days]")
    ylabel(ylab)
    savefig(filename)
    save_txt_value_vs_time(filename, values, times)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--output_dir', type=str, help='diretory containing the snapshots', default="output")
    parser.add_argument('--snapshot_name', type=str, help='snapshots initilas', default="snapshot_")
    parser.add_argument('--beginStep', type=int,  help='first step', default=0)
    parser.add_argument('--lastStep', type=int,  help='last step', default=-1)
    parser.add_argument('--skipStep', type=int, help='number of steps to skip', default=1)
    parser.add_argument('--plotting_dir', type=str, help='diretory to plot to', default="plots")
    parser.add_argument('--value', type=str, help='value to plot and compare profile', default="inner_separation")
    parser.add_argument('--logplot', type=lambda x: (str(x).lower() in ['true', '1', 'yes']), help='logplot',
                        default=True)
    parser.add_argument('--core_id', type=int, help='', default=1e9+1)
    parser.add_argument('--secondary_id', type=int, help='', default=1e9)
    parser.add_argument('--tertiary_id', type=int, help='', default=1e9+2)
    parser.add_argument('--take_inner_mass', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot according to all inner mass or just the core?',
                        default=False)
    parser.add_argument('--surrounding_radius', type=float,  help='radius around the object of interest to calculate for',
                        default=10*rsol)
    parser.add_argument('--around_object_id', type=int, help='id of the object to plot surrounding of', default=1e9+2)

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

    plot_value_range(snapshot_number_list, args.output_dir, args.plotting_dir, args.value, args.core_id,
                     args.secondary_id, args.tertiary_id, args.take_inner_mass, args.surrounding_radius,
                     args.around_object_id)

