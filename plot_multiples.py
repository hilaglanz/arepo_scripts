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

    return ((snapshot.pos[obj_index] - center)**2).sum()**0.5

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

def plot_value_range(snapshot_list, snapshot_dir, plotting_dir, value, core_id=1e9+2, secondary_id=1e9,
                     tertiary_id=1e9+1, take_inner_mass=True):
    times = []
    values = []

    for snapshot_num in snapshot_list:
        snapshot = gadget_readsnap(snapshot_num, snapshot_dir)
        times.append(snapshot.time)
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
        plot_vs_time(value, values, times, False)
        filename = get_times_filename(snapshot_list, plotting_dir, value)
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
    parser.add_argument('--core_id', type=int, help='', default=1e9+2)
    parser.add_argument('--secondary_id', type=int, help='', default=1e9)
    parser.add_argument('--tertiary_id', type=int, help='', default=1e9+1)
    parser.add_argument('--take_inner_mass', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot according to all inner mass or just the core?',
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

    plot_value_range(snapshot_number_list, args.output_dir, args.plotting_dir, args.value, args.core_id,
                     args.secondary_id, args.tertiary_id, args.take_inner_mass)

