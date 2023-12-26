import pickle
import argparse
from make_plots import plot_single_value, regularize_time_units, basic_units, restore_basic_units, copy_current_units
from loadmodules import *

def plot_from_pickle(heatmap_pickle,stream_pickle, sink_pickle, vmin, vmax, cmap, logplot):
    with open(heatmap_pickle,'rb') as opened_heatmap_pickle:
        heatmap = pickle.load(opened_heatmap_pickle)
    with open(stream_pickle, 'rb') as opened_stream_pickle:
        stream = pickle.load(opened_stream_pickle)
    with open(sink_pickle, 'rb') as opened_scatter_pickle:
        sink = pickle.load(opened_scatter_pickle)

    pylab.pcolormesh(heatmap.pos_x, heatmap.pos_y, heatmap.slice, shading='flat',
                          norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax), rasterized=False,
                          cmap=cmap, logplot=logplot)

    streamplot(stream.pos_x, stream.pos_y, stream.vel_x, stream.vel_y, density=2, color='black')

    circ = Circle((sink.pos_x, sink.pos_y), sink.radius, fill=False, color='white', linestyle='dashed', linewidth=3.0)
    gca().add_patch(circ)


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--snapshot_file', type=str,  help='path to snapshot file', default= "output/snapshot_000.hdf5")
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "plots")
    parser.add_argument('--from_pickle', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should we plot from a saved pickle file? if not- plot and create pickle files',
                        default=False)
    parser.add_argument('--value', type=str,  help='value to be plotted', default= "rho")
    parser.add_argument('--cmap', type=str,  help='cmap for each subplot', default= "hot")
    parser.add_argument('--axes0', type=int,  help='horizonal axes to plot in', default= None)
    parser.add_argument('--axes1', type=int,  help='vertical axes to plot in', default= None)
    parser.add_argument('--vmin', type=float,  help='minimal range plotting', default=None)
    parser.add_argument('--vmax', type=float,  help='maximum range plotting', default=None)
    parser.add_argument('--boxsize', type=float,  nargs='+', help='boxsize', default=None)
    parser.add_argument('--logplot', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='logplot',
                        default=True)
    parser.add_argument('--res', type=int, help='plotting resolution', default=1024)
    parser.add_argument('--numthreads', type=int, help='threads for plotting', default=1)
    parser.add_argument('--center_x', type=float, help='point on x axis to be the center of the plot', default=16)
    parser.add_argument('--center_y', type=float, help='point on y axis to be the center of the plot', default=16)
    parser.add_argument('--center_z', type=float, help='point on z axis to be the center of the plot', default=16)
    parser.add_argument('--relative_to_sink_id', nargs='+', type=int,  help='id of sink particle to use as a reference point', default= None)
    parser.add_argument('--units_length', type=str,  help='name of the length units', default= "Ra")
    parser.add_argument('--units_velocity', type=str,  help='name of the velocity units default is cgs', default= "cs")

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    loaded_snap = gadget_readsnapname(args.snapshot_file)
    snap= args.snapshot_file.split('_')[-1].split('.')[0]
    print("doing snap ", snap)
    saving_file = args.saving_dir + "/pickled_" + snap
    val = args.value
    print(val)
    old_basic_units = copy_current_units()
    if args.from_pickle:
        plot_from_pickle(args.saving_file, args.saving_file+"_stream", args.saving_file+"_scatter", args.vmin,
                         args.vmax, args.cmap, args.logplot)
    else:
        plot_single_value(loaded_snap, value=val, cmap=args.cmap, box=[args.boxsize,args.boxsize,args.boxsize],
                          vrange=[args.vmin, args.vmax], logplot=args.logplot, res=args.res,
                          numthreads=args.numthreads, center=[args.center_x,args.center_y,args.center_z],
                          relative_to_sink_id=args.relative_to_sink_id,
                          plot_points=True,
                          unit_length=args.units_length,
                          unit_velocity=args.units_velocity, unit_density=None,
                          plot_velocities=True, axes=[args.axes0,args.axes1], saving_file=saving_file)

    regularize_time_units(loaded_snap)
    title('time : {:.2g}'.format(loaded_snap.time * basic_units["time"].factor) +
          " [" + basic_units["time"].unit + "]")
    restore_basic_units(old_basic_units)
    filename = args.saving_dir + "/Aslice_" + val + "_{0}.png".format(snap)
    print("saving to: ", filename)
    savefig(filename)
    print("saved fig")