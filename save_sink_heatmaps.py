import numpy as np
import pickle
from pickle import dump, load

from make_plots import plot_single_value, regularize_time_units, basic_units, restore_basic_units, copy_current_units
from loadmodules import *

class plotted_stream:
    pos_x = None
    pos_y = None
    vel_x = None
    vel_y = None
    def __init__(self, pos_x, pos_y, vel_x, vel_y):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.vel_x = vel_x
        self.vel_y = vel_y


class plotted_heatmap:
    pos_x = None
    pos_y = None
    slice = None
    def __init__(self, pos_x, pos_y, slice):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.slice = slice

class plotted_scatter:
    pos_x = None
    pos_y = None
    radius = None
    def __init__(self, pos_x, pos_y, radius):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.radius = radius


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
    parser.add_argument('--center_x', type=float, help='point on x axis to be the center of the plot', default=None)
    parser.add_argument('--center_y', type=float, help='point on y axis to be the center of the plot', default=None)
    parser.add_argument('--center_z', type=float, help='point on z axis to be the center of the plot', default=0)
    parser.add_argument('--relative_to_sink_id', nargs='+', type=int,  help='id of sink particle to use as a reference point', default= None)

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    loaded_snap = gadget_readsnapname(args.snapshot_file)
    snap=int(args.snapshot_file.split('_')[-1].split('.')[0])
    print("doing snap ", snap)
    saving_file = args.saving_dir + "/pickled_" + snap
    val = args.value
    print(val)
    old_basic_units = copy_current_units()
    plot_single_value(loaded_snap, value=val, cmap=args.cmap, box=args.boxsize,
                      vrange=[args.vmin, args.vmax], logplot=args.logplot, res=args.res,
                      numthreads=args.numthreads, center=args.center,
                      relative_to_sink_id=args.relative_to_sink_id,
                      plot_points=True,
                      unit_length="Ra",
                      unit_velocity=None, unit_density=None,
                      plot_velocities=True, axes=[args.axes0,args.axes1], saving_file=saving_file)

    regularize_time_units(loaded_snap)
    title('time : {:.2g}'.format(loaded_snap.time * basic_units["time"].factor) +
          " [" + basic_units["time"].unit + "]")
    restore_basic_units(old_basic_units)
    filename = args.saving_dir + "/Aslice_" + val + "_{0}.png".format(snap)
    print("saving to: ", filename)
    savefig(filename)
    print("saved fig")