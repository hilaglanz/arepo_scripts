import os
import glob
import argparse
import numpy as np
from loadmodules import *

name_and_units = {"rho":("density","g/cm^3"), "temp":("Temperature","K"), "vel":("Velocity","cm/s"), "mass":("Mass","g")}

def plot_single_value(loaded_snap, snap_num, value='rho', snapshotDir= "output", plottingDir="plots", firstSnap=0,lastSnap=-1,skipSteps=1,box=False,
               vrange=False,logplot=True, res=1024, numthreads=1, center=True,plot_points=True,
               additional_points_size=30,additional_points_shape='X', additional_points_color='w', units_length = 'cm',
               plot_velocities=False, newfig=True):
    label = value
    if value in name_and_units.keys():
        label = name_and_units[value][0]
        label += "[" + name_and_units[value][1] + "]"
    else:
        '''for val in name_and_units.keys():
            if value in val:
                label += "[" + name_and_units[value][1] + "]"
                break'''
    if "xnuc" in value:
        loaded_snap.data["rho"+value] = loaded_snap.rho * loaded_snap.data[value]
        value = "rho" + value
        label += " [" + name_and_units["rho"][1] + "]"
        print(value)

    print(value)
    loaded_snap.plot_Aslice(value, logplot=logplot, colorbar=True, cblabel=label, center=center, vrange=vrange,
                                  box=box, res=res, numthreads=numthreads, newfig=newfig)
    if box == False:
        box = [loaded_snap.boxsize, loaded_snap.boxsize]
    if plot_points:
        points, = np.where(loaded_snap.type > 0)
        if len(points) > 0:
            print("plotting points")
            for point in points:
                point_pos = loaded_snap.pos[point]

                scatter(point_pos[0], point_pos[1], additional_points_size, additional_points_color,
                        additional_points_shape)

                if loaded_snap.type[point] == 5:
                    print("plotting accretion radius of: ", loaded_snap.parameters['SinkFormationRadius'])
                    circ = Circle((point_pos[0], point_pos[1]), loaded_snap.parameters['SinkFormationRadius']
                                  , fill=False, color='white', linestyle='dashed', linewidth=3.0)
                    print(circ)
                    gca().add_patch(circ)
    if plot_velocities:
        loaded_snap.data['velx'] = loaded_snap.vel[:, 0]
        loaded_snap.data['vely'] = loaded_snap.vel[:, 1]
        slice_velx = loaded_snap.get_Aslice('velx', box=box, res=res, numthreads=numthreads)
        slice_vely = loaded_snap.get_Aslice('vely', box=box, res=res, numthreads=numthreads)
        posx = slice_velx['x'][:-1]
        posy = slice_velx['y'][:-1]
        velx = slice_velx['grid']
        vely = slice_vely['grid']
        streamplot(posx, posy, velx, vely, density=2, color='black')
        # quiver(loaded_snap.pos[:,0],loaded_snap.pos[:,1],loaded_snap.vel[:,0], loaded_snap.vel[:,1],
        # scale=50)#*loaded_snap.parameters['BoxSize']/box[0])

    xlabel('x [' + units_length + ']')
    ylabel('y [' + units_length + ']')

def get_single_value(value,index=0):
    if value is None:
        return value

    if len(value) > index:
        return value[index]

    return value[0]

def plot_range(value='rho', snapshotDir= "output", plottingDir="plots", firstSnap=0,lastSnap=-1,skipSteps=1,box=False,
               vrange=False,logplot=True, res=1024, numthreads=1, center=True,plot_points=True,
               additional_points_size=30,additional_points_shape='X', additional_points_color='w', units_length = 'cm',
               plot_velocities=False):
    snapshots = glob.glob(snapshotDir + '/./snapshot_*')
    print("found snapshots: ", snapshots)
    maxSnap=len(snapshots)
    if lastSnap == -1:
        lastSnap = maxSnap - 1
    else:
        lastSnap = min(lastSnap, maxSnap)
    print("doing from snapshot ", firstSnap, " to snapshot ", lastSnap)

    if firstSnap > lastSnap:
        print("cannot do firstSnap > lastSnap")
        return

    if not os.path.exists(plottingDir):
        os.mkdir(plottingDir)

    for snap in range(firstSnap,lastSnap+1,skipSteps):
        print("doing snapshot ",snap)
        loaded_snap = gadget_readsnap(snap, snapshotDir)
        print(type(value), value)
        if len(value) == 1:
            val = value[0]
            print(val)
            plot_single_value(loaded_snap, snap_num=snap, value=val, snapshotDir=snapshotDir, plottingDir=plottingDir,
                              firstSnap=firstSnap,
                              lastSnap=lastSnap, skipSteps=skipSteps, box=get_single_value(box),
                              vrange=get_single_value(vrange), logplot=get_single_value(logplot), res=res,
                              numthreads=numthreads, center=center, plot_points=plot_points,
                              additional_points_size=additional_points_size,
                              additional_points_shape=additional_points_shape,
                              additional_points_color=additional_points_color, units_length=units_length,
                              plot_velocities=plot_velocities)
            title('time : {:.2f} [s]'.format(loaded_snap.time))
            filename = plottingDir + "/Aslice_" + val + "_{0}.png".format(snap)
            print("saving to: ", filename)
            savefig(filename)
            print("saved fig")
        else:
            fig = figure(figsize=(36,20))
            fig.subplots_adjust(hspace=0.4,wspace=0.4)
            pyplot.rcParams.update({'font.size': 60, 'font.family': 'Serif'})
            num_figures = int(ceil(len(value)/2))
            for index,val in enumerate(value):
                curr_subplot = int(num_figures*100 + 21 + index)
                print("curr subplot: ", curr_subplot)
                subplot(curr_subplot)
                plot_single_value(loaded_snap, snap_num=snap, value=val, snapshotDir=snapshotDir, plottingDir=plottingDir,
                                  firstSnap=firstSnap,
                                  lastSnap=lastSnap, skipSteps=skipSteps, box=get_single_value(box,index),
                                  vrange=get_single_value(vrange,index), logplot=get_single_value(logplot,index),
                                  res=res,
                                  numthreads=numthreads, center=center, plot_points=plot_points,
                                  additional_points_size=additional_points_size,
                                  additional_points_shape=additional_points_shape,
                                  additional_points_color=additional_points_color, units_length=units_length,
                                  plot_velocities=plot_velocities, newfig=False)
                pyplot.rcParams.update({'font.size': 60, 'font.family': 'Serif'})
            #title('time : {:.2f} [s]'.format(loaded_snap.time))
            suptitle('time : {:.2f} [s]'.format(loaded_snap.time), fontsize='x-large')
            pyplot.rcParams.update({'font.size': 60, 'font.family': 'Serif'})
            filename = plottingDir + "/Aslice_" + "_".join(value) + "_{0}.png".format(snap)
            print("saving to: ", filename)
            savefig(filename)
            print("saved fig")


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--beginStep', type=int,  help='first step', default=0)
    parser.add_argument('--lastStep', type=int,  help='last step', default=-1)
    parser.add_argument('--skipStep', type=int, help='number of steps to skip', default=1)
    parser.add_argument('--source_dir', type=str,  help='path to snapshot files directory', default= sys.argv[0])
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "plots")
    parser.add_argument('--value', nargs='+', type=str,  help='value to be plotted', default= ["rho"])
    parser.add_argument('--vmin', type=float,  nargs='+', help='minimal range plotting', default=None)
    parser.add_argument('--vmax', type=float,  nargs='+', help='maximum range plotting', default=None)
    parser.add_argument('--boxsize', type=float,  nargs='+', help='boxsize', default=None)
    parser.add_argument('--logplot', nargs='+', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='logplot',
                        default=[True])
    parser.add_argument('--res', type=int, help='plotting resolution', default=1024)
    parser.add_argument('--numthreads', type=int, help='threads for plotting', default=1)
    parser.add_argument('--center_x', type=float, help='point on x axis to be the center of the plot', default=None)
    parser.add_argument('--center_y', type=float, help='point on y axis to be the center of the plot', default=None)
    parser.add_argument('--center_z', type=float, help='point on z axis to be the center of the plot', default=0)
    parser.add_argument('--plot_points', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='should plot other than gas?',
                        default=True)
    parser.add_argument('--plot_velocities', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='plot velocity field',
                        default=False)
    parser.add_argument('--additional_points_size', type=float,  help='point sizes in plot', default = 30)
    parser.add_argument('--additional_points_shape', type=str,  help='point shapes in plot', default= "X")
    parser.add_argument('--additional_points_color', type=str,  help='point colors in plot', default= "w")
    parser.add_argument('--units_length', type=str,  help='name of the length units', default= "cm")

    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    box = False
    if args.boxsize is not None:
        box = [[args.boxsize[i], args.boxsize[i]] for i in range(len(args.boxsize))]

    vrange = False
    if args.vmin is not None and args.vmax is not None:
        vrange = [[args.vmin[i], args.vmax[i]] for i in range(len(args.vmin))]

    center = False
    if args.center_x is not None and args.center_y is not None:
        center = [args.center_x, args.center_y,args.center_z]

    plot_range(args.value, args.source_dir, args.saving_dir, args.beginStep, args.lastStep, args.skipStep, box=box,
               vrange=vrange, logplot=args.logplot, res=args.res, numthreads= args.numthreads, center=center, plot_points=args.plot_points,
               additional_points_size=args.additional_points_size, additional_points_shape=args.additional_points_shape,
               additional_points_color=args.additional_points_color, units_length=args.units_length,
               plot_velocities=args.plot_velocities)
