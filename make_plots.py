import os
import glob
import argparse
import numpy as np
from loadmodules import *

name_and_units = {"rho":[r'$\rho$',r'$g/cm^3$', 1.0], "temp":["Temperature","K", 1.0], "vel":["Velocity","$cm/s$", 1.0],
                  "mass":["Mass","g", 1.0], "time":["time", "s", 1.0], "length": ["length", "cm",1.0],
                  "acc":["acceleration", "$cm/s^2$", 1.0], "pres":["Pressure", "Ba", 1.0], "entr":["Entropy","", 1.0]}

species = ['n', 'p', '^{4}He', '^{11}B', '^{12}C', '^{13}C', '^{13}N', '^{14}N', '^{15}N', '^{15}O',
           '^{16}O', '^{17}O', '^{18}F', '^{19}Ne', '^{20}Ne', '^{21}Ne', '^{22}Ne', '^{22}Na',
           '^{23}Na', '^{23}Mg', '^{24}Mg', '^{25}Mg', '^{26}Mg', '^{25}Al', '^{26}Al',
           '^{27}Al', '^{28}Si', '^{29}Si', '^{30}Si', '^{29}P', '^{30}P', '^{31}P', '^{31}S',
           '^{32}S', '^{33}S', '^{33}Cl', '^{34}Cl', '^{35}Cl', '^{36}Ar', '^{37}Ar', '^{38}Ar',
           '^{39}Ar', '^{39}K', '^{40}Ca', '^{43}Sc', '^{44}Ti', '^{47}V', '^{48}Cr', '^{51}Mn',
           '^{52}Fe', '^{56}Fe', '^{55}Co', '^{56}Ni', '^{58}Ni', '^{59}Ni']

def add_name_and_unit(value, name, unit, factor = 1.0):
    name_and_units[value] = [name, unit, factor]
def change_snap_units(loaded_snap):
    loaded_snap.time *= name_and_units["time"][2]
    loaded_snap.data["pos"][:, 0] *= name_and_units["length"][2]
    loaded_snap.data["pos"][:, 1] *= name_and_units["length"][2]
    loaded_snap.data["pos"][:, 2] *= name_and_units["length"][2]

    loaded_snap.center *= name_and_units["length"][2]
    loaded_snap.data["vel"][:, 0] *= name_and_units["vel"][2]
    loaded_snap.data["vel"][:, 1] *= name_and_units["vel"][2]
    loaded_snap.data["vel"][:, 2] *= name_and_units["vel"][2]

    loaded_snap.data["mass"] *= name_and_units["mass"][2]
    loaded_snap.data["rho"] *= name_and_units["rho"][2]
    loaded_snap.data["vol"] *= (name_and_units["length"][2] ** 3)
    # TODO: convert also temperature
def project_vector(v,r):
    dist = np.sqrt((r*r).sum(axis=1))
    return ((r*v).sum(axis=1)) / dist

def plot_stream(loaded_snap, value='vel', xlab='x', ylab='y', axes=[0,1], box=False, res=1024, numthreads=1):
    loaded_snap.data[value + xlab] = loaded_snap.data[value][:, axes[0]]
    loaded_snap.data[value + ylab] = loaded_snap.data[value][:, axes[1]]
    slice_velx = loaded_snap.get_Aslice(value + xlab, box=box, res=res, center=center, numthreads=numthreads)
    slice_vely = loaded_snap.get_Aslice(value + ylab, box=box, res=res, center=center, numthreads=numthreads)
    posx = slice_velx['x'][:-1]
    posy = slice_velx['y'][:-1]
    velx = pylab.transpose(slice_velx['grid'])
    vely = pylab.transpose(slice_vely['grid'])
    streamplot(posx, posy, velx, vely, density=2, color='black')
    # quiver(loaded_snap.pos[:,0],loaded_snap.pos[:,1],loaded_snap.vel[:,0], loaded_snap.vel[:,1],
    # scale=50)#*loaded_snap.parameters['BoxSize']/box[0])

def plot_single_value(loaded_snap, value='rho', cmap="hot", box=False, vrange=False,logplot=True, res=1024, numthreads=1,
                      relative_to_sink_id=None, center=True,plot_points=True, additional_points_size=30,
                      additional_points_shape='X', additional_points_color='w', unit_length='cm', unit_velocity="$cm/s$",
                      unit_density=r'$g/cm^3$', plot_velocities=False, plot_bfld=False,
                      newfig=True, axes=[0,1], modified_units = False, ignore_types=[]):
    label = value
    convert_to_cgs = True
    if unit_velocity is not None:
        name_and_units["vel"][1] = unit_velocity
        name_and_units["time"][1] = r'$' + unit_length + "/" + unit_velocity + '$'
        convert_to_cgs = False

    if unit_density is not None:
        name_and_units["rho"][1] = unit_density
        convert_to_cgs = False

    if not modified_units:
        if convert_to_cgs:
            name_and_units["vel"][2] *= float(loaded_snap.parameters["UnitVelocity_in_cm_per_s"])
            name_and_units["length"][2] *= float(loaded_snap.parameters["UnitLength_in_cm"])
            name_and_units["acc"][2] *= float(loaded_snap.parameters["UnitVelocity_in_cm_per_s"])**2 / float(loaded_snap.parameters["UnitLength_in_cm"])
            name_and_units["time"][2] *= float(loaded_snap.parameters["UnitLength_in_cm"]) / float(loaded_snap.parameters["UnitVelocity_in_cm_per_s"])
            name_and_units["mass"][2] *= float(loaded_snap.parameters["UnitMass_in_g"])
            name_and_units["rho"][2] *= float(loaded_snap.parameters["UnitMass_in_g"]) / float(loaded_snap.parameters["UnitLength_in_cm"])**3
            print("converting to cgs units")

    change_snap_units(loaded_snap)
    modified_units = True

    print("units: ")
    for val in name_and_units.values():
        print(val[0], val[2])

    loaded_snap, value = calculate_label_and_value(loaded_snap, value, relative_to_sink_id)

    if value in name_and_units.keys():
        label = name_and_units[value][0]
        label += " [" + name_and_units[value][1] + "]"

    print(value)
    xlab = chr(ord('x') + axes[0])
    ylab = chr(ord('x') + axes[1])

    loaded_snap.plot_Aslice(value, logplot=logplot, colorbar=True, cblabel=label, cmap=cmap, center=center, vrange=vrange,
                                  box=box, res=res, numthreads=numthreads, newfig=newfig, axes=axes,
                            minimum=min(1e-8, 0.1*vrange[0]))
    if box == False:
        box = [loaded_snap.boxsize, loaded_snap.boxsize]
    if plot_points:
        points, = np.where((loaded_snap.type > 0) & (loaded_snap.type not in ignore_types))
        if len(points) > 0:
            print("plotting points")
            for point in points:
                point_pos = loaded_snap.data["pos"][point]

                scatter(point_pos[axes[0]], point_pos[axes[1]], additional_points_size, additional_points_color,
                        additional_points_shape)

                if loaded_snap.type[point] == 5:
                    print("plotting accretion radius of: ", loaded_snap.parameters['SinkFormationRadius']*name_and_units["length"][2])
                    circ = Circle((point_pos[axes[0]], point_pos[axes[1]]), loaded_snap.parameters['SinkFormationRadius']*name_and_units["length"][2]
                                  , fill=False, color='white', linestyle='dashed', linewidth=3.0)
                    print(circ)
                    gca().add_patch(circ)
    if plot_velocities:
        plot_stream(loaded_snap, value='vel', xlab=xlab, ylab=ylab, axes=axes, box=box, res=res, numthreads=numthreads)
    elif plot_bfld:
        plot_stream(loaded_snap, value='bfld', xlab=xlab, ylab=ylab, axes=axes, box=box, res=res, numthreads=numthreads)

    xlabel(xlab + ' [' + unit_length + ']', loc="left")
    ylabel(ylab + ' [' + unit_length + ']')


def calculate_label_and_value(loaded_snap, value, relative_to_sink_id):
    if value in loaded_snap.data.keys():
        if len(loaded_snap.data[value].shape) == 1:
            add_computed_value_to_name_and_unit_dict(loaded_snap, value)
            return loaded_snap, value

    if loaded_snap.computeValueGas(value) == 0:
        if len(loaded_snap.data[value].shape) == 1:
            add_computed_value_to_name_and_unit_dict(loaded_snap, value)
            return loaded_snap, value

    if "xnuc" in value:
        loaded_snap.data["rho" + value] = loaded_snap.rho * loaded_snap.data[value]
        value = "rho" + value
        add_name_and_unit(value, r'$\rho \left(' + species[int(value.split("xnuc")[-1])] + r"\right)$",
                          name_and_units["rho"][1])

    if value == "mean_a":
        loaded_snap.calculate_mean_a()
        add_name_and_unit(value, "Mean Atomic Weight", "")

    if value == "bfld" or value == "B":
        loaded_snap.data["B"] = np.sqrt((loaded_snap.data['bfld'] * loaded_snap.data['bfld']).sum(axis=1))
        value = "B"

    if "vel" in value:
        loaded_snap.data['velx'] = loaded_snap.data["vel"][:, 0]
        loaded_snap.data['vely'] = loaded_snap.data["vel"][:, 1]
        loaded_snap.data['velz'] = loaded_snap.data["vel"][:, 2]

    if value == "vel":
        loaded_snap.data['vel_size'] = np.sqrt((loaded_snap.vel ** 2).sum(axis=1))
        value = "vel_size"
        add_name_and_unit(value, "Velocity", name_and_units["vel"][1], name_and_units["vel"][2])

    if "vort" in value:
        loaded_snap.data['vortx'] = loaded_snap.data["vort"][:, 0]
        loaded_snap.data['vorty'] = loaded_snap.data["vort"][:, 1]
        loaded_snap.data['vortz'] = loaded_snap.data["vort"][:, 2]

    if value == "vort":
        loaded_snap.data['vort_size'] = np.sqrt((loaded_snap.vort ** 2).sum(axis=1))
        value = "vort_size"

    if value == "mach":
        loaded_snap.computeMach()
        add_name_and_unit(value, "Mach number", name_and_units["vel"][1], name_and_units["vel"][2])

    if value == "cs" or "sound" in value:
        loaded_snap.computeMach()
        value = "sound"
        add_name_and_unit(value, "$c_s$", name_and_units["vel"][1], name_and_units["vel"][2])

    if "grap" in value:
        loaded_snap.data['grapx'] = loaded_snap.data["grap"][:, 0]
        loaded_snap.data['grapy'] = loaded_snap.data["grap"][:, 1]
        loaded_snap.data['grapz'] = loaded_snap.data["grap"][:, 2]

    if value == "grap":
        if relative_to_sink_id is None:
            loaded_snap.data['grap_size'] = np.sqrt((loaded_snap.grap ** 2).sum(axis=1))
            value = "grap_size"
        else:
            dist, r, sink_idk = calculate_sink_properties(loaded_snap, relative_to_sink_id)
            loaded_snap.data['grap_r'] = project_vector(loaded_snap.data["grap"], r)
            value = "grap_r"

    if value == "g_sink":
        dist, r, sink_idk = calculate_sink_properties(loaded_snap, relative_to_sink_id)
        loaded_snap.data['g_sink'] = G * loaded_snap.mass[sink_idk] / dist ** 2
        add_name_and_unit(value, "g_sink", name_and_units["acc"][1], name_and_units["acc"][2])

    if value == "HSE" and relative_to_sink_id is not None:
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "g_sink", relative_to_sink_id)
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "grap", relative_to_sink_id)
        loaded_snap.data["HSE"] = -1.0 * loaded_snap.data["grap_r"] / (loaded_snap.data["g_sink"] * loaded_snap.rho)

    if "grav" in value:
        loaded_snap.data['gravx'] = loaded_snap.data["grav"][:, 0]
        loaded_snap.data['gravy'] = loaded_snap.data["grav"][:, 1]
        loaded_snap.data['gravz'] = loaded_snap.data["grav"][:, 2]

    if value == "grav":
        if relative_to_sink_id is None:
            loaded_snap.data['grav_size'] = np.sqrt((loaded_snap.grav ** 2).sum(axis=1))
            value = "grav_size"
        else:
            dist, r, sink_idk = calculate_sink_properties(loaded_snap, relative_to_sink_id)
            loaded_snap.data['grav_r'] = project_vector(loaded_snap.data["grav"], r)
            value = "grav_r"

    if value == "grap_r_over_rho" and relative_to_sink_id is not None:
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "grap_r", relative_to_sink_id)
        loaded_snap.data[value] = -1.0 * loaded_snap.data["grap_r"] / loaded_snap.rho[np.where(loaded_snap.type == 0)]

    if value == "v_grav":
        loaded_snap.data[value+"_x"] = (loaded_snap.vel[np.where(loaded_snap.type == 0)] * loaded_snap.grav[:,[0,3,6]]).sum(axis=1)
        loaded_snap.data[value+"_y"] = (loaded_snap.vel[np.where(loaded_snap.type == 0)] * loaded_snap.grav[:,[1,4,7]]).sum(axis=1)
        loaded_snap.data[value+"_z"] = (loaded_snap.vel[np.where(loaded_snap.type == 0)] * loaded_snap.grav[:,[2,5,8]]).sum(axis=1)
        loaded_snap.data[value] = np.array([[loaded_snap.data[value+"_x"][i], loaded_snap.data[value+"_y"][i],
                                    loaded_snap.data[value+"_z"][i]] for i in range(loaded_snap.nparticlesall[0])])
        if relative_to_sink_id is not None:
            dist, r, sink_idk = calculate_sink_properties(loaded_snap, relative_to_sink_id)
            value = "v_grav_r"
            loaded_snap.data[value] = -1.0 * project_vector(loaded_snap.data["v_grav"], r)

    if value == "momentum_vdot":
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "g_sink", relative_to_sink_id)
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "grap_r_over_rho", relative_to_sink_id)
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "v_grav", relative_to_sink_id)
        loaded_snap.data[value] = loaded_snap.data["g_sink"] + loaded_snap.data["grap_r_over_rho"] + loaded_snap.data["v_grav_r"]
        add_name_and_unit(value, r"$g_{sink} - \nabla P /\rho - v\cdot \nabla v$", name_and_units["acc"][1])

    return loaded_snap, value


def add_computed_value_to_name_and_unit_dict(loaded_snap, value):
    if value not in name_and_units.keys():
        name = value
        if value in loaded_snap.hdf5_name_conversion.keys():
            name = loaded_snap.hdf5_name_conversion.keys()
        add_name_and_unit(value, name, "")


def calculate_sink_properties(loaded_snap, relative_to_sink_id):
    sink_idk = get_sink_idk(loaded_snap, relative_to_sink_id)
    sink_pos = loaded_snap.pos[sink_idk]
    r = loaded_snap.pos[np.where(loaded_snap.type == 0)] - sink_pos
    dist = np.sqrt((r * r).sum(axis=1))
    return dist, r, sink_idk


def get_sink_idk(loaded_snap, relative_to_sink_id):
    sink_idks = np.where(loaded_snap.type == 5)
    sink_idk = sink_idks[0][relative_to_sink_id]
    print("doing for sink particle")
    return sink_idk


def get_single_value(value,index=0):
    if value is None:
        return value

    if len(value) > index:
        return value[index]

    return value[0]

def get_snapshot_number_list(snapshotDir="outupt", snapshotName="snapshot_", firstSnap=0, lastSnap=-1, skipSteps=1):
    snapshots = glob.glob(snapshotDir + '/./' + snapshotName + '*')
    snapshots.sort(key=lambda st: int(st.split(snapshotName)[1].split(".")[0]))
    sorted(snapshots)
    print("found snapshots: ", snapshots[0], " to ", snapshots[-1])
    maxSnap = int((snapshots[-1].split('snapshot_')[-1]).split('.hdf5')[0])
    if lastSnap == -1:
        lastSnap = maxSnap
    else:
        lastSnap = min(lastSnap, maxSnap)
    print("doing from snapshot ", firstSnap, " to snapshot ", lastSnap)

    if firstSnap > lastSnap:
        print("cannot do firstSnap > lastSnap")
        return []

    return range(firstSnap, lastSnap+1, skipSteps)

def plot_range(value=['rho'], snapshotDir= "output", plottingDir="plots", firstSnap=0,lastSnap=-1,skipSteps=1,box=False,
               vrange=False, cmap=["hot"], logplot=True, res=1024, numthreads=1, center=True, relative_to_sink_id=None,
               plot_points=True,
               additional_points_size=30,additional_points_shape='X', additional_points_color='w', units_length = 'cm',
               units_velocity="$cm/s$", units_density=r'$g/cm^3$', plot_velocities=False, plot_bfld=False,
               axes_array=[[0,1]], ignore_types=[]):

    if not os.path.exists(plottingDir):
        os.mkdir(plottingDir)
    convert_to_cgs = False
    if units_velocity is None and units_density is None:
        convert_to_cgs = True
    modified_units = False
    curr_cmap = cmap[0]
    for snap in get_snapshot_number_list(snapshotDir, "snapshot_", firstSnap, lastSnap, skipSteps):
        print("doing snapshot ",snap)
        loaded_snap = gadget_readsnap(snap, snapshotDir)
        if len(value) == 1:
            val = value[0]
            print(val)
            plot_single_value(loaded_snap, value=val, cmap=curr_cmap, box=get_single_value(box),
                              vrange=get_single_value(vrange), logplot=get_single_value(logplot), res=res,
                              numthreads=numthreads, center=center, relative_to_sink_id=relative_to_sink_id,
                              plot_points=plot_points,
                              additional_points_size=additional_points_size,
                              additional_points_shape=additional_points_shape,
                              additional_points_color=additional_points_color, unit_length=units_length,
                              unit_velocity= units_velocity, unit_density= units_density,
                              plot_velocities=plot_velocities, plot_bfld= plot_bfld, axes=get_single_value(axes_array),
                              modified_units=modified_units, ignore_types=ignore_types)
            title('time : {:.2g}'.format(loaded_snap.time) + " [" + name_and_units["time"][1] + "]")
            filename = plottingDir + "/Aslice_" + val + "_{0}.png".format(snap)
            print("saving to: ", filename)
            savefig(filename)
            print("saved fig")
        else:
            fig = figure(figsize=(36,20))
            fig.subplots_adjust(hspace=0.4,wspace=0.4)
            rcParams.update({'font.size': 40, 'font.family': 'Serif'})
            rcParams['text.usetex'] = True
            num_figures = int(ceil(len(value)/2))
            for index,val in enumerate(value):
                if num_figures >= 1:
                    curr_subplot = int(num_figures*100 + 21 + index)
                print("curr subplot: ", curr_subplot)
                subplot(curr_subplot)
                curr_cmap = cmap[index % len(cmap)]
                plot_single_value(loaded_snap,  value=val, cmap=curr_cmap, box=get_single_value(box,index),
                                  vrange=get_single_value(vrange,index), logplot=get_single_value(logplot,index),
                                  res=res,
                                  numthreads=numthreads, center=center, relative_to_sink_id=relative_to_sink_id,
                                  plot_points=plot_points,
                                  additional_points_size=additional_points_size,
                                  additional_points_shape=additional_points_shape,
                                  additional_points_color=additional_points_color, unit_length=units_length,
                                  unit_velocity= units_velocity, unit_density= units_density,
                                  plot_velocities=plot_velocities, plot_bfld= plot_bfld, newfig=False,
                                  axes=get_single_value(axes_array, index), ignore_types=ignore_types)
                rcParams.update({'font.size': 40, 'font.family': 'Serif'})
                rcParams['text.usetex'] = True

            #title('time : {:.2f} [s]'.format(loaded_snap.time))
            suptitle('time : {:.2g}'.format(loaded_snap.time) + " [" + name_and_units["time"][1] + "]", fontsize='x-large')
            rcParams.update({'font.size': 40, 'font.family': 'Serif'})
            rcParams['text.usetex'] = True
            filename = plottingDir + "/Aslice_" + "_".join(value) + "_{0}.png".format(snap)
            print("saving to: ", filename)
            savefig(filename)
            print("saved fig")
        close('all')
        modified_units = True


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--beginStep', type=int,  help='first step', default=0)
    parser.add_argument('--lastStep', type=int,  help='last step', default=-1)
    parser.add_argument('--skipStep', type=int, help='number of steps to skip', default=1)
    parser.add_argument('--source_dir', type=str,  help='path to snapshot files directory', default= sys.argv[0])
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "plots")
    parser.add_argument('--value', nargs='+', type=str,  help='value to be plotted', default= ["rho"])
    parser.add_argument('--cmap', nargs='+', type=str,  help='cmap for each subplot', default= ["hot"])
    parser.add_argument('--axes0', nargs='+', type=int,  help='horizonal axes to plot in', default= None)
    parser.add_argument('--axes1', nargs='+', type=int,  help='vertical axes to plot in', default= None)
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
    parser.add_argument('--relative_to_sink_id', nargs='+', type=int,  help='id of sink particle ro use as a reference point', default= None)
    parser.add_argument('--plot_points', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='should plot other than gas?',
                        default=True)
    parser.add_argument('--ignore_types', nargs='+', type=int,  help='particle types to ignore', default=[])
    parser.add_argument('--plot_velocities', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='plot velocity field',
                        default=False)
    parser.add_argument('--plot_bfld', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),  help='plot magnetic field stream',
                        default=False)
    parser.add_argument('--additional_points_size', type=float,  help='point sizes in plot', default = 30)
    parser.add_argument('--additional_points_shape', type=str,  help='point shapes in plot', default= "X")
    parser.add_argument('--additional_points_color', type=str,  help='point colors in plot', default= "w")
    parser.add_argument('--units_length', type=str,  help='name of the length units', default= "cm")
    parser.add_argument('--units_velocity', type=str,  help='name of the velocity units default is cgs', default= None)
    parser.add_argument('--units_density', type=str,  help='name of the density units default is cgs', default= None)
    parser.add_argument('--factor_mass', type=float,  help='multiply mass unit by this factor', default=1.0)
    parser.add_argument('--factor_length', type=float,  help='multiply length unit by this factor', default=1.0)
    parser.add_argument('--factor_velocity', type=float,  help='multiply velocity unit by this factor', default=1.0)

    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    box = False
    if args.boxsize is not None:
        box = [[args.boxsize[i], args.boxsize[i], args.boxsize[i]] for i in range(len(args.boxsize))]

    vrange = False
    if args.vmin is not None and args.vmax is not None:
        vrange = [[args.vmin[i], args.vmax[i]] for i in range(len(args.vmin))]

    center = False
    if args.center_x is not None and args.center_y is not None:
        center = [args.center_x, args.center_y,args.center_z]

    axes_array = [[0,1]]
    if args.axes0 is not None and args.axes1 is not None:
        axes_array = [[args.axes0[i],args.axes1[i]] for i in range(len(args.axes0))]

    name_and_units["rho"][2] *= (args.factor_mass/args.factor_length**3)
    name_and_units["length"][2] *= args.factor_length
    name_and_units["vel"][2] *= (args.factor_velocity)
    name_and_units["mass"][2] *= (args.factor_mass)
    name_and_units["time"][2] *= (args.factor_length/args.factor_velocity)
    #TODO: add conversion to temperature

    plot_range(args.value, args.source_dir, args.saving_dir, args.beginStep, args.lastStep, args.skipStep, box=box,
               vrange=vrange, logplot=args.logplot, cmap=args.cmap, res=args.res, numthreads= args.numthreads, center=center,
               relative_to_sink_id=args.relative_to_sink_id,
               plot_points=args.plot_points, additional_points_size=args.additional_points_size,
               additional_points_shape=args.additional_points_shape,
               additional_points_color=args.additional_points_color,
               units_length=args.units_length, units_velocity=args.units_velocity, units_density= args.units_density,
               plot_velocities=args.plot_velocities, plot_bfld= args.plot_bfld, axes_array=axes_array,
               ignore_types=args.ignore_types)
