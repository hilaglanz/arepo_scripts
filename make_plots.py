import os
import glob
import argparse
import pickle

from loadmodules import *
from save_sink_heatmaps import save_heatmap,save_scatter,save_stream


species = ['n', 'p', '^{4}He', '^{11}B', '^{12}C', '^{13}C', '^{13}N', '^{14}N', '^{15}N', '^{15}O',
           '^{16}O', '^{17}O', '^{18}F', '^{19}Ne', '^{20}Ne', '^{21}Ne', '^{22}Ne', '^{22}Na',
           '^{23}Na', '^{23}Mg', '^{24}Mg', '^{25}Mg', '^{26}Mg', '^{25}Al', '^{26}Al',
           '^{27}Al', '^{28}Si', '^{29}Si', '^{30}Si', '^{29}P', '^{30}P', '^{31}P', '^{31}S',
           '^{32}S', '^{33}S', '^{33}Cl', '^{34}Cl', '^{35}Cl', '^{36}Ar', '^{37}Ar', '^{38}Ar',
           '^{39}Ar', '^{39}K', '^{40}Ca', '^{43}Sc', '^{44}Ti', '^{47}V', '^{48}Cr', '^{51}Mn',
           '^{52}Fe', '^{56}Fe', '^{55}Co', '^{56}Ni', '^{58}Ni', '^{59}Ni']

class UnitConversion:
    def __init__(self, unit, factor=1.0):
        self.unit = unit
        self.factor = factor

class UnitName:
    def __init__(self, name, unit):
        self.name = name
        self.unit_name = unit

basic_units = {"rho":UnitConversion(r'$g/cm^3$'), "temp":UnitConversion("K"), "vel":UnitConversion("$cm/s$"),
               "mass":UnitConversion("g"), "time":UnitConversion("s"), "length": UnitConversion("cm"),
               "vol": UnitConversion(r"$cm^3$"), "acce":UnitConversion("$cm/s^2$"), "pres":UnitConversion("Ba"),
               "force":UnitConversion("dyne"),
               "u":UnitConversion("erg"), "ang_mom": UnitConversion(r"$cm^2 / s$") ,"none": UnitConversion("")}

name_and_units = {"rho":UnitName(r'$\rho$',"rho"), "temp":UnitName("Temperature","temp"),
                  "vel":UnitName("Velocity","vel"), "sound":UnitName("c_s","vel"),
                  "mach":UnitName(r'$\mathcal{M}$', "none"), "mass":UnitName("Mass","mass"),
                  "time":UnitName("Time", "time"), "length": UnitName("Length", "length"),
                  "pos":UnitName("Position", "length"), "vol": UnitName("Volume","vol"),
                  "acce":UnitName("Acceleration", "acce"), "pres": UnitName("Pressure", "pres"),
                  "force": UnitName("Force", "force"),
                  "u": UnitName("Energy","u"), "entr": UnitName("Entropy", "none")}
def add_name_and_unit(value, name, unit):
    if value not in name_and_units.keys():
        if unit not in basic_units.keys():
            basic_units[unit] = UnitConversion(unit)
        name_and_units[value] = UnitName(name, unit)

def change_snap_units(loaded_snap):
    for key in loaded_snap.data.keys():
        add_computed_value_to_name_and_unit_dict(loaded_snap, key)
        if name_and_units[key].unit_name == "none":
            continue
        print("changing units of ", key, " by a factor of ", basic_units[name_and_units[key].unit_name].factor)
        loaded_snap.data[key] *= basic_units[name_and_units[key].unit_name].factor

def copy_basic_units(old_dictionary, new_dictionary=None):
    if new_dictionary is None:
        new_dictionary = {}
    print("replacing existing dictionary")
    for key in basic_units.keys():
        new_dictionary[key] = old_dictionary[key]

    return new_dictionary
def copy_current_units():
    return copy_basic_units(basic_units)
def restore_basic_units(old_basic_units):
    copy_basic_units(old_basic_units, new_dictionary=basic_units)


def regularize_length_units(boxsize):
    if basic_units["length"].unit != "cm":
        return

    au = 1.495978707*10**13
    if basic_units["length"].factor <= 1.0/parsec or boxsize > 2 * parsec:
        basic_units["length"].factor /= parsec
        basic_units["length"].unit = r'$pc'
        print("changed length units to pc")

    elif basic_units["length"].factor <= 1.0/au or boxsize > 2 * au:
        basic_units["length"].factor /= au
        basic_units["length"].unit = r'$AU$'
        print("changed length units to AU")

    elif basic_units["length"].factor <= 1.0/rsol or boxsize > 2 * rsol:
        basic_units["length"].factor /= rsol
        basic_units["length"].unit = r'$R_\odot$'
        print("changed length units to Rsun")
def regularize_time_units(snapshot):
    if basic_units["time"].unit != "s":
        return

    if basic_units["time"].factor <= 1/yr or snapshot.time * basic_units["time"].factor > yr:
        basic_units["time"].factor /= yr
        basic_units["time"].unit = "years"
    elif basic_units["time"].factor <= 1/day or snapshot.time * basic_units["time"].factor > day:
        basic_units["time"].factor /= day
        basic_units["time"].unit = "days"
    elif basic_units["time"].factor <= 1/hour or snapshot.time * basic_units["time"].factor > hour:
        basic_units["time"].factor /= hour
        basic_units["time"].unit = "hours"


def change_unit_conversion(factor_length, factor_velocity, factor_mass):
    basic_units["rho"].factor *= (factor_mass / (factor_length ** 3))
    basic_units["length"].factor *= factor_length
    basic_units["vol"].factor *= factor_length ** 3
    basic_units["vel"].factor *= factor_velocity
    basic_units["ang_mom"].factor *= factor_velocity * factor_length
    basic_units["acce"].factor *= ((factor_velocity ** 2) / factor_length)
    basic_units["mass"].factor *= factor_mass
    basic_units["time"].factor *= (factor_length/factor_velocity)
    basic_units["pres"].factor *= ((factor_mass * factor_velocity ** 2) / (factor_length ** 3)) # mass*acc/area
    basic_units["force"].factor *= ((factor_mass * factor_velocity ** 2) / factor_length) # mass*acc
    #basic_units["entr"].factor = (basic_units["pres"].factor / (basic_units["rho"].factor ** (5.0 / 3)))
    basic_units["u"].factor *= (factor_mass * (factor_velocity ** 2))
    # TODO: convert also temperature?
def project_vector(v,r):
    dist = np.sqrt((r*r).sum(axis=1))
    return ((r*v).sum(axis=1)) / dist


def plot_stream(loaded_snap, value='vel', xlab='x', ylab='y', axes=[0,1], box=False, res=1024, numthreads=1, center=None,
                saving_file=None):
    if center is None:
        center = loaded_snap.center

    loaded_snap.data[value + xlab] = loaded_snap.data[value][:, axes[0]]
    loaded_snap.data[value + ylab] = loaded_snap.data[value][:, axes[1]]
    slice_velx = loaded_snap.get_Aslice(value + xlab, box=box, res=res, center=center, numthreads=numthreads)
    slice_vely = loaded_snap.get_Aslice(value + ylab, box=box, res=res, center=center, numthreads=numthreads)
    posx = slice_velx['x'][:-1]
    posy = slice_velx['y'][:-1]
    velx = pylab.transpose(slice_velx['grid'])
    vely = pylab.transpose(slice_vely['grid'])
    streamplot(posx, posy, velx, vely, density=2, color='black')
    save_stream(posx, posy, velx, vely, saving_file)
    # quiver(loaded_snap.pos[:,0],loaded_snap.pos[:,1],loaded_snap.vel[:,0], loaded_snap.vel[:,1],
    # scale=50)#*loaded_snap.parameters['BoxSize']/box[0])





def plot_single_value(loaded_snap, value='rho', cmap="hot", box=False, vrange=False,logplot=True, res=1024, numthreads=1,
                      relative_to_sink_id=None, central_id=None, center=True,plot_points=True, additional_points_size=30,
                      additional_points_shape='X', additional_points_color='w', unit_length='cm', unit_velocity="$cm/s$",
                      unit_density=r'$g/cm^3$', plot_velocities=False, plot_bfld=False,
                      newfig=True, axes=[0,1], modified_units = False, ignore_types=[], colorbar=True,
                      plot_xlabel=True, plot_ylabel=True, factor_value=1.0, units_value=None, saving_file=None, contour=False):

    if box == False:
        box = [loaded_snap.boxsize, loaded_snap.boxsize]

    change_basic_units(loaded_snap, unit_density, unit_length, unit_velocity)
    change_value_units(value, units_value, factor_value)
    change_snap_units(loaded_snap)

    print("units: ")
    for val in name_and_units.values():
        print(val.name, basic_units[val.unit_name].factor)
    loaded_snap, value = calculate_label_and_value(loaded_snap, value, relative_to_sink_id, central_id=central_id,)

    print(value)

    label = extract_label(value)
    xlab = chr(ord('x') + axes[0])
    ylab = chr(ord('x') + axes[1])

    if central_id is not None:
        center = loaded_snap.pos[np.where(loaded_snap.id == central_id)][0]
        print("centralizing around id ", central_id, " at ", center)

    if contour:
        levels=10
        proj=False
    else:
        levels= [0.99]
        proj=False
    loaded_snap.plot_Aslice(value, logplot=logplot, colorbar=colorbar, cblabel=label, cmap=cmap, center=center, vrange=vrange,
                                  box=box, res=res, numthreads=numthreads, newfig=newfig, axes=axes,
                            minimum=min(1e-8, 0.1*vrange[0]), contour=contour, levels=levels, proj=proj)
    stream_saving_file = None
    if saving_file is not None:
        stream_saving_file = saving_file + "_stream"
    save_heatmap(axes, box, center, loaded_snap, numthreads, res, saving_file, value)

    if plot_velocities:
        plot_stream(loaded_snap, value='vel', xlab=xlab, ylab=ylab, axes=axes, box=box, res=res, numthreads=numthreads,
                    center=center, saving_file=stream_saving_file)
    elif plot_bfld:
        plot_stream(loaded_snap, value='bfld', xlab=xlab, ylab=ylab, axes=axes, box=box, res=res, numthreads=numthreads,
                    center=center)
    if plot_points:
        points = [idx for idx in range(loaded_snap.npart) if (loaded_snap.type[idx] not in ignore_types) and
                  (loaded_snap.type[idx] > 0)]
        if len(points) > 0:
            print("plotting points")
            for point in points:
                point_pos = loaded_snap.data["pos"][point]

                scatter(point_pos[axes[0]], point_pos[axes[1]], additional_points_size, additional_points_color,
                        additional_points_shape)

                if loaded_snap.type[point] == 5:
                    print("plotting accretion radius of: ",
                          loaded_snap.parameters['SinkFormationRadius']*basic_units["length"].factor)
                    circ = Circle((point_pos[axes[0]], point_pos[axes[1]]),
                                  loaded_snap.parameters['SinkFormationRadius']*basic_units["length"].factor
                                  , fill=False, color='white', linestyle='dashed', linewidth=3.0)
                    print(circ)
                    gca().add_patch(circ)

                    save_scatter(axes, loaded_snap, point_pos, basic_units["length"].factor, saving_file)
    '''
    regularize_length_units(max(box))
    change_ticks(xaxis=True)
    change_ticks(xaxis=False)
    '''
    if plot_xlabel:
        xlabel(xlab + ' [' + basic_units["length"].unit + ']', loc="left")
    if plot_ylabel:
        ylabel(ylab + ' [' + basic_units["length"].unit + ']')




def change_ticks(xaxis=True):
    ticklabels = []
    if xaxis:
        ticks = xticks()
    else:
        ticks = yticks()

    for tick in ticks:
        if (tick == 0):
            ticklabels += [r'$0.0$']
        else:
            ticklabels += [r'$%.2f \cdot 10^{%d}$' % (tick * basic_units["length"].factor /
                                                       10 ** (ceil(log10(abs(tick * basic_units["length"].factor)))),
                                                       ceil(log10(abs(tick * basic_units["length"].factor))))]
    if xaxis:
        xticks(ticks, ticklabels, size=24, y=-0.1, va='baseline')
    else:
        yticks(ticks, ticklabels, size=24, ha='right')


def extract_label(value):
    label = value
    if value in name_and_units.keys():
        label = name_and_units[value].name
        if name_and_units[value].unit_name != "none":
            label += " [" + basic_units[name_and_units[value].unit_name].unit + "]"
    return label


def change_basic_units(loaded_snap, unit_density, unit_length, unit_velocity):
    convert_to_cgs = True
    if unit_velocity is not None:
        basic_units["vel"].unit = unit_velocity
        basic_units["acce"].unit = r'$' + unit_velocity + '^2 /' + unit_length + '$'
        basic_units["time"].unit = r'$' + unit_length + "/" + unit_velocity + '$'
        basic_units["ang_mom"].unit = r'$' + unit_length + "\cdot " + unit_velocity + '$'
        basic_units["vol"].unit = r'$' + unit_length + '^3$'
        basic_units["length"].unit = r'$' + unit_length + '$'

        if unit_density is not None:
            basic_units["pres"].unit = r'$' + unit_density + '\cdot' + unit_velocity + '^2$'
        else:
            basic_units["pres"].unit = r'$\rho_\infty \cdot ' + unit_velocity + '^2$'

        convert_to_cgs = False
    if unit_density is not None:
        basic_units["rho"].unit = unit_density
        convert_to_cgs = False
    if convert_to_cgs:
        change_unit_conversion(factor_length=float(loaded_snap.parameters["UnitLength_in_cm"]),
                               factor_velocity=float(loaded_snap.parameters["UnitVelocity_in_cm_per_s"]),
                               factor_mass=float(loaded_snap.parameters["UnitMass_in_g"]))
        print("converting to cgs units")


def change_value_units(value, units_value, factor_value):
    if factor_value is None:
        return

    if value in basic_units.keys():
        basic_units[value].factor = factor_value
        print("changing ", value, " factor to ", factor_value)
        if units_value is not None and units_value != "None":
            basic_units[value].unit = units_value
            print("changing ", value, " units to ", units_value)
    elif "size" in value and value.split("_size")[0] in basic_units.keys():
        basic_units[value.split("_size")[0]].factor = factor_value
        print("changing ", value.split("_size")[0], " factor to ", factor_value)
        if units_value is not None and units_value != "None":
            basic_units[value.split("_size")[0]].unit = units_value
            print("changing ", value.split("_size")[0], " units to ", units_value)


def get_value_at_inf(value, data):
    if value not in data:
        print("data does not have ", value)
        return None
    border_ids = np.where(data['pos'][:,0] < 1)
    n = 2
    while len(border_ids[0]) < 1:
        border_ids = np.where(data['pos'][:,0] < n)
        n += 1

    val_inf = data[value][border_ids].mean()
    print("calculated ", value, " at infinity= ", val_inf)

    return val_inf

def calculate_label_and_value(loaded_snap, value, relative_to_sink_id, central_id=None):
    if "xnuc" in value:
        loaded_snap.data["rho" + value] = loaded_snap.rho * loaded_snap.data[value]
        value = "rho" + value
        add_name_and_unit(value, r'$\rho \left(' + species[int(value.split("xnuc")[-1])] + r"\right)$", "rho")

    if value in loaded_snap.data.keys():
        if len(loaded_snap.data[value].shape) == 1:
            add_computed_value_to_name_and_unit_dict(loaded_snap, value)
            return loaded_snap, value

    if loaded_snap.computeValueGas(value) == 0:
        if len(loaded_snap.data[value].shape) == 1:
            add_computed_value_to_name_and_unit_dict(loaded_snap, value)
            return loaded_snap, value

    if value == "mean_a":
        loaded_snap.calculate_mean_a()
        add_name_and_unit(value, "Mean Atomic Weight", "none")

    if value == "bfld" or value == "B":
        loaded_snap.data["B"] = np.sqrt((loaded_snap.data['bfld'] * loaded_snap.data['bfld']).sum(axis=1))
        value = "B"

    if "vel" in value:
        loaded_snap.data['vel_x'] = loaded_snap.data["vel"][:, 0]
        loaded_snap.data['vel_y'] = loaded_snap.data["vel"][:, 1]
        loaded_snap.data['vel_z'] = loaded_snap.data["vel"][:, 2]
        add_name_and_unit(value, "Velocity" + value.split('_')[-1], "vel")

    if "drag" in value:
        loaded_snap.data[value+'_x'] = loaded_snap.data["acce"] * loaded_snap.mass[:, None]
        loaded_snap.data[value+'_y'] = loaded_snap.data["acce"] * loaded_snap.mass[:, None]
        loaded_snap.data[value+'_z'] = loaded_snap.data["acce"] * loaded_snap.mass[:, None]
        loaded_snap.data["drag"] = loaded_snap.data["acce"] * loaded_snap.mass[:, None]
        if "drag" == value:
            value += "_size"
        add_name_and_unit(value, "Drag Force" + value.split('_')[-1], "force")

    if "_size" in value:
        loaded_snap.data[value] = np.sqrt((loaded_snap.data[value.split('_size')[0]] ** 2).sum(axis=1))

    if "vort" in value:
        loaded_snap.data['vort_x'] = loaded_snap.data["vort"][:, 0]
        loaded_snap.data['vort_y'] = loaded_snap.data["vort"][:, 1]
        loaded_snap.data['vort_z'] = loaded_snap.data["vort"][:, 2]

    if value == "mach":
        loaded_snap.computeMach()

    if value == "cs" or "sound" in value:
        loaded_snap.computeMach()
        value = "sound"


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
        add_name_and_unit(value, "g_sink", "acce")

    if value == "grap_r_over_rho" and relative_to_sink_id is not None:
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "grap", relative_to_sink_id)
        loaded_snap.data[value] = -1.0 * loaded_snap.data["grap_r"] / loaded_snap.rho
        add_name_and_unit(value, r"$\nabla P /\rho$", "acce")

    if value == "HSE" and relative_to_sink_id is not None:
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "grap_r_over_rho", relative_to_sink_id)
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "g_sink", relative_to_sink_id)
        loaded_snap.data["HSE"] = loaded_snap.data["grap_r_over_rho"] / loaded_snap.data["g_sink"]

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

    if value == "g-grap_r_over_rho":
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "grap_r_over_rho", relative_to_sink_id)
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "g_sink", relative_to_sink_id)
        loaded_snap.data[value] = loaded_snap.data["g_sink"] - loaded_snap.data["grap_r_over_rho"]
        add_name_and_unit(value, r"$g_{sink} - \nabla P /\rho$", "acce")

    if value == "momentum_vdot":
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "g-grap_r_over_rho", relative_to_sink_id)
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "v_grav", relative_to_sink_id)
        loaded_snap.data[value] = loaded_snap.data["g-grap_r_over_rho"] - loaded_snap.data["v_grav_r"]
        add_name_and_unit(value, r"$g_{sink} - \nabla P /\rho - v\cdot \nabla v$", "acce")

    if value == "entr_ratio":
        loaded_snap, temp_value = calculate_label_and_value(loaded_snap, "entr", relative_to_sink_id)
        entr_inf = get_value_at_inf("entr", loaded_snap.data)
        loaded_snap.data[value] = (loaded_snap.data["entr"] - entr_inf) / entr_inf
        add_name_and_unit(value, r"$\delta s / s_\infty$", "none")

    if value == "angular_momentum":
        ind = loaded_snap.data['type'] == 0
        if relative_to_sink_id is not None:
            dist, r, sink_idk = calculate_sink_properties(loaded_snap, relative_to_sink_id)
            rcm = loaded_snap.data['pos'][sink_idk,:]
            vcm = loaded_snap.data['vel'][sink_idk,:]
        else:
            rcm = loaded_snap.center
            vcm = np.zeros((1,3))
        print("calculating angular momentum from ", rcm, vcm)
        print((loaded_snap.pos[ind].astype(np.float64)[:, 0] - rcm[0]))
        print((loaded_snap.vel[ind].astype(np.float64)[:, 1] - vcm[0]))
        print((loaded_snap.pos[ind].astype(np.float64)[:, 0] - rcm[0])
                * (loaded_snap.vel[ind].astype(np.float64)[:, 1] - vcm[0]))
        loaded_snap.data[value] = (loaded_snap.mass[ind].astype(np.float64) * (
                (loaded_snap.pos[ind].astype(np.float64)[:, 0] - rcm[0])
                * (loaded_snap.vel[ind].astype(np.float64)[:, 1] - vcm[0]) -
                (loaded_snap.pos[ind].astype(np.float64)[:, 1] - rcm[1])
                * (loaded_snap.vel[ind].astype(np.float64)[:, 0] - vcm[1]))
                ).sum()
        add_name_and_unit(value, "Angular momentum", "ang_mom")

    return loaded_snap, value

def calculate_value_relative_to_vector(loaded_snap, value, vector):
    if value not in loaded_snap.data.keys() or len(loaded_snap.data[value].shape) <= 1:
        print("cannot compute ", value, " relative to the motion axis")
        return
    else:
        print("computing ", value, "relative to ", vector)
        vector_size = np.sqrt((vector ** 2).sum())
        loaded_snap.data[value+"_v"] = (loaded_snap.data[value] * vector).sum(axis=1) / vector_size
        loaded_snap.data[value+"_u"] = np.sqrt((loaded_snap.data[value] ** 2).sum(axis=1) -
                                               loaded_snap.data[value+"_v"]**2)
        add_name_and_unit(value+"_v", name_and_units[value].name + "_v", name_and_units[value].unit_name)
        add_name_and_unit(value+"_u", name_and_units[value].name + "_u", name_and_units[value].unit_name)

        return loaded_snap
def add_computed_value_to_name_and_unit_dict(loaded_snap, value):
    if value not in name_and_units.keys():
        name = value
        if value in loaded_snap.hdf5_name_conversion.keys():
            name = loaded_snap.hdf5_name_conversion.keys()
        add_name_and_unit(value, name, "none")
        print("added ", value, " into the name and units dictionary with no units")


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
    snapshots = glob.glob(snapshotDir + '/./' + snapshotName + '*.hdf5')
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

def plot_single_value_evolutions(value=['rho'], snapshotDir= "output", plottingDir="plots", firstSnap=0,lastSnap=-1,skipSteps=1,box=False,
               vrange=False, cmap=["hot"], logplot=True, res=1024, numthreads=1, center=True, relative_to_sink_id=None,
               central_id=None, plot_points=True,
               additional_points_size=30,additional_points_shape='X', additional_points_color='w', units_length = 'cm',
               units_velocity="$cm/s$", units_density=r'$g/cm^3$', plot_velocities=False, plot_bfld=False,
               axes_array=[[0,1]], ignore_types=[], horizontal=True, relative_to_motion=False, snapshots_list=None):
    if not os.path.exists(plottingDir):
        os.mkdir(plottingDir)
    convert_to_cgs = False
    if units_velocity is None and units_density is None:
        convert_to_cgs = True
    modified_units = False
    if snapshots_list is None:
        snapshots_list = get_snapshot_number_list(snapshotDir, "snapshot_", firstSnap, lastSnap, skipSteps)
    num_figures = len(snapshots_list)

    for index, val in enumerate(value):
        print(val)
        fig = figure(figsize=(num_figures*15, 17))
        rcParams.update({'font.size': 70, 'font.family': 'Serif', 'axes.formatter.useoffset':False})
        rcParams['text.usetex'] = True
        curr_cmap = cmap[index % len(cmap)]
        for snap_i, snap in enumerate(snapshots_list):
            print("doing snapshot ", snap)
            if horizontal:
                curr_subplot = int(100 + 10*num_figures + (snap_i+1))
                ax = subplot(1, num_figures, 1)
                curr_ax = subplot(curr_subplot, sharey=ax)
            else:
                curr_subplot = int(num_figures * 100 + 10 + (snap_i + 1))
                ax = subplot(num_figures, 1, 1)
                curr_ax = subplot(curr_subplot, sharex=ax)
            loaded_snap = gadget_readsnap(snap, snapshotDir,
                                          loadonlytype=[t for t in range(6) if t not in ignore_types])

            old_basic_units = copy_current_units()
            print("curr snapshot: ", snap_i + 1)
            plot_single_value(loaded_snap,  value=val, cmap=curr_cmap, box=get_single_value(box,index),
                                  vrange=get_single_value(vrange,index), logplot=get_single_value(logplot,index),
                                  res=res,
                                  numthreads=numthreads, center=get_single_value(center,index),
                              relative_to_sink_id=get_single_value(relative_to_sink_id,index),
                                  central_id=get_single_value(central_id, index), plot_points=plot_points,
                                  additional_points_size=additional_points_size,
                                  additional_points_shape=additional_points_shape,
                                  additional_points_color=additional_points_color, unit_length=units_length,
                                  unit_velocity= units_velocity, unit_density= units_density,
                                  plot_velocities=plot_velocities, plot_bfld= plot_bfld, newfig=False,
                                  axes=get_single_value(axes_array, index), ignore_types=ignore_types, colorbar=False,
                              plot_xlabel=(horizontal is True or ((not horizontal) and (snap_i == num_figures))),
                              plot_ylabel=(not horizontal or ((horizontal) and (snap_i == 0))))
            #subplot(curr_subplot)
            regularize_time_units(loaded_snap)
            ax.tick_params(axis='x',labelrotation=45)
            curr_ax.set_title('{:.3g}'.format(loaded_snap.time * basic_units["time"].factor) +
                              " [" + basic_units["time"].unit + "]", fontsize='70',loc='right')
            restore_basic_units(old_basic_units)

            rcParams.update({'font.size': 70, 'font.family': 'Serif', 'axes.formatter.useoffset':False})
            rcParams['text.usetex'] = True
            if horizontal is True and snap_i!=0:
                curr_ax.set_axis_off()
            if horizontal is False and snap_i+1!=num_figures:
                curr_ax.set_axis_off()
        if horizontal:
            fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.002, hspace=0.2)
        else:
            fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.2, hspace=0.002)
        cax = fig.add_axes([0.905, 0.1946, 0.04/num_figures, 0.702])
        if "xnuc" in val:
            val = "rho" + val
        colorbar(cax=cax, label= name_and_units[val].name + " [" + basic_units[name_and_units[val].unit_name].unit + "]",
                 aspect=15, pad=0, shrink=1)
        tight_layout(pad=0, h_pad=0, w_pad=0, rect=(0.01, 0, 0.9, 1))
        #title('time : {:.2f} [s]'.format(loaded_snap.time))
        rcParams.update({'font.size': 70, 'font.family': 'Serif', 'axes.formatter.useoffset':False})
        rcParams['text.usetex'] = True
        filename = plottingDir + "/Aslice_" + val + "_" + "_".join([str(s) for s in snapshots_list]) + ".png".format(snap)
        print("saving to: ", filename)
        savefig(filename)
        print("saved fig")
        close('all')
        modified_units = True

def plot_range(value=['rho'], snapshotDir= "output", plottingDir="plots", firstSnap=0,lastSnap=-1,skipSteps=1,box=False,
               vrange=False, cmap=["hot"], logplot=True, res=1024, numthreads=1, center=True, relative_to_sink_id=None,
               central_id=None, plot_points=True,
               additional_points_size=30,additional_points_shape='X', additional_points_color='w', units_length = 'cm',
               units_velocity="$cm/s$", units_density=r'$g/cm^3$', plot_velocities=False, plot_bfld=False,
               axes_array=[[0,1]], ignore_types=[], per_value_evolution=False, relative_to_motion=False,
               factor_value=[1.0], units_value=[None], contour=False, snapshots_list=None):

    if per_value_evolution:
        return plot_single_value_evolutions(value, snapshotDir, plottingDir, firstSnap, lastSnap, skipSteps, box,
               vrange, cmap, logplot, res, numthreads, center, relative_to_sink_id, central_id,
               plot_points,
               additional_points_size,additional_points_shape, additional_points_color, units_length,
               units_velocity, units_density, plot_velocities, plot_bfld,
               axes_array, ignore_types, snapshots_list=snapshots_list)

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
            old_basic_units = copy_current_units()
            plot_single_value(loaded_snap, value=val, cmap=curr_cmap, box=get_single_value(box),
                              vrange=get_single_value(vrange), logplot=get_single_value(logplot), res=res,
                              numthreads=numthreads, center=get_single_value(center),
                              relative_to_sink_id=get_single_value(relative_to_sink_id),
                              central_id=get_single_value(central_id), plot_points=plot_points,
                              additional_points_size=additional_points_size,
                              additional_points_shape=additional_points_shape,
                              additional_points_color=additional_points_color, unit_length=units_length,
                              unit_velocity= units_velocity, unit_density= units_density,
                              plot_velocities=plot_velocities, plot_bfld= plot_bfld, axes=get_single_value(axes_array),
                              modified_units=modified_units, ignore_types=ignore_types,
                              factor_value=factor_value[0], units_value=units_value[0], contour=contour)

            regularize_time_units(loaded_snap)
            title('time : {:.2g}'.format(loaded_snap.time * basic_units["time"].factor) +
                  " [" + basic_units["time"].unit + "]")
            restore_basic_units(old_basic_units)
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
            old_basic_units = copy_current_units()
            for index,val in enumerate(value):
                if num_figures >= 1:
                    curr_subplot = int(num_figures*100 + 21 + index)
                print("curr subplot: ", curr_subplot)
                subplot(curr_subplot)
                curr_cmap = cmap[index % len(cmap)]
                plot_single_value(loaded_snap,  value=val, cmap=curr_cmap, box=get_single_value(box,index),
                                  vrange=get_single_value(vrange,index), logplot=get_single_value(logplot,index),
                                  res=res,
                                  numthreads=numthreads, center=get_single_value(center,index),
                                  relative_to_sink_id=get_single_value(relative_to_sink_id, index),
                                  central_id=get_single_value(central_id, index), plot_points=plot_points,
                                  additional_points_size=additional_points_size,
                                  additional_points_shape=additional_points_shape,
                                  additional_points_color=additional_points_color, unit_length=units_length,
                                  unit_velocity= units_velocity, unit_density= units_density,
                                  plot_velocities=plot_velocities, plot_bfld= plot_bfld, newfig=False,
                                  axes=get_single_value(axes_array, index), ignore_types=ignore_types,
                                  factor_value=factor_value[index % len(units_value)],
                                  units_value=units_value[index % len(units_value)], contour=contour)
                if index < len(value) - 1:
                    restore_basic_units(old_basic_units)
                rcParams.update({'font.size': 40, 'font.family': 'Serif'})
                rcParams['text.usetex'] = True

            #title('time : {:.2f} [s]'.format(loaded_snap.time))
            regularize_time_units(loaded_snap)
            suptitle('time : {:.2g}'.format(loaded_snap.time * basic_units["time"].factor) +
                     " [" + basic_units["time"].unit + "]", fontsize='x-large')
            restore_basic_units(old_basic_units)
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
    parser.add_argument('--center_x', nargs='+', type=float, help='point on x axis to be the center of the plot', default=None)
    parser.add_argument('--center_y', nargs='+', type=float, help='point on y axis to be the center of the plot', default=None)
    parser.add_argument('--center_z', nargs='+', type=float, help='point on z axis to be the center of the plot', default=0)
    parser.add_argument('--relative_to_sink_id', nargs='+', type=int,  help='id of sink particle to use as a reference point', default= None)
    parser.add_argument('--relative_to_id', nargs='+', type=int,  help='id of centeral particle', default= None)
    parser.add_argument('--plot_per_value_evolution', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot value evolution figure from the different snapshot of each value?',
                        default=False)
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
    parser.add_argument('--factor_value', nargs='+', type=float,  help='multiply value unit by this factor', default=[None])
    parser.add_argument('--units_value', nargs='+', type=str,  help='name of the value units', default=[None])
    parser.add_argument('--plot_contours', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='should plot contours?',
                        default=False)
    parser.add_argument('--snapshot_list', nargs='+', type=int,  help='list of snapshots to plot for '
                                                                      '(currently only for evolution)', default=[None])
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
        center = [[args.center_x, args.center_y,args.center_z]  for i in range(len(args.center_x))]

    axes_array = [[0,1]]
    if args.axes0 is not None and args.axes1 is not None:
        axes_array = [[args.axes0[i],args.axes1[i]] for i in range(len(args.axes0))]

    if args.snapshot_list[0] is None:
        snapshots_list = None
    else:
        snapshots_list = args.snapshot_list

    change_unit_conversion(args.factor_length, args.factor_velocity, args.factor_mass)
    #TODO: add conversion to temperature

    plot_range(args.value, args.source_dir, args.saving_dir, args.beginStep, args.lastStep, args.skipStep, box=box,
               vrange=vrange, logplot=args.logplot, cmap=args.cmap, res=args.res, numthreads= args.numthreads, center=center,
               relative_to_sink_id=args.relative_to_sink_id, central_id=args.relative_to_id,
               plot_points=args.plot_points, additional_points_size=args.additional_points_size,
               additional_points_shape=args.additional_points_shape,
               additional_points_color=args.additional_points_color,
               units_length=args.units_length, units_velocity=args.units_velocity, units_density= args.units_density,
               plot_velocities=args.plot_velocities, plot_bfld= args.plot_bfld, axes_array=axes_array,
               ignore_types=args.ignore_types, per_value_evolution=args.plot_per_value_evolution,
               factor_value=args.factor_value, units_value=args.units_value, contour=args.plot_contours,
               snapshots_list=snapshots_list)
