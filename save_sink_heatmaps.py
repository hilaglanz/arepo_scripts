
import pickle

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



def save_heatmap(axes, box, center, loaded_snap, numthreads, res, saving_file, value):
    if saving_file is not None:
        slice = loaded_snap.get_Aslice(value, center=center,
                                       box=box, res=res, numthreads=numthreads, axes=axes)
        posx = slice['x']
        posy = slice['y']
        slice_to_save = pylab.transpose(slice['grid'])
        with open(saving_file, 'wb') as pickle_file:
            pickle.dump(plotted_heatmap(posx, posy, slice_to_save), pickle_file)

def save_stream(posx, posy, velx, vely, saving_file):
    if saving_file is not None:
        with open(saving_file, 'wb') as pickle_file:
            pickle.dump(plotted_stream(posx, posy, velx, vely), pickle_file)
def save_scatter(axes, loaded_snap, point_pos, length_factor, saving_file):
    if saving_file is not None:
        scatter_saving_file = saving_file + "_scatter"
        with open(scatter_saving_file, 'wb') as pickle_file:
            pickle.dump(plotted_scatter(point_pos[axes[0]], point_pos[axes[1]],
                                        loaded_snap.parameters['SinkFormationRadius'] *
                                        length_factor), pickle_file)