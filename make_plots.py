import os
import glob
import numpy as np
from loadmodules import *


def plot_range(value='rho',outputDir="output", firstSnap=0,lastSnap=-1,skipSteps=1,box=False,vrange=False,logplot=True, res=1024,
               center=True,plot_points=True, additional_points_size=30,additional_points_shape='X', additional_points_color='w', units_length = 'cm'):
    snapshots = glob.glob('./snapshot_*')
    maxSnap=len(snapshots)
    if lastSnap == -1:
        lastSnap = maxSnap
    else:
        lastSnap = min(lastSnap, maxSnap)
    if firstSnap > lastSnap:
        print("cannot do firstSnap > lastSnap")
        return

    for snap in range(firstSnap,lastSnap + 1,skipSteps):
        loaded_snap = gadget_readsnap(snap,outputDir)
        loaded_snap.plot_Aslice(value,logplot=logplot,colorbar=True, center= center, vrange=vrange, box=box, res=res)

        if box == False:
            box=[loaded_snap.boxsize,loaded_snap.boxsize]
        if plot_points:
            points = np.where(loaded_snap.type > 0)
            print("plotting points")
            for point in points:
                point_pos = loaded_snap.pos[point]
                np.scatter(point_pos[0], point_pos[1],additional_points_size, additional_points_color, additional_points_shape)

                if loaded_snap.type[point] == 5:
                    np.Circle((point_pos[0], point_pos[1]), loaded_snap.parameters['SinkFormationRadius']*res/box[0],
                              fill=False, color='white', linestyle='dashed', linewidth=3.0)


        np.x_label('x [' + units_length + ']' )
        np.y_label('y [' + units_length + ']' )
        np.title('time : '+ str(loaded_snap.parameters['TimeBetSnapshot'] * skipSteps * snap) + ' [s]' )
        np.savefig(outputDir + 'Aslice_' + value + '_' + str(snap) + '.jpg')
        print("saved fig")



def main(args):
    return

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    main(sys.argv)