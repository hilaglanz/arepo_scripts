from load_modules import *
import os
import glob
import argparse
import numpy as np


def plot_profile_test(output_dir,snapshot_name,plotting_dir,testing_value="rho",log=True,snapshot_number_array=[0,8,10]):
    fig = figure(figsize=(36,20))
    rcParams.update({'font.size': 40, 'font.family': 'Serif'})
    rcParams['text.usetex'] = True
    rcParams['lines.linewidth'] = 3.0
    evenly_spaced_interval = np.linspace(0, 1, len(snapshot_number_array))
    line_colors = [cm.rainbow(x) for x in evenly_spaced_interval]
    labels=[]
    for snapshot_number in snapshot_number_array:
        s = gadget_readsnap(snapshot_number,output_dir,snapshot_name)
        s.plot_radprof(testing_value, log=log,color=line_colors[snapshot_number])
        labels.append("snap " + str(snapshot_number) + "," + str(s.time)+ " [s]")
    if len(snapshot_number_array) > 1:
        legend(labels)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--output_dir', type=str, help='diretory containing the snapshots', default="output")
    parser.add_argument('--snapshot_name', type=str, help='snapshots initilas', default="snapshot_")
    parser.add_argument('--plotting_dir', type=str, help='diretory to plot to', default="plots")



if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()