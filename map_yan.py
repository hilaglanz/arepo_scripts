import argparse
import sys
from loadmodules import *

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--explosion_time', type=float,  help='time between ignition and last snapshot in seconds')
    parser.add_argument('--output', type=str,  help='path to output directory', default="../output")
    parser.add_argument('--snapshot_base', type=str,  help='snapshot base', default="snapshot_")
    parser.add_argument('--output_grid', type=str,  help='path to output directory for the new model', default="./grid/")

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    y = yag.yag(arepodir=args.output, arepobase=args.snapshot_base, outdir=args.output_grid)
    y.map3D(explosiontime=args.explosion_time, centering=True, map1DSNU=True, map2DSNU=False, maxVel=3e9,
            normaliseAbundances=True, map1D=False, removeHeShell=False,boxsize=y.get_boxsize())
    print("new model save in ", args.output_grid, "/model_SNU_1D.txt")
