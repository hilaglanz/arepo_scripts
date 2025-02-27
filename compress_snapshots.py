import gzip
import sys,os
from make_plots import *

def compress_one_by_one(snapshotDir, firstSnap, lastSnap, skipSteps, output_dir):
    snapshots_list = get_snapshot_number_list(snapshotDir, "snapshot_", firstSnap, lastSnap, skipSteps)
    for snap_i, snap in enumerate(snapshots_list):
        print("doing snapshot ", snap)
        with open(snapshotDir+"/snapshot_"+snap+".hdf5",'rb') as fin:
            with gzip.open(output_dir+"/snapshot_"+snap+".gz",'wb') as fout:
                fout.writelines(fin)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--beginStep', type=int,  help='first step', default=0)
    parser.add_argument('--lastStep', type=int,  help='last step', default=-1)
    parser.add_argument('--skipStep', type=int, help='number of steps to skip', default=1)
    parser.add_argument('--source_dir', type=str,  help='path to snapshot files directory', default= sys.argv[0])
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "zipped")


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    compress_one_by_one(args.source_dir, args.beginStep, args.lastStep, args.skipStep, args.saving_dir)