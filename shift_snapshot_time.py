import argparse
from loadmodules import *
from make_plots import get_snapshot_number_list

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--beginStep', type=int,  help='first step', default=0)
    parser.add_argument('--lastStep', type=int,  help='last step', default=-1)
    parser.add_argument('--skipStep', type=int, help='number of steps to skip', default=1)
    parser.add_argument('--source_dir', type=str,  help='path to snapshot files directory', default= sys.argv[0])
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "plots")
    parser.add_argument('--timeshift', type=float,  help='time in s to shift all snapshots', default= 0)

    return parser

def change_time(snapshotDir, savingDir, firstSnap, lastSnap, skipSteps, timeshift=0):
    if not os.path.exists(savingDir):
        os.mkdir(savingDir)
    for snap in get_snapshot_number_list(snapshotDir, "snapshot_", firstSnap, lastSnap, skipSteps):
        print("doing snapshot ", snap)
        loaded_snap = gadget_readsnap(snap, snapshotDir)
        loaded_snap.time += timeshift
        gadget_write_ics(savingDir+"/snapshot_"+snap, loaded_snap.data, double=True,format="hdf5",time=loaded_snap.time)

    print("Done shifting by ", timeshift, "[s].")

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    change_time(args.source_dir, args.saving_dit, args.beginStep, args.lastStep, args.timeshift)
