from loadmodules import *
import argparse

def remove_tracers(snapshot_name, new_snapshot_name):
    snap = gadget_readsnapname(snapshot_name, loadonlytype=[0,1,3,4,5])
    gadget_write_ics(new_snapshot_name, snap.data)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--snapshot_file', type=str, help='path to snapshot file containing the tracersr', default="")
    parser.add_argument('--new_snapshot_file', type=str, help='path to new snapshot file without the tracersr', default="")

    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    remove_tracers(args.snapshot_file, args.new_snapshot_file)

