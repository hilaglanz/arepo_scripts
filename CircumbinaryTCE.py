import argparse
import sys, os

from loadmodules import *

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--giant_snapshot_file', type=str, help='path to snapshot file containing the giant star', default="")
    parser.add_argument('--load_types', type=int, nargs='+', help='load only these types, '
                                                                  'if there is a point mass companion one '
                                                                  'should also load type 1', default=[0])
    parser.add_argument('--rhocut', type=float, help='lower cutoff of density from the snapshot', default=1)

    parser.add_argument('--relative_to_RL', type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
                        help='is the distance should be relative to RL size?',
                        default=True)
    parser.add_argument('--RL_factor', type=float, help='if relative to RL, by what factor?', default=2)

    parser.add_argument('--outer_separation', type=float, help='initial separation between the binary objects', default=None)
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="tce.ic.dat")
    return parser


if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
