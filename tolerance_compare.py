import argparse
from loadmodules import *

def compare_tolerance(snapshot_num, value="ka_r", output1="output-6", output2="output-8", plots_dir="plots", xrange=700,
                      val_label=None, val_units=""):
    run_1e6 = gadget_readsnap(snapshot_num, output1, lazy_load=True,loadonlytype=[0])
    run_1e8 = gadget_readsnap(snapshot_num, output2, lazy_load=True,loadonlytype=[0])
    if val_label is None:
        val_label = value

    r_vec6 = run_1e6.pos - run_1e6.center
    r_mag6 = np.linalg.norm(r_vec6, axis=1)
    v_rad6 = np.sum(run_1e6.vel * r_vec6, axis=1) / r_mag6
    run_1e6.data["vr"] = v_rad6

    r_vec8 = run_1e8.pos - run_1e8.center
    r_mag8 = np.linalg.norm(r_vec8, axis=1)
    v_rad8 = np.sum(run_1e8.vel * r_vec8, axis=1) / r_mag8
    run_1e8.data["vr"] = v_rad8

    # 1. Extract the radial profiles
    # We use 'vel' for the radial velocity check
    dr_in_rsum = 5
    dr_val = dr_in_rsum * rsol
    dynamic_nshells = int(xrange / dr_in_rsun) + 10
    prof8= run_1e8.get_radprof("vr", nshells=dynamic_nshells, dr=dr_val)
    r_1e8 = prof8[1,:]
    v_1e8  = prof8[0,:]
    prof6= run_1e6.get_radprof("vr", nshells=dynamic_nshells, dr=dr_val)
    r_1e6 = prof6[1,:]
    v_1e6  = prof6[0,:]

    r_1e8 /= rsol
    r_1e6 /= rsol
    v_1e8 /= 10**5
    v_1e6 /= 10**5

    # 2. Extract Opacity (to see where the dust starts)
    # This helps correlate error with the dust formation region
    prof_value = run_1e8.get_radprof(value, nshells=dynamic_nshells, dr=dr_val)
    value_1e8 = prof_value[0,:]

    # np.interp MUST have increasing x-values
    sort_idx = np.argsort(r_1e6)
    r_1e6_sorted = r_1e6[sort_idx]
    v_1e6_sorted = v_1e6[sort_idx]

    # Remove NaNs if they exist (common at the center/sink boundary)
    mask = ~np.isnan(v_1e6_sorted)
    r_clean = r_1e6_sorted[mask]
    v_clean = v_1e6_sorted[mask]

    # 3. Handle interpolation
    # (In case the radial bins aren't exactly the same between runs)
    v_1e6_interp = np.interp(r_1e8, r_clean, v_clean)

    # 4. Calculate Relative Residual
    # We add a tiny epsilon to avoid division by zero in static regions
    epsilon = 1e-5
    residual = np.abs(v_1e8 - v_1e6_interp) / (np.abs(v_1e8) + epsilon)
    # --- Plotting ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 11), sharex=True,
                                   gridspec_kw={'height_ratios': [2, 1]})

    # Top: Physics Comparison
    ax1.plot(r_1e8, v_1e8, 'k-', label='Tolerance 1e-8', alpha=0.8)
    ax1.plot(r_1e8, v_1e6_interp, 'r--', label='Tolerance 1e-6', alpha=0.8)
    ax1.set_ylabel('Radial Velocity [km/s]')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)

    # Overlay Opacity (Secondary Axis) to show the Dust Zone
    ax1_tw = ax1.twinx()
    ax1_tw.plot(r_1e8, value_1e8, 'g:', alpha=0.5, label=val_label)
    ax1_tw.set_yscale('log')
    ax1_tw.set_ylabel(val_label + " " + val_units, color='g')

    # Bottom: Numerical Error (Residual)
    ax2.plot(r_1e8, residual * 100, color='blue', lw=1.5)
    ax2.axhline(y=1.0, color='darkgreen', linestyle='-', alpha=0.6, label='1% Error')
    ax2.axhline(y=0.1, color='orange', linestyle=':', alpha=0.6, label='0.1% Error')

    ax2.set_yscale('log')
    ax2.set_ylim(1e-3, 1000) # Show from 0.001% to 100%
    ax2.set_ylabel('Relative Difference [%]')
    ax2.set_xlabel('Radius [Rsun]')
    ax2.legend(loc='lower right', fontsize='small')
    ax2.grid(True, which="both", ls="-", alpha=0.2)

    plt.suptitle('IDORT Tolerance Convergence Test: 1e-8 vs 1e-6', fontsize=14)
    plt.tight_layout()
    plt.xlim(0, xrange)
    savefig(plots_dir + "/Tolerance_comp_{0}_".format(snapshot_num)+"_"+value+".jpg")
    print("saved!")

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--snapshot_num', type=int,  help='snapshot number', default=0)
    parser.add_argument('--output1', type=str,  help='path to snapshot files directory of the low tolerance', default= "output-6")
    parser.add_argument('--output2', type=str,  help='path to snapshot files directory of the high tolerance', default= "output-8")
    parser.add_argument('--saving_dir', type=str,  help='path to output directory', default= "plots")
    parser.add_argument('--value', type=str,  help='value to be plotted', default="ka_r")
    parser.add_argument('--val_label', type=str,  help='value label to be plotted', default="Opacity (Dust)")
    parser.add_argument('--val_units', type=str,  help='value units to be plotted', default="[cm^2/g]")
    parser.add_argument('--xrange', type=float, help='x range to plot (in Rsun)', default=1000)
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()

    compare_tolerance(snapshot_num=args.snapshot_num, output1=args.output1, output2=args.output2, xrange=args.xrange,
                      value=args.value, val_label=args.val_label, val_units=args.val_units, plots_dir=args.saving_dir)

