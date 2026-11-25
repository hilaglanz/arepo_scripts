# Diagnose the initial-shock reservoir and Bowen response

`analyze_bowen_launch.py` reads existing single-star, noncosmological CGS HDF5
snapshots with the `loadmodules` installed in the interpreter you use. It does not
modify ICs, snapshots, AREPO source, or run a simulation.

Snapshot 291 is the initial reservoir baseline, **not an assumed relaxed state**.
Include it along with 413, 479 and 496 in each branch. The recently supplied
no-Bowen CSV includes all four numbers (its cadence becomes every snapshot from
477 through 496). Use the actual snapshot directories, not the directory of CSVs.

```bash
arepoPy=/usr/local/ph_hagai/anaconda3/envs/arepo-env/bin/python3.10

"$arepoPy" analyze_bowen_launch.py --source-dir /path/to/bowen/output \
  --label bowen --snapshots 291 413 479 496 --output-dir launch_check --plot

"$arepoPy" analyze_bowen_launch.py --source-dir /path/to/no_bowen/output \
  --label no_bowen --snapshots 291 413 479 496 --output-dir launch_check --plot
```

The branches can be processed on different computers; bring their CSV and JSON
outputs together afterward. Existing output files with the same label are
overwritten, so use distinct labels. Missing requested snapshots are reported as
reader errors rather than silently substituted. Missing physical fields are
reported explicitly. Missing added scattering is treated as zero only when
the snapshot Config metadata establishes that Bowen was disabled.

## Outputs and resource use

For each label, the script saves:

- `*_profiles.csv`: 48 logarithmic radial bins from 300 to 3000 solar radii.
- `*_shells.csv`: 300-500, 500-700, 700-1000, 1000-1500, 1500-2000 and 2000-3000
  solar radii, useful for comparing reservoir inventories and fluxes directly.
- `*_metadata.json`: physical times, core reference frame, Bowen parameters,
  parameter provenance, relevant flags and the actual Python/reader paths.
- `*_profiles.png` when `--plot` is used; plotting requires matplotlib already
  present in the AREPO Python environment.

This reads one snapshot at a time and lazily loads one gas field at a time. A
monolithic field still occupies memory in full while it is being reduced: for
28 million cells a float64 vector is about 672 MB. Additional radii, indices and
selected-region arrays remain resident. `--block-size` limits arithmetic
temporaries, not the reader's complete field allocation. No thermodynamic fields
unrelated to this diagnostic are requested, and no gas-density filter is applied.

## Read these columns first

1. `gas_mass_Msun` and `delta_gas_mass_from_291_Msun`: how much gas occupies each
   fixed-width shell following the core, compared with the initial shocked state.
   A decrease alone cannot distinguish inward drainage from outward clearing;
   read it together with signed fluxes and the neighbouring shells.
2. `outward_Msun_per_yr`, `inward_Msun_per_yr`, `net_Msun_per_yr`: positive
   magnitudes for each direction, with net = outward - inward. These are radial
   averages of flux across finite shell widths, not escape measurements.
3. `density_eligible_mass_fraction` and `cool_and_dense_mass_fraction`: the
   fractions above the strict density gate and also at/below the gas-temperature
   condensation threshold. Hotter gas still has a nonzero logistic tail.
4. `stored_added_kappa_mass_mean_cm2_g` and
   `halfmax_added_opacity_mass_fraction`: the actual saved extra opacity and the
   mass fraction where it reaches half of the chosen maximum. This is an opacity
   threshold, not a measured physical condensed-dust fraction.
5. `gamma_base_mass_mean`, `gamma_added_mass_mean`, `gamma_total_mass_mean`:
   separate the stored base and added contributions. The base may already
   contain grain opacity and the code's Rosseland/Planck switching.
6. `fraction_of_inward_mass_gamma_total_gt1`: whether the infalling reservoir
   currently experiences a radiation-force proxy larger than spherical gravity.
   Gas can remain inward-moving while it decelerates. Values below one do not
   exclude partial radiative support.

Means are gas-mass weighted, except `density_volume_mean_g_cm3`, which is total
mass divided by the sum of cell volumes m/rho. Empty-bin means and fractions
with no infalling mass are NaN; empty-bin mass and flux are zero. Plot signed
quantities on linear or symmetric-log axes.

## What the trial calculation means

`trial_kappa_mass_mean_cm2_g`, `gamma_trial_total_mass_mean`, and
`fraction_of_inward_mass_gamma_trial_total_gt1` apply the **current gas-temperature
Bowen law** to the snapshot without changing its radiation field or gas state:

    rho > rhomin: kappa_trial = kappamax / (1 + exp((Tgas - Tcond)/dT))
    otherwise:  kappa_trial = 0
    Gamma_trial_total = (kappa_base + kappa_trial) * F_radial / (c*g)

This replaces the stored added opacity in the diagnostic; it is never added
twice. On the no-Bowen branch it asks whether its existing gas would be eligible
for a significant force under this law. It does not predict the dynamically
coupled result of enabling Bowen.

The script reads recorded parameters when available. For files lacking them,
the documented CLI fallbacks are `--temperature 1100 --width 50
--kappa-max 14.821 --rho-min 1e-15`. Match these fallbacks to the Bowen run and
check the printed/JSON parameter sources before comparing branches.
`IrtStarDustMinRadius` is not a gate in the current Bowen implementation; it
controls a separate opacity-table blend, so this script does not impose it.

Gravity includes the core plus all interior gas, including gas inside 300 solar
radii. It is a spherical approximation. Vector radiation flux is projected
about the actual core; the stored scalar `fradr` is not used. Velocities are
relative to the moving core, or a stationary shell for DMFIXED snapshots.
Gamma uses the laboratory-frame flux and omits pressure and other acceleration
terms; it is not the exact conservative momentum update used by IDORT.

The initial shock can supply a reservoir that Bowen partly supports or
re-accelerates. A response of that reservoir is physically distinct from a
later continuously replenished wind. These Eulerian measurements do not identify
the same material parcels across snapshots, especially with refinement and mesh
mass exchange. No particle-ID tracking or escape classification is inferred.

## Local validation

The numerical routines can be checked with synthetic arrays without importing
loadmodules or opening simulation snapshots:

```bash
python -m unittest discover -s tests -p 'test_analyze_bowen_launch.py' -v
```

The production reader and optional plotting still need to run in the actual
AREPO Python environment. The CSV metadata records which reader was used.
