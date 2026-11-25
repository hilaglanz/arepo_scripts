# Bowen opacity: gas or grey radiation temperature

`OPACITY_MESA` and `BOWEN_OPACITY` enable the additional effective scattering
opacity. The optional integer parameter `BowenTemperatureMode` selects its
activation temperature:

- `0` (default): the existing gas-temperature law. Omitting the parameter keeps
  the old behavior and records `0` in `parameters-usedvalues` and snapshot
  parameter metadata. The numerical expression is unchanged.
- `1`: use `T = (E_rad / RAD_CONST)^(1/4)` for the additional Bowen term only.
  Gas EOS and table opacities still use gas temperature. This mode requires
  IDORT, non-comoving CGS units (mass, length and velocity units all equal to 1),
  and a production build without `IDORT_TEST`.

Both modes retain the strict `Density > BowenOpacityRhomin` gate. The threshold
is in code density units; these equal g/cm^3 in the supported mode-1 setup.
Opacity is in cm^2 per gram of gas. Condensation temperature and width must be
positive and finite; maximum opacity and the density threshold must be
nonnegative and finite. The new radiation mode uses an overflow-safe logistic.

## Combined sensitivity example

Add/replace these entries in the intended single-star parameter file; do not
append duplicate entries:

```text
BowenTemperatureMode       1
BowenOpacityKappamax       4.5363373545
DustCondensationTemp       1100
BowenOpacitydT             50
BowenOpacityRhomin         1e-15
```

This illustrative Mie normalization assumes spherical, amorphous Fe-free
Mg2SiO4 grains of radius 0.3 micron, solid density 3.27 g/cm^3, dust/gas mass
ratio 0.001, and a 3060.752 K incident blackbody. It is the stellar-spectrum
mean of transport scattering, proportional to `Q_sca * (1 - g)`, not raw
scattering or total extinction. It is not a universal replacement for 14.821.
The 0.5-micron example instead gives 4.383205 cm^2/g. Grain size and abundance
are assumptions, not predicted dust-formation variables. No parameter file is
automatically modified and no Mie calculation occurs per cell or timestep.

## Initialization and restarts

For an HDF5 snapshot continuation (restart mode 2), build with `OUTPUT_IDORT`
and without `OUTPUT_IDORT_NO_IV`. The snapshot must contain gas `E_rad` and `IV`.
Missing datasets abort instead of being silently padded with zeros. The saved
radiation field sets the new opacity immediately; no new initial RT relaxation
is requested for a snapshot restart.

For a fresh IC (mode 0), zero/uninitialized radiation may temporarily use gas
temperature for the first opacity estimate. Positive radiation, if actually
read, is used immediately. In the EOS_OPAL input path radiation fields are
normally skipped for fresh ICs unless `READ_IN_ALL_IC_FIELDS` is enabled.
After the initial radiation solve, opacity and optical depth are refreshed.
Consequently mode 1 with a fresh IC requires `IDORT_NO_INI` to be disabled.

Outside this explicit bootstrap, density-eligible cells require positive finite
`E_rad`. Negative, nonfinite, or zero radiation is reported after leaving the
OpenMP loop. Strict positivity is a deliberate illuminated-stellar-atmosphere
contract; it is not a general statement that physically dark cells cannot exist.
This adds no per-cell stored fields or MPI collectives to the opacity loop.

All five Bowen settings can be changed through the parameter file when resuming
a raw checkpoint written by this build. A change refreshes the added scattering
and optical depth of all local gas cells, including inactive cells, before the
next radiation solve. An unchanged prescription preserves the stored state.
For `ONEDIMS_SPHERICAL`, change the prescription through snapshot mode 2 instead:
its cell radii require mesh faces that are not yet available during raw restart.
Adding the new global parameter changes
the binary checkpoint layout: start this build from an HDF5 snapshot, not an
older executable's raw restart file. Keep the original simulation time and
core treatment when comparing continuations from snapshot 291.

## Numerical and physical scope

The existing timestep order is retained: solve radiation, then update opacity
from the resulting moments for the next solve. This is a lagged-opacity
approximation. The inner IDORT iteration does not reconverge dust opacity;
`IDORT_IT_COEVOL` does not change that. Check a short continuation with a smaller
timestep before an expensive run. A fully coupled dust-opacity iteration would
need to rebuild the RT coefficients, not simply call `set_opacity` again.

`E_rad` is already radiation energy density: do not divide it by c or replace
it with `abs(F)/c`. The local energy incorporates the angular transport and
trapping represented by the grey solver. It does not contain the spectrum needed
for material-specific, absorption-weighted Mie grain equilibrium. This mode is
therefore a grey Bowen sensitivity model, not a full silicate dust model.

The additional `KappaScatter` affects both radiation transport and momentum.
It adds no true dust absorption/emission term; adding such a term to the current
gas-temperature thermal source alone would not describe separately heated
grains. Grain opacity already present in the baseline tables is also not
removed. These changes do not resolve that potential overlap.

## Focused validation

From the repository root, with Python 3 and a C99 compiler available:

```sh
python3 tests/check_bowen_opacity.py
```

The small harness exercises the actual opacity functions extracted from the C
source using synthetic cells. It needs no MPI, opacity tables or snapshots.
It is not a substitute for compiling the intended full AREPO configuration and
checking a short MPI continuation. No initial conditions or simulations are
generated by this test.
