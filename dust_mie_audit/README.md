# Bowen/Mie numerical audit

This folder contains offline calculations and an **unapplied AREPO implementation
proposal**. No AREPO source, IC, snapshot, or existing run configuration was changed.
No package was installed and no simulation was launched.

The calculation ran through the available Node JavaScript execution tool because
Python and Node executables were not on the shell PATH. It did **not** execute
miepython itself. `mie_audit.mjs` implements the homogeneous-sphere Mie series,
using downward logarithmic derivatives and Miller recurrence, and checks against
seven published miepython reference cases plus the small-particle Rayleigh limit.

## Model assumptions

- Optical data: Jena amorphous Mg2SiO4, Jager et al. (2003), `mg2.txt`.
- Each case uses one spherical grain radius; this is not a size distribution.
- Adopted solid density: 3.27 g/cm3 (an explicit normalization assumption).
- Stellar blackbody temperature: 3060.752 K; photospheric radius 1.748342799e13 cm.
- Scattering opacity uses the flux mean of Qsca*(1-g), per gram of gas.
- Reference dust/gas mass ratio: 0.001, an illustrative input, not measured here.
- Grain temperature uses Qabs-weighted heating/cooling and finite-disc dilution,
  with unattenuated stellar illumination. Gas collisions, diffuse/reprocessed
  radiation, grain growth and drift are absent.

## Calculated values

| Radius (micron) | kappa_tr for dust/gas=0.001 (cm2/g gas) | Td at 2 Rstar (K) | Radius for Td=1100 K (Rstar) | Dust/gas needed for kappa=14.821 |
|---:|---:|---:|---:|---:|
| 0.1 | 1.288064 | 792.04 | 1.45855 | 0.0115064 |
| 0.3 | 4.536337 | 865.69 | 1.58553 | 0.00326717 |
| 0.5 | 4.383205 | 903.64 | 1.64670 | 0.00338132 |
| 1.0 | 2.489396 | 925.71 | 1.67483 | 0.00595365 |

Opacity scales linearly with dust/gas and inversely with adopted solid density
at fixed radius/optical constants. The equilibrium temperature is independent
of dust abundance only under the prescribed transparent illumination model.
The 1100 K crossing is **not proof that grains actually form there**.

All seven reference cases passed their stated tolerances. Fourfold refinement
of the optical-data intervals changed these opacities by at most 3.6e-6 relative;
Td at 2-3 Rstar changed by less than 0.03 K. The full unweighted stellar Planck
integral is used as denominator, so missing spectral flux is not renormalized.
The refined quadrature covers that integral to approximately 3e-7.

`lookup_checks.json` records midpoint equilibrium residuals and exact-temperature
spot checks for the sample 513-row table. Checked interpolation errors near the
condensation transition are below 0.005 K. Long-wave truncation is checked as a
sensitivity test, not a rigorous bound on unmeasured optical constants.

## Critical physical limitation

Follow-up table-provenance audit: the local Rosseland file exactly preserves
all 703 original Ferguson GS98 X=0.7, Z=0.02 values on the shared grid,
including all 323 shared values at logT <= 3.5. The archived source is
`C:/Users/hila/Dropbox/RMHD/IDORT/AGB/f05.gs98.tar.gz`, member `g98.7.02.tron`.
The published source includes molecular opacity and grain absorption plus
scattering. Adding Bowen therefore cannot be justified as simply supplying
absent molecular scattering. Its independence from existing grain opacity
is unverified; the exact Planck-table construction remains unresolved.

The computed temperature is a transparent-environment estimate, not the actual
grain temperature in the AREPO radiation field. Applying the new temperature
logistic to the frozen 500-571 mean density profile yields approximately 88 g/cm2
of dust-weighted gas column: transport scattering depth about 400 and
stellar-spectrum absorption depth about 1.9 for the 0.3-micron example.
These are toy mean-profile values, not measured snapshot optical depths.
They nevertheless show that making this unattenuated-radius prescription the
production default is not justified by this calculation.

Large scattering depth changes mean intensity and angular transport; it does
not by itself imply exponentially extinguished bolometric flux. Reprocessed
radiation, backwarming, gas-grain coupling and finite growth can change the
temperature/formation result. A consistent production extension should determine
absorbed heating from the radiation field and invert a precomputed grain cooling
function. Grey bolometric radiation alone does not specify its required spectral
weighting. Existing grain contributions in the gas tables must also be resolved.

The old single-star donor remains an initial hydrodynamic model or legacy control.
It has not equilibrated to a different dust prescription. This audit does not
predict a mass-loss multiplier or require rebuilding the stellar interior.

## Files and reproduction

- `mg2.txt`: downloaded optical constants, positive first column in cm^-1;
  the zero-wavenumber row is excluded. Positive wavelengths span 0.1957-948.5 micron.
- `mie_audit.mjs`: numerical solver, reference checks, integrated models and
  sample table generator.
- `results_subdivisions_1.json`, `results_subdivisions_4.json`: numerical results.
- `mean_profile_optical_depth.mjs`, `.json`, `.md`: reproducible frozen-profile
  optical-depth calculation, numerical results/provenance, and readable report.
- `mie_Td_Mg2SiO4_a0p3um_T3061.dat`: sample Td lookup for 1-1024 Rstar, with
  two numeric columns (r/Rstar, Td/K) and comment metadata.
- `AREPO_mie_temperature_implementation.md`: exact manual integration proposal,
  using existing BOWEN_OPACITY and runtime mode 0=legacy / 1=radial table.
  It has not been compiled or run in AREPO.

With Node available, run this in the audit directory:

```sh
node --input-type=module -e "const m=await import('./mie_audit.mjs'); console.log(await m.runAudit('.',4)); console.log(await m.writeTemperatureTable('.'));"
```

Data SHA256: `C011A5C7C283C36F257143F41E75F29FC18A99D21CE54AA634CEEE0C068E597C`.

Sources:

- https://www2.astro.uni-jena.de/Laboratory/OCDB/data/silicate/amorph/mg2.txt
- https://www2.astro.uni-jena.de/Laboratory/OCDB/amsilicates.html
- https://doi.org/10.1051/0004-6361:20030916
- https://github.com/scottprahl/miepython/blob/main/miepython/mie_nojit.py
- https://github.com/scottprahl/miepython/blob/main/tests/test_mie.py
- https://miepython.readthedocs.io/en/latest/07_algorithm.html
- https://arxiv.org/pdf/1209.2112
