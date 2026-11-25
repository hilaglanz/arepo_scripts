# Actual opacity-table provenance and overlap audit

The accessible Rosseland table retains Ferguson GS98 opacity values at every tested original grid point. A separate `KappaScatter` field initially equal to zero does **not** establish that scattering is absent from the tabulated opacity. Adding an assumed large-grain population cannot currently be certified free of overlap with the baseline grains.

## Direct local evidence

The AGB parameter file names `gs98_z0.02_x0.7_ma_ross.data` and `gs98_z0.02_x0.7_ma_planck.data` at [param.txt:124](C:/Users/hila/Dropbox/RMHD/IDORT/AGB/1.5/param.txt:124).

Accessible table directory:

`C:\Users\hila\Dropbox\RMHD\IDORT\RSGg_m10lgl4.6_res6ib03n3s_r6s_new4_kap\RSGg_m10lgl4.6_res6ib03n3s_r6s_new4_kap`

Both table headers, line 1, state:

> fixed metal distribution -- based on Ferguson, Alexander, et al (2005) with GS98 photospheric metals

Line 3 specifies X=0.700000, Z=0.020000, 37 logR columns and 180 logT rows. Neither table includes component-resolved grain absorption/scattering columns.

I also found the local archive `C:\Users\hila\Dropbox\RMHD\IDORT\AGB\f05.gs98.tar.gz`. Streaming member `g98.7.02.tron` to memory yields the header:

> Grevesse & Sauval 1998 with X= 0.700000 and Z= 0.020000

Comparison with the used `*_ma_ross.data` at common temperature/density grid coordinates gives:

| Comparison | Shared cells | Exactly equal | Maximum absolute log-opacity difference |
|---|---:|---:|---:|
| All common coordinates, logT 2.7–4.5 | 703 | 703 | 0 |
| Common coordinates with logT <= 3.5 | 323 | 323 | 0 |

There are 37 shared temperature values and 19 original density coordinates. The archived member has 85 temperature rows; the used table has a different interpolated/extended grid. This numerical match confirms retention of the original tabulated values at the tested points, including the cool grain regime. It does not recover every processing step, authenticate the archive's download history, or mathematically exclude custom modifications between those points.

No generation script or component-resolved provenance for the custom `*_ma_planck.data` was located in the inspected repository/table directories. Its matching Ferguson header supports its stated origin but does not establish whether its original extinction opacity was changed to absorption-only. **No hashes of current cluster copies were verified.** Findings apply to the accessible local files.

Reproduce the read-only comparison and local file hashes with [audit_ferguson_table_provenance.ps1](audit_ferguson_table_provenance.ps1). It streams archive content with `tar -xOzf`; it never extracts files or loads snapshots. The recorded run is [ferguson_table_provenance_result.json](ferguson_table_provenance_result.json).

## What the primary documentation says

MESA explicitly describes the Ferguson low-temperature opacity family as including molecules and grains. Its opacity module returns a Rosseland radiative mean, with separate higher-temperature data and blends: [MESA kap overview](https://docs.mesastar.org/en/latest/kap/overview.html).

[Ferguson et al. (2005)](https://arxiv.org/pdf/astro-ph/0502045), sections 3.3–3.4, includes molecular lines and collision-induced absorption, and computes grain opacity using total extinction, Qext=Qabs+Qsca, for an MRN size distribution of spherical grains. Section 5 defines its Planck mean using the same monochromatic-opacity notation as its Rosseland mean. Thus the published family is not a gas-only or molecule-free opacity source; the local custom Planck file still needs its own processing history to establish its absorption/scattering interpretation.

## How the inspected AREPO source uses it

- [mesa_opacity.c:679](C:/Users/hila/GitRepository/arepo/src/mesa_opacity.c:679) evaluates the Rosseland table at gas density and gas temperature; lines 684–687 evaluate the Planck table. It does not decompose either into molecules, grains, absorption and scattering.
- [mesa_opacity.c:689](C:/Users/hila/GitRepository/arepo/src/mesa_opacity.c:689) blends Ross toward Planck beyond the configured radius, using gas-temperature thresholds 1200/1500 K. This changes the mean used; it does not remove table dust.
- [mesa_opacity.c:712](C:/Users/hila/GitRepository/arepo/src/mesa_opacity.c:712) implements the density-dependent photospheric switch; with `IRT_STAR_PH_OPACITY_SWITCH_TO_ROSS`, it assigns Planck toward the already modified Ross value at low density.
- [mesa_opacity.c:731](C:/Users/hila/GitRepository/arepo/src/mesa_opacity.c:731) initializes a *separate* scattering coefficient to zero, then adds the Bowen term. The legacy term uses gas temperature and a density gate. It contains no molecule-specific chemistry or cross-section calculation.
- [IDORT_run.c:733](C:/Users/hila/GitRepository/arepo/src/idort/IDORT_run.c:733) uses `KappaRoss + KappaScatter` as the total extinction-like coefficient, `KappaPlanck` in the thermal source, and their difference in the redistribution term. [IDORT_run.c:1410](C:/Users/hila/GitRepository/arepo/src/idort/IDORT_run.c:1410) uses the same combination in the implicit solver. The added coefficient therefore changes radiation transport as well as radiation momentum coupling; the source comment 'purely as momentum' is not a complete description of its numerical effects.
- [IDORT_run.c:1187](C:/Users/hila/GitRepository/arepo/src/idort/IDORT_run.c:1187) includes the Ross-plus-scattering combination in a radiation momentum-source branch; [IDORT_implicit_finite_volume_solver.c:654](C:/Users/hila/GitRepository/arepo/src/idort/IDORT_implicit_finite_volume_solver.c:654) passes that combination into the transport state.

The difference of independently weighted Rosseland and Planck means is not an independently measured physical scattering component. Nor does subtracting a generic grain constant from a Rosseland mean reconstruct a gas-only table: the harmonic mean is not additive in species opacities.

## Example actual coefficients

At rho=1e-13 g/cm³, beyond the dust-switch radius, in the low-density photospheric-switch regime, bilinear log-table interpolation gives the following cm²/g values. “Effective table” is both stored Ross and Planck after the enabled switches; it is still evaluated at the gas temperature.

| Gas T (K) | Raw Ross | Raw Planck | Effective table | Legacy Bowen scattering |
|---:|---:|---:|---:|---:|
| 1100 | 0.114085 | 1.740383 | 1.740383 | 7.410500 |
| 1500 | 0.00411988 | 1.165183 | 0.00411988 | 0.00497022 |
| 2000 | 0.0000369742 | 0.259969 | 0.0000369742 | 0.000000225724 |

These aggregate coefficients do not reveal how much of each table value is molecular versus grain opacity. A small Rosseland value does not prove that molecular opacity or dust scattering is absent.

The defensible statement is that an added large-grain term represents an explicit phenomenological addition or sensitivity model until a consistent baseline/population decomposition is supplied. It is not established as a correction for missing molecules, and a changed size distribution is not by itself proof that the same condensible material has not already contributed to the baseline tables.

Read-only audit: no source edits, installations, simulation runs or snapshot reads.
