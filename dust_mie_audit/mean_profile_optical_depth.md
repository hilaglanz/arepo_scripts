# Frozen mean-profile optical-depth sensitivity audit

This is an offline toy calculation on the saved snapshots-500–571 averaged atmosphere, not an actual snapshot optical depth or a relaxed new donor. No snapshots or AREPO source are modified.

For dust/gas=0.001, integrate rho_bar times the logistic grain activation and the strict rho_bar>1e-15 gate. The logistic width is 50 K around 1100 K. Partial activation **inside** the nominal 1.5855 Rstar condensation radius is retained; there is no hard cutoff there. The upper integration radius is 49.724 Rstar.

| Outward from r/Rstar | Weighted gas column (g/cm²) | Scattering transport tau | Stellar-spectrum absorption tau | Raw scattering tau |
|---:|---:|---:|---:|---:|
| 1 | 88.462 | 401.29 | 1.9044 | 988.72 |
| 1.5855 | 47.224 | 214.22 | 1.0167 | 527.81 |
| 2 | 3.6591 | 16.599 | 0.078774 | 40.897 |
| 3 | 0.39047 | 1.7713 | 0.0084062 | 4.3642 |

The old gas-temperature logistic on the same mean profile has weighted column 0.0268456 g/cm² and tau=0.397878 with kappa_max=14.821 cm²/g. This also is opacity evaluated on averages, not mean cellwise opacity.

Linear-density interpolation gives total weighted column 88.46231 g/cm²; log-density interpolation gives 85.73229 g/cm² (-3.09%).

| Simpson subintervals per union-knot interval | Total weighted column (g/cm²) | Total transport tau | Legacy scattering tau |
|---:|---:|---:|---:|
| 1 | 88.46230835 | 401.29487384 | 0.397877466 |
| 4 | 88.46231054 | 401.29488375 | 0.397877901 |
| 16 | 88.46231054 | 401.29488379 | 0.397877903 |

The input profile spacing is about 0.0711 Rstar. Integration splits at input profile and lookup knots and density-gate crossings before Simpson refinement; convergence here establishes quadrature accuracy only, not physical resolution.

- Opacity evaluated on averaged gas temperature/density is not the average cellwise opacity; shocks, clumps and correlations are absent.
- The atmosphere is frozen while the new grain-temperature assignment is imposed; dynamics and grain growth could change the density and dust distribution.
- Stellar-spectrum flux-mean coefficients are held fixed; spectral transfer, diffuse radiation and backwarming are not solved.
- Large scattering optical depth changes angular intensity; conserved spherical bolometric flux is not exponentially attenuated by scattering alone.
- The thin-temperature lookup fails this frozen-profile thinness check. Use it as a diagnostic closure, not an already validated production dust-temperature model.

## Reproduction and provenance

Run `node mean_profile_optical_depth.mjs [profile_directory] [output_directory]`. Defaults use the original Dropbox profile directory and this script's directory. The adjacent JSON records constants, input paths and SHA-256 hashes, detailed inward/outward columns, and convergence.

- density: `C:/Users/hila/Dropbox/RMHD/IDORT/AGB/1.5/1e-8/dust/rho_profile__500_to_571_time_avg.txt`; SHA-256 `1533f5dace6e57beaea62565dfdf62d2193797897f03c2e4e678034dc6bcf84f`.
- gasTemperature: `C:/Users/hila/Dropbox/RMHD/IDORT/AGB/1.5/1e-8/dust/temp_profile__500_to_571_time_avg.txt`; SHA-256 `7e5ee0849d5dd68d47ddaba2010e8b561a296883abb3ac18d31ada993fd7d483`.
- grainTemperature: `C:/Users/hila/GitRepository/arepo_scripts/dust_mie_audit/mie_Td_Mg2SiO4_a0p3um_T3061.dat`; SHA-256 `f6df3ef3a40c29ff279051da82f3d0e1ce098a9e89d56f2c9d34f86f9a31a3cd`.
- mieCoefficients: `C:/Users/hila/GitRepository/arepo_scripts/dust_mie_audit/results_subdivisions_4.json`; SHA-256 `3c568edbd790f8bfb07a2f6222cd4597768991ea7e13922b333415c8ad9ba17d`.
