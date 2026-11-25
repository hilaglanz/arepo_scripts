# Proposed AREPO integration: a precomputed Mie grain-temperature table

This is a read-only implementation proposal. It has not been applied to AREPO,
compiled, or run. It targets the source currently in
`C:/Users/hila/GitRepository/arepo`.

**Scope after the numerical audit:** mode 1 below is a controlled sensitivity
experiment, not a recommended production default. Applying the 0.3-micron
unattenuated-temperature table and the current condensation logistic to the
saved 500-571 mean atmosphere gives a toy transport scattering depth about 400
and stellar-spectrum absorption depth about 1.9 (dust/gas=0.001). These are not
measured snapshot optical depths, but the proposed dust assignment fails its
own transparent-illumination check on that frozen atmosphere. A production
grain-temperature model must account for the radiation field and its spectral
weighting. Keep mode 0 for the existing physical prescription; mode 1 deliberately
tests the simpler radial closure.

**Opacity overlap identified in the follow-up audit:** the local Rosseland file
retains all 703 shared values of the original Ferguson GS98 X=0.7, Z=0.02 table
in `f05.gs98.tar.gz` (`g98.7.02.tron`), including all 323 shared values at
logT <= 3.5. The published Ferguson construction includes molecular opacity
and grain absorption plus scattering. This implementation adds an extra dust
term to that baseline; it does not remove the existing tabulated grains.
Consequently it is not a demonstrated non-overlapping gas-plus-dust model.
Planck-file generation details still require verification. Turning off an
opacity-blending switch cannot remove grains already in the table.

No additional compile-time flag is needed. Keep `BOWEN_OPACITY` and use an explicit
runtime mode: `0` retains the existing gas-temperature prescription, `1` reads a
grain-temperature table. The existing parameter reader treats registered tags as
required; the three new tags below must therefore appear even in mode 0. This
avoids silently switching old configurations to a new physical model.

The table supplies only the temperature used in the added Bowen scattering
factor. It does not change gas temperature, gas opacity lookups, the energy
equation, or the existing photospheric opacity switches.

## Table and parameter contract

The table has two whitespace-separated numeric columns:

```text
# r_over_Rstar  Td_K
# Record Teff, composition, grain size/distribution, optical-constant source,
# absorption-equilibrium solver, and the independently computed opacity here.
```

Numeric rows must start at `r/Rstar = 1`, have strictly increasing radius,
positive finite temperature, and non-increasing temperature. Blank lines and
comments starting with `#` or `%` are accepted, including trailing comments.
Use a grid extending through all density-eligible cells; `1` to `1024` is a
reasonable supplied grid for the present box. Out-of-range lookups stop with a
clear error rather than silently inventing an extrapolation. The table is for
one fixed incident stellar spectrum and one selected grain population.

```text
BowenGrainTemperatureMode   1
BowenStarRadius_cm         1.748342799e13
BowenGrainTemperatureTable /absolute/path/to/mie_Td_Mg2SiO4_a0p3um_T3061.dat
```

For one **combined sensitivity trial**, rather than separate temperature and
normalization trials, also use the matched 0.3-micron Mie normalization:

```text
BowenOpacityKappamax        4.5363373545
DustCondensationTemp       1100
BowenOpacitydT              50
BowenOpacityRhomin          1e-15
```

This assumes amorphous Fe-free Mg2SiO4, effective grain radius 0.3 micron,
solid density 3.27 g/cm3, saturated dust/gas mass ratio 0.001 and the spectrum
recorded in the table. It does not predict grain sizes or abundances. A combined
trial tests the response to the whole altered prescription; it cannot attribute
that response independently to temperature and opacity normalization. The former
14.821 value has not been proven universally incorrect by this assumed example.

Snapshot 291 is a common comparison start only if it is the saved state from
which the old Bowen run began. Preserve the same restart/time/core treatment.
Check the illumination and opacity-overlap assumptions before interpreting a
long run as a physical dust-wind prediction.

To retain the legacy gas-temperature prescription, use mode `0`, radius `0`, and
table value `unused`. The stable logistic expression below is mathematically
equivalent to the original expression, but is not promised to be bitwise identical.
The path is resolved relative to the executable's working directory when it is
not absolute. These are new required tags, not built-in defaults.

`BowenStarRadius_cm` is the photospheric radius. Do not substitute
`IrtStarCoreRadius`, which describes the much smaller radiation-injection/core
region. The table's stellar spectrum must match the intended luminosity and
photospheric radius; regenerate the table when that spectrum or grain population
changes. Comments document provenance but this minimal loader cannot verify the
physics encoded in them.

## 1. Add parameters and a prototype

In `src/allvars.h`, inside the existing `OPACITY_MESA` / `BOWEN_OPACITY` block near
line 3429, after `BowenOpacityRhomin`, add:

```c
    int BowenGrainTemperatureMode;   /* 0: gas T; 1: Mie equilibrium table */
    double BowenStarRadius_cm;      /* photospheric radius, physical cm */
    char BowenGrainTemperatureTable[MAXLEN_PATH];
```

In `src/parameters.c`, within the existing Bowen parameter registrations near
line 2590, after the `BowenOpacityRhomin` registration, add:

```c
      strcpy(tag[nt], "BowenGrainTemperatureMode");
      addr[nt] = &All.BowenGrainTemperatureMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "BowenStarRadius_cm");
      addr[nt] = &All.BowenStarRadius_cm;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BowenGrainTemperatureTable");
      addr[nt] = All.BowenGrainTemperatureTable;
      id[nt++] = PARAM_STRING;
```

In `check_parameters()` in the same file, near line 4599, add:

```c
#if defined(OPACITY_MESA) && defined(BOWEN_OPACITY)
  if(All.BowenGrainTemperatureMode != 0 && All.BowenGrainTemperatureMode != 1)
    terminate("BowenGrainTemperatureMode must be 0 or 1");
  if(!isfinite(All.DustCondensationTemp) || All.DustCondensationTemp <= 0 ||
     !isfinite(All.BowenOpacitydT) || All.BowenOpacitydT <= 0 ||
     !isfinite(All.BowenOpacityKappamax) || All.BowenOpacityKappamax < 0 ||
     !isfinite(All.BowenOpacityRhomin) || All.BowenOpacityRhomin < 0)
    terminate("Invalid Bowen temperature, width, opacity, or density threshold");
  if(All.BowenGrainTemperatureMode == 1)
    {
#ifndef IDORT
      terminate("Bowen grain-temperature table mode requires IDORT in this implementation");
#endif
      if(All.ComovingIntegrationOn)
        terminate("Bowen grain-temperature table mode requires non-comoving stellar coordinates");
      if(!isfinite(All.BowenStarRadius_cm) || All.BowenStarRadius_cm <= 0 ||
         !isfinite(All.UnitLength_in_cm) || All.UnitLength_in_cm <= 0)
        terminate("BowenStarRadius_cm and UnitLength_in_cm must be positive and finite");
    }
#endif
```

In `src/mesa_opacity.h`, beside `update_opacity()` near line 123, inside its
existing `OPACITY_MESA` guard, add:

```c
#ifdef BOWEN_OPACITY
void bowen_grain_temperature_init(void);
#endif
```

## 2. Add the table loader and interpolation in mesa_opacity.c

Add `#include <string.h>` with the standard includes. Insert the following near
the top of `src/mesa_opacity.c`, after the existing static table declarations and
inside its existing `OPACITY_MESA` guard. The arrays are process-local, immutable
after startup, and have no mutable interpolation cache. No allocation, MPI,
file I/O, or fatal-error call occurs in the threaded lookup.

```c
#ifdef BOWEN_OPACITY
#define BOWEN_TD_MAX_ROWS 8192
static int bowen_td_n = 0;
static double bowen_td_logx[BOWEN_TD_MAX_ROWS];
static double bowen_td_logt[BOWEN_TD_MAX_ROWS];

void bowen_grain_temperature_init(void)
{
  bowen_td_n = 0;
  if(All.BowenGrainTemperatureMode == 0)
    return;

  if(ThisTask == 0)
    {
      FILE *fp = fopen(All.BowenGrainTemperatureTable, "r");
      if(!fp)
        terminate("Cannot open Bowen grain-temperature table: %s",
                  All.BowenGrainTemperatureTable);

      char line[2048];
      int lineno = 0;
      double previous_x = 0.0, previous_t = 0.0;
      while(fgets(line, sizeof(line), fp))
        {
          lineno++;
          if(!strchr(line, '\n') && !feof(fp))
            terminate("Bowen temperature table line %d is too long", lineno);

          char *comment = strpbrk(line, "#%");
          if(comment)
            *comment = '\0';
          char *p = line;
          while(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
            p++;
          if(*p == '\0')
            continue;

          double x, t;
          char extra;
          if(sscanf(p, "%lf %lf %c", &x, &t, &extra) != 2 ||
             !isfinite(x) || !isfinite(t) || x < 1.0 || t <= 0.0)
            terminate("Invalid Bowen temperature table row at line %d", lineno);
          if(bowen_td_n == 0 && fabs(x - 1.0) > 1e-12)
            terminate("Bowen temperature table must start at r/Rstar=1");
          if(bowen_td_n > 0 && (x <= previous_x || t > previous_t))
            terminate("Bowen table radius must increase and temperature must not increase (line %d)", lineno);
          if(bowen_td_n == BOWEN_TD_MAX_ROWS)
            terminate("Bowen temperature table exceeds %d rows", BOWEN_TD_MAX_ROWS);

          /* Set the first coordinate exactly to the photosphere. */
          bowen_td_logx[bowen_td_n] = bowen_td_n == 0 ? 0.0 : log(x);
          bowen_td_logt[bowen_td_n] = log(t);
          bowen_td_n++;
          previous_x = x;
          previous_t = t;
        }
      if(ferror(fp))
        terminate("Error reading Bowen temperature table: %s", All.BowenGrainTemperatureTable);
      fclose(fp);
      if(bowen_td_n < 2)
        terminate("Bowen temperature table requires at least two numeric rows");
    }

  MPI_Bcast(&bowen_td_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(bowen_td_logx, bowen_td_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(bowen_td_logt, bowen_td_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  mpi_printf("BOWEN: loaded %d Td rows from %s; r/Rstar in [1,%g], Rstar=%g cm\n",
             bowen_td_n, All.BowenGrainTemperatureTable,
             exp(bowen_td_logx[bowen_td_n - 1]), All.BowenStarRadius_cm);
}

/* Returns NaN for invalid/out-of-range input; caller handles errors after OpenMP. */
static double bowen_grain_temperature_lookup(double radius_ratio)
{
  if(bowen_td_n < 2 || !isfinite(radius_ratio) || radius_ratio < 1.0)
    return NAN;
  const double lx = log(radius_ratio);
  if(lx < bowen_td_logx[0] || lx > bowen_td_logx[bowen_td_n - 1])
    return NAN;

  int lo = 0, hi = bowen_td_n - 1;
  while(hi - lo > 1)
    {
      const int mid = lo + (hi - lo) / 2;
      if(lx < bowen_td_logx[mid])
        hi = mid;
      else
        lo = mid;
    }
  const double f = (lx - bowen_td_logx[lo]) /
                  (bowen_td_logx[hi] - bowen_td_logx[lo]);
  return exp(bowen_td_logt[lo] + f * (bowen_td_logt[hi] - bowen_td_logt[lo]));
}
#endif
```

## 3. Initialize before the first opacity calculation

In `src/begrun.c`, add this inside the `OPACITY_MESA` initialization block near
lines 495-506, after the existing Planck initialization and before the closing
`#endif /* #ifdef OPACITY_MESA */`:

```c
#ifdef BOWEN_OPACITY
  bowen_grain_temperature_init();
#endif
```

This location is in `begrun1()`, not `begrun2()`. `begrun1()` runs before `init()`,
which can already call initial IDORT relaxation, and before checkpoint loading.
The new table is process-local static storage, not checkpointed pointers. On a
checkpoint restart the runtime fields must be recopied from the new parameter
file as in section 5, so the parameters continue to match the table loaded here.
Mode 0 performs no table I/O. Each program startup reloads the selected table.

## 4. Use the grain temperature only in the added scattering term

Keep `SphP[i].KappaScatter = 0.0;` and replace only the existing
`#ifdef BOWEN_OPACITY` block near `src/mesa_opacity.c:733` with:

```c
#ifdef BOWEN_OPACITY
    if(SphP[i].Density > All.BowenOpacityRhomin)
      {
        double grain_temperature = SphP[i].EOSTemperature;
        if(All.BowenGrainTemperatureMode == 1)
          {
#ifdef IDORT
            const double r_cm = get_distance_from_core(i) * All.UnitLength_in_cm;
            /* No added dust at or inside the stellar photosphere. */
            if(r_cm <= All.BowenStarRadius_cm)
              return;
            grain_temperature = bowen_grain_temperature_lookup(r_cm / All.BowenStarRadius_cm);
#else
            /* check_parameters() rejects this runtime mode in non-IDORT builds. */
            grain_temperature = NAN;
#endif
          }

        if(!isfinite(grain_temperature) || grain_temperature <= 0.0)
          {
            SphP[i].KappaScatter = NAN;
            return;
          }
        const double temperature_offset =
            (grain_temperature - All.DustCondensationTemp) / All.BowenOpacitydT;
        const double exp_factor = exp(-fabs(temperature_offset));
        const double condensed_fraction = temperature_offset >= 0.0
            ? exp_factor / (1.0 + exp_factor)
            : 1.0 / (1.0 + exp_factor);
        SphP[i].KappaScatter += All.BowenOpacityKappamax * condensed_fraction;
      }
#endif
```

`set_opacity()` currently ends immediately after this block, so the early
returns preserve all gas-opacity assignments already made. The original density
threshold retains its existing code-density units. In this run all base units
are 1, so it equals the cgs density threshold. No extra radius cut tied to
`IrtStarDustMinRadius` is introduced.

To report a lookup failure outside the parallel region, add
`int bowen_errors = 0;` in `update_opacity()` near line 747, change its OpenMP
pragma to:

```c
#pragma omp parallel for private(idx, i) reduction(+:bowen_errors)
```

Immediately after `set_opacity(SphP, i);` inside that loop, add:

```c
#ifdef BOWEN_OPACITY
        if(!isfinite(SphP[i].KappaScatter))
          bowen_errors++;
#endif
```

Immediately after the loop, add:

```c
#ifdef BOWEN_OPACITY
    if(bowen_errors)
      terminate("Bowen dust-temperature evaluation failed in %d active cells on this task; "
                "check temperatures/source geometry and extend the table radius range if needed",
                bowen_errors);
#endif
```

The interpolation writes no shared state. Errors are marked per cell, counted by
an OpenMP reduction, and reported only after threads have joined. No MPI
collective is added to the opacity-update loop.

## 5. Copy all Bowen runtime settings on checkpoint restart

Inside `reread_params_after_loading_restart()` in `src/restart.c`, within its
`OPACITY_MESA` block near lines 747-756, add:

```c
#ifdef BOWEN_OPACITY
  All.DustCondensationTemp = all.DustCondensationTemp;
  All.BowenOpacitydT = all.BowenOpacitydT;
  All.BowenOpacityKappamax = all.BowenOpacityKappamax;
  All.BowenOpacityRhomin = all.BowenOpacityRhomin;
  All.BowenGrainTemperatureMode = all.BowenGrainTemperatureMode;
  All.BowenStarRadius_cm = all.BowenStarRadius_cm;
  strcpy(All.BowenGrainTemperatureTable, all.BowenGrainTemperatureTable);
#endif
```

Here `all` is the saved new parameter-file configuration; `All` was overwritten
by the checkpoint. The existing four Bowen parameters also need these copies:
the current source does not permit them to change on checkpoint restart.

Adding fields changes the binary checkpoint layout. Do not assume an old
executable's raw restart files can be read by this modified executable. Start
the modified build from an IC or compatible HDF5 snapshot first. Later raw
checkpoints written by this same layout can use the runtime mode and table
changes above. Preserve/version the table file for reproducibility.

## Geometry, scope, and verification still required

- `get_distance_from_core()` in `src/giant_stars/rg.c:81` uses the cached donor
  position when `ID_RGCORE` is defined; otherwise it uses box center. Existing
  ID_RGCORE initialization and timestep updates populate that cache before RT.
  A moving or binary donor must have the correct source ID in its build.
- The table assumes the incident spectrum and dilution model used in its
  generator. It does not dynamically include shielding, backwarming, companions,
  or changes of stellar spectrum. It is not a replacement for multigroup dust
  radiative equilibrium in an arbitrary radiation field.
- Changing only this block separates the extra Bowen component's temperature
  from gas temperature. Grain opacity already inside the Ferguson tables and
  the existing gas-temperature opacity blends are still present.
- No new per-cell hydrodynamic variable or dust energy equation is introduced.
  The added scattering uses the existing IDORT coupling. Optional Td output
  should be labelled as an equilibrium estimate, not an evolved dust state.
- Before production use, compile the relevant configurations and verify table
  endpoints/interpolation, mode-0 equivalence, density/inside-star gates, MPI
  consistency, OpenMP consistency, and restart parameter persistence. None of
  those executable checks has been performed for this proposal.
