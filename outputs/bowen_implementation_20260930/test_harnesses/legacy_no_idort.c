#define OPACITY_PLANCK
#define BOWEN_OPACITY

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define RAD_CONST 7.5657e-15
struct sph_particle_data {
  double Density, EOSTemperature, E_rad;
  double KappaRoss, KappaPlanck, KappaScatter, tau;
};
struct {
  int BowenTemperatureMode;
  double BowenOpacityRhomin, BowenOpacityKappamax;
  double DustCondensationTemp, BowenOpacitydT, UnitDensity_in_cgs;
} All;
static struct sph_particle_data cells[3], *SphP = cells;
static int NumGas = 3;
static int active[4] = {-1, 0, 1, 2};
static struct { int NActiveParticles; int *ActiveParticleList; }
  TimeBinsHydro = {4, active};
struct opal_opacity_result { double kappa_ross; };
struct planck_opacity_result { double kappa_planck; };
static int opal_calls, terminate_calls, invalid_count, terminated_in_parallel;
static int calls_at_termination;
static void get_opal_opacity(double rho, double temperature,
                             struct opal_opacity_result *result) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
  opal_calls++;
  result->kappa_ross = rho + temperature * 0.001;
}
static void get_planck_opacity(double rho, double temperature,
                               struct planck_opacity_result *result) {
  result->kappa_planck = rho + temperature * 0.002;
}
static double get_cell_radius(int i) { return 2.0 + i; }
static void terminate(const char *format, ...) {
  va_list arguments;
  (void)format;
  va_start(arguments, format);
  invalid_count = va_arg(arguments, int);
  va_end(arguments);
  terminate_calls++;
  calls_at_termination = opal_calls;
#ifdef _OPENMP
  terminated_in_parallel = omp_in_parallel();
#endif
}
static int checks;
#define CHECK(condition) do { checks++; if(!(condition)) { \
  fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #condition); return 1; \
} } while(0)
static int close_enough(double a, double b) {
  return fabs(a-b) <= 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));
}
static void reset(void) {
  int i;
  All.BowenTemperatureMode = 0;
  All.BowenOpacityRhomin = 1e-15;
  All.BowenOpacityKappamax = 4.5363373545;
  All.DustCondensationTemp = 1100.0;
  All.BowenOpacitydT = 50.0;
  All.UnitDensity_in_cgs = 1.0;
  opal_calls = terminate_calls = invalid_count = terminated_in_parallel = 0;
  calls_at_termination = 0;
  TimeBinsHydro.NActiveParticles = 4;
  for(i=0; i<3; i++) {
    cells[i].Density = 1e-12;
    cells[i].EOSTemperature = 2000.0;
    cells[i].E_rad = RAD_CONST * pow(1100.0, 4.0);
    cells[i].KappaScatter = 999.0;
    cells[i].KappaRoss = 100.0 + i;
    cells[i].KappaPlanck = 200.0 + i;
    cells[i].tau = -99.0;
  }
}
#ifdef BOWEN_OPACITY
static double get_bowen_opacity(const struct sph_particle_data *cell, int allow_radiation_bootstrap)
{
  if(!(cell->Density > All.BowenOpacityRhomin))
    return 0.0;

  if(All.BowenTemperatureMode == 0)
    /* Preserve the legacy expression, including its rounding. */
    return All.BowenOpacityKappamax /
        (1.0 + exp((cell->EOSTemperature - All.DustCondensationTemp) / All.BowenOpacitydT));

  double temperature_for_dust = cell->EOSTemperature;
#ifdef IDORT
  const double radiation_energy = cell->E_rad;
  if(isfinite(radiation_energy) && radiation_energy > 0.0)
    {
      /* E_rad is already an energy density in physical CGS units.
       * This is a grey proxy, not a spectral Mie grain temperature. */
      temperature_for_dust = sqrt(sqrt(radiation_energy / RAD_CONST));
    }
  else if(!(allow_radiation_bootstrap && radiation_energy == 0.0))
    return NAN;
#else
  /* check_parameters() rejects radiation mode without IDORT. */
  return NAN;
#endif
  if(!isfinite(temperature_for_dust) || temperature_for_dust <= 0.0)
    return NAN;

  const double temperature_offset =
      (temperature_for_dust - All.DustCondensationTemp) / All.BowenOpacitydT;
  const double exponential = exp(-fabs(temperature_offset));
  const double condensation_fraction = temperature_offset >= 0.0
      ? exponential / (1.0 + exponential) : 1.0 / (1.0 + exponential);
  return All.BowenOpacityKappamax * condensation_fraction;
}
#endif
void set_opacity(struct sph_particle_data *SphP, int i, int allow_radiation_bootstrap)
{
//    double Temperature = SphP[i].EOSTemperature;
//#ifdef IRT_OPACITY_T_RAD
//    double rad_const;
//#ifdef IDORT_TEST
//    rad_const = All.IdortRadConst;
//#else
//    rad_const = RAD_CONST;
//#endif
//    Temperature = sqrt(sqrt(SphP[i].Erad / rad_const));
//#endif
    
    struct opal_opacity_result res;
    get_opal_opacity(SphP[i].Density, SphP[i].EOSTemperature, &res);
    SphP[i].KappaRoss = res.kappa_ross;
    SphP[i].KappaPlanck = res.kappa_ross;
    
#ifdef OPACITY_PLANCK
    struct planck_opacity_result res2;
    get_planck_opacity(SphP[i].Density, SphP[i].EOSTemperature, &res2);
    SphP[i].KappaPlanck = res2.kappa_planck;
    double factor;
#ifdef IRT_STAR_DUST_OPACITY_SWITCH
    double r = get_distance_from_core(i);
    double Tl = 1200.0;
    double Th = 1500.0;
    double Temp = SphP[i].EOSTemperature;
    factor = 1.0;
    if(r > All.IrtStarDustMinRadius)
    {
        if(Temp > Th)
        {
            factor = 1.0;
        }
        else if(Temp < Tl)
        {
            factor = 0.0;
        }
        else
        {
            factor = (Temp - Tl) / (Th - Tl);
        }
    }
    SphP[i].KappaRoss = factor * res.kappa_ross + (1.0-factor) * res2.kappa_planck;
#endif
#ifdef IRT_STAR_PH_OPACITY_SWITCH
      double rho = SphP[i].Density;
      factor = 0.0;
      if(rho > 10.0 * All.IrtStarPhKapDensity)
      {
          factor = 1.0;
      }
      else if(rho > 0.1 * All.IrtStarPhKapDensity && rho <= 10.0 * All.IrtStarPhKapDensity)
      {
          factor = (rho - 0.1 * All.IrtStarPhKapDensity) / (10.0 * All.IrtStarPhKapDensity - 0.1 * All.IrtStarPhKapDensity);
      }
#ifdef IRT_STAR_PH_OPACITY_SWITCH_TO_ROSS
    SphP[i].KappaPlanck = factor * res2.kappa_planck + (1.0-factor) * SphP[i].KappaRoss;
#else
    SphP[i].KappaRoss = factor * res.kappa_ross + (1.0-factor) * res2.kappa_planck;
#endif
#endif
#endif
    
    SphP[i].KappaScatter = 0.0;//0.34;
    
#ifdef BOWEN_OPACITY
    /* Effective transport scattering: changes radiation transport and momentum. */
    SphP[i].KappaScatter += get_bowen_opacity(&SphP[i], allow_radiation_bootstrap);
#endif
}
void update_opacity(int allow_radiation_bootstrap)
{
    int err;
//    struct opal_opacity_result res;
    int idx, i;
#ifdef OPACITY_MESA_DEBUG
    mpi_printf("IDORT: Updating opacity... \n");
#endif
    
#ifdef BOWEN_OPACITY
    int invalid_dust_temperatures = 0;
#pragma omp parallel for private(idx, i) reduction(+:invalid_dust_temperatures)
#else
#pragma omp parallel for private(idx, i)
#endif
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
        set_opacity(SphP, i, allow_radiation_bootstrap);
#ifdef BOWEN_OPACITY
        if(All.BowenTemperatureMode == 1 && !isfinite(SphP[i].KappaScatter))
          invalid_dust_temperatures++;
#endif
    }
#ifdef BOWEN_OPACITY
    if(invalid_dust_temperatures)
      terminate("BOWEN: invalid dust temperature in %d active cells on this task. "
                "Radiation mode requires finite positive E_rad; gas-temperature bootstrap "
                "is allowed only for zero E_rad before a fresh-IC radiation solve.",
                invalid_dust_temperatures);
#endif
#ifdef OPACITY_MESA_DEBUG
    mpi_printf("IDORT: Finish updating opacity... \n");
#endif
}
#ifdef BOWEN_OPACITY
void refresh_bowen_opacity(void)
{
  int invalid_dust_temperatures = 0;
#pragma omp parallel for reduction(+:invalid_dust_temperatures)
  for(int i = 0; i < NumGas; i++)
    {
      SphP[i].KappaScatter = get_bowen_opacity(&SphP[i], 0);
      if(!isfinite(SphP[i].KappaScatter))
        invalid_dust_temperatures++;
      else
        SphP[i].tau = (SphP[i].KappaRoss + SphP[i].KappaScatter) *
            SphP[i].Density * All.UnitDensity_in_cgs * get_cell_radius(i);
    }
  if(invalid_dust_temperatures)
    terminate("BOWEN: cannot refresh changed checkpoint opacity: %d cells have invalid dust temperatures.",
              invalid_dust_temperatures);
}
#endif

int main(void) {
  int i;
  reset();
#ifdef BOWEN_OPACITY
  {
    const double temperatures[] = {50.0, 1000.0, 1100.0, 2000.0, 1e8};
    for(i=0; i<5; i++) {
      cells[0].EOSTemperature = temperatures[i];
      cells[0].E_rad = NAN;
      CHECK(get_bowen_opacity(&cells[0], 0) ==
        All.BowenOpacityKappamax /
        (1.0 + exp((temperatures[i] - All.DustCondensationTemp) / All.BowenOpacitydT)));
    }
    cells[0].Density = All.BowenOpacityRhomin;
    CHECK(get_bowen_opacity(&cells[0], 0) == 0.0);
    cells[0].Density *= 0.5;
    CHECK(get_bowen_opacity(&cells[0], 0) == 0.0);
  }
#ifdef IDORT
  reset();
  All.BowenTemperatureMode = 1;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0), All.BowenOpacityKappamax / 2));
  cells[0].EOSTemperature = 1e7;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0), All.BowenOpacityKappamax / 2));
  CHECK(close_enough(get_bowen_opacity(&cells[0], 1), All.BowenOpacityKappamax / 2));
  cells[0].E_rad = 0.0;
  cells[0].EOSTemperature = 1100.0;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 1), All.BowenOpacityKappamax / 2));
  CHECK(isnan(get_bowen_opacity(&cells[0], 0)));
  {
    const double invalid[] = {-1.0, NAN, INFINITY};
    for(i=0; i<3; i++) {
      cells[0].E_rad = invalid[i];
      CHECK(isnan(get_bowen_opacity(&cells[0], 0)));
      CHECK(isnan(get_bowen_opacity(&cells[0], 1)));
    }
  }
  cells[0].E_rad = RAD_CONST * pow(1e8, 4.0);
  CHECK(get_bowen_opacity(&cells[0], 0) == 0.0);
  cells[0].E_rad = RAD_CONST * pow(1.0, 4.0);
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0),
    All.BowenOpacityKappamax / (1.0 + exp((1.0-1100.0)/50.0))));
  reset();
  All.BowenTemperatureMode = 1;
  update_opacity(0);
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 3); /* Negative active-list entry was skipped. */
  for(i=0; i<3; i++) {
    CHECK(cells[i].KappaRoss == cells[i].Density + cells[i].EOSTemperature * 0.001);
    CHECK(cells[i].KappaPlanck == cells[i].Density + cells[i].EOSTemperature * 0.002);
    CHECK(close_enough(cells[i].KappaScatter, All.BowenOpacityKappamax / 2));
    CHECK(cells[i].tau == -99.0); /* Startup-only tau refresh is in init.c. */
  }
  reset();
  All.BowenTemperatureMode = 1;
  cells[0].E_rad = -1.0;
  cells[1].E_rad = NAN;
  update_opacity(0);
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(calls_at_termination == 3); /* Errors aggregated after all cells. */
  CHECK(!terminated_in_parallel);
  reset();
  All.BowenTemperatureMode = 1;
  All.BowenOpacityKappamax = 12.0;
  TimeBinsHydro.NActiveParticles = 1; /* Only entry is -1; all cells inactive. */
  refresh_bowen_opacity();
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 0);
  for(i=0; i<3; i++) {
    CHECK(close_enough(cells[i].KappaScatter, 6.0));
    CHECK(cells[i].KappaRoss == 100.0 + i);
    CHECK(cells[i].KappaPlanck == 200.0 + i);
    CHECK(cells[i].tau == (cells[i].KappaRoss + cells[i].KappaScatter) *
          cells[i].Density * All.UnitDensity_in_cgs * get_cell_radius(i));
  }
  cells[0].E_rad = -1.0;
  cells[1].E_rad = NAN;
  cells[2].KappaScatter = 999.0;
  refresh_bowen_opacity();
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(close_enough(cells[2].KappaScatter, 6.0));
  CHECK(!terminated_in_parallel);
#endif
  reset();
  All.BowenOpacityKappamax = 17.0;
  TimeBinsHydro.NActiveParticles = 1;
  refresh_bowen_opacity();
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 0);
  for(i=0; i<3; i++) {
    CHECK(cells[i].KappaScatter == 17.0 /
          (1.0 + exp((cells[i].EOSTemperature-1100.0)/50.0)));
    CHECK(cells[i].KappaRoss == 100.0 + i);
    CHECK(cells[i].KappaPlanck == 200.0 + i);
    CHECK(cells[i].tau == (cells[i].KappaRoss + cells[i].KappaScatter) *
          cells[i].Density * All.UnitDensity_in_cgs * get_cell_radius(i));
  }
#endif
  reset();
  update_opacity(0);
  CHECK(terminate_calls == 0);
  CHECK(opal_calls == 3);
  for(i=0; i<3; i++) {
    CHECK(cells[i].tau == -99.0); /* Legacy/non-Bowen update is unchanged. */
    CHECK(cells[i].KappaRoss == cells[i].Density + cells[i].EOSTemperature * 0.001);
    CHECK(cells[i].KappaPlanck == cells[i].Density + cells[i].EOSTemperature * 0.002);
#ifndef BOWEN_OPACITY
    CHECK(cells[i].KappaScatter == 0.0);
#else
    CHECK(cells[i].KappaScatter == All.BowenOpacityKappamax /
          (1.0 + exp((cells[i].EOSTemperature-All.DustCondensationTemp)/All.BowenOpacitydT)));
#endif
  }
  printf("PASS: %d checks\n", checks);
  return 0;
}
