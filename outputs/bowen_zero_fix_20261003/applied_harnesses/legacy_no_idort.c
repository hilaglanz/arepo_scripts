#define OPACITY_PLANCK
#define BOWEN_OPACITY

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define RAD_CONST 7.5657e-15
#define IDORT_NUM_ANGLES 4
struct sph_particle_data {
  double Density, EOSTemperature, E_rad;
  double KappaRoss, KappaPlanck, KappaScatter, tau;
  double Volume, Intensity[IDORT_NUM_ANGLES], IV[IDORT_NUM_ANGLES];
};
static struct particle_data {
  unsigned long long ID;
  int Type;
  double Mass, Pos[3];
} P[3];
static int ThisTask = 29, RestartFlag = 1;
struct {
  int BowenTemperatureMode;
  double BowenOpacityRhomin, BowenOpacityKappamax;
  double DustCondensationTemp, BowenOpacitydT, UnitDensity_in_cgs, Time;
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
static char diagnostic_text[8192], terminate_text[1024];
static int diagnostic_calls, diagnostic_in_parallel;
static int capture_fprintf(FILE *stream, const char *format, ...) {
  va_list arguments;
  const size_t used = strlen(diagnostic_text);
  int written;
  (void)stream;
#ifdef _OPENMP
  diagnostic_in_parallel |= omp_in_parallel();
#endif
  diagnostic_calls++;
  va_start(arguments, format);
  written = vsnprintf(diagnostic_text + used, sizeof(diagnostic_text) - used,
                      format, arguments);
  va_end(arguments);
  return written;
}
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
  va_start(arguments, format);
  vsnprintf(terminate_text, sizeof(terminate_text), format, arguments);
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
  All.Time = 200050555.8;
  opal_calls = terminate_calls = invalid_count = terminated_in_parallel = 0;
  calls_at_termination = 0;
  diagnostic_calls = diagnostic_in_parallel = 0;
  diagnostic_text[0] = terminate_text[0] = '\0';
  TimeBinsHydro.NActiveParticles = 4;
  for(i=0; i<3; i++) {
    cells[i].Density = 1e-12;
    cells[i].EOSTemperature = 2000.0;
    cells[i].E_rad = RAD_CONST * pow(1100.0, 4.0);
    cells[i].KappaScatter = 999.0;
    cells[i].KappaRoss = 100.0 + i;
    cells[i].KappaPlanck = 200.0 + i;
    cells[i].tau = -99.0;
    cells[i].Volume = 16.0 + i;
    P[i].ID = 10000000000ULL + i;
    P[i].Type = 0;
    P[i].Mass = 123.0 + i;
    P[i].Pos[0] = 1.0 + i;
    P[i].Pos[1] = 2.0 + i;
    P[i].Pos[2] = 3.0 + i;
    for(int angle=0; angle<IDORT_NUM_ANGLES; angle++) {
      cells[i].Intensity[angle] = 1.0;
      cells[i].IV[angle] = cells[i].Volume;
    }
  }
}
#define fprintf capture_fprintf
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
  if(!isfinite(radiation_energy) || radiation_energy < 0.0)
    return NAN;
  if(allow_radiation_bootstrap && radiation_energy == 0.0)
    {
      /* Only fresh ICs may use gas temperature before radiation is initialized. */
      if(!isfinite(temperature_for_dust) || temperature_for_dust <= 0.0)
        return NAN;
    }
  else
    {
      /* E_rad is an energy density in physical CGS units. Zero radiation is
       * allowed (e.g. IDORT_OUTFLOW) and gives zero grey temperature.
       * Separate roots avoid overflow in E_rad / RAD_CONST. */
      temperature_for_dust = sqrt(sqrt(radiation_energy)) / sqrt(sqrt(RAD_CONST));
    }
#else
  /* check_parameters() rejects radiation mode without IDORT. */
  return NAN;
#endif
  if(!isfinite(temperature_for_dust) || temperature_for_dust < 0.0)
    return NAN;

  const double temperature_offset =
      (temperature_for_dust - All.DustCondensationTemp) / All.BowenOpacitydT;
  const double exponential = exp(-fabs(temperature_offset));
  const double condensation_fraction = temperature_offset >= 0.0
      ? exponential / (1.0 + exponential) : 1.0 / (1.0 + exponential);
  return All.BowenOpacityKappamax * condensation_fraction;
}
static void report_bowen_cell(const char *context, int i, int count, int allow_radiation_bootstrap, int zero_radiation)
{
  static int reported_zero_radiation = 0;
  if(zero_radiation && reported_zero_radiation)
    return;
  if(zero_radiation)
    reported_zero_radiation = 1;

  fprintf(stderr, "BOWEN DIAGNOSTIC: %s; context=%s task=%d local_count=%d Time=%.17g "
          "RestartFlag=%d bootstrap=%d i=%d ID=%llu type=%d mass=%.17g "
          "pos=(%.17g,%.17g,%.17g) rho=%.17g Tgas=%.17g KappaScatter=%.17g\n",
          zero_radiation ? "finite zero radiation accepted; first report on this rank" : "invalid dust temperature",
          context, ThisTask, count, All.Time, RestartFlag, allow_radiation_bootstrap, i,
          (unsigned long long)P[i].ID, P[i].Type, P[i].Mass,
          P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].Density, SphP[i].EOSTemperature, SphP[i].KappaScatter);
#ifdef IDORT
  double intensity_min = INFINITY, intensity_max = -INFINITY;
  int zero_angles = 0, negative_angles = 0, nonfinite_angles = 0, nonfinite_iv = 0;
  for(int angle = 0; angle < IDORT_NUM_ANGLES; angle++)
    {
      const double intensity = SphP[i].Intensity[angle];
      if(!isfinite(intensity))
        nonfinite_angles++;
      else
        {
          intensity_min = fmin(intensity_min, intensity);
          intensity_max = fmax(intensity_max, intensity);
          zero_angles += (intensity == 0.0);
          negative_angles += (intensity < 0.0);
        }
      nonfinite_iv += !isfinite(SphP[i].IV[angle]);
    }
  if(nonfinite_angles == IDORT_NUM_ANGLES)
    intensity_min = intensity_max = NAN;
  fprintf(stderr, "BOWEN RADIATION: task=%d ID=%llu E_rad=%.17g Volume=%.17g "
          "stored_I_min_finite=%.17g stored_I_max_finite=%.17g angles=%d "
          "zero_angles=%d negative_angles=%d nonfinite_angles=%d nonfinite_IV=%d\n",
          ThisTask, (unsigned long long)P[i].ID, SphP[i].E_rad, SphP[i].Volume,
          intensity_min, intensity_max, IDORT_NUM_ANGLES,
          zero_angles, negative_angles, nonfinite_angles, nonfinite_iv);
#endif
  fflush(stderr);
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
    int invalid_dust_temperatures = 0, zero_radiation_cells = 0;
    int first_invalid = NumGas, first_zero = NumGas;
#pragma omp parallel for private(idx, i) reduction(+:invalid_dust_temperatures,zero_radiation_cells) reduction(min:first_invalid,first_zero)
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
          {
            invalid_dust_temperatures++;
            first_invalid = i < first_invalid ? i : first_invalid;
          }
#ifdef IDORT
        else if(All.BowenTemperatureMode == 1 && !allow_radiation_bootstrap &&
                SphP[i].Density > All.BowenOpacityRhomin && SphP[i].E_rad == 0.0)
          {
            zero_radiation_cells++;
            first_zero = i < first_zero ? i : first_zero;
          }
#endif
#endif
    }
#ifdef BOWEN_OPACITY
    if(invalid_dust_temperatures)
      {
        report_bowen_cell("update_opacity", first_invalid, invalid_dust_temperatures, allow_radiation_bootstrap, 0);
        terminate("BOWEN: invalid dust temperature in %d active cells on this task; see BOWEN DIAGNOSTIC. "
                  "Radiation must be finite and nonnegative. Fresh-IC bootstrap requires positive finite gas temperature.",
                  invalid_dust_temperatures);
      }
    if(zero_radiation_cells)
      report_bowen_cell("update_opacity", first_zero, zero_radiation_cells, allow_radiation_bootstrap, 1);
#endif
#ifdef OPACITY_MESA_DEBUG
    mpi_printf("IDORT: Finish updating opacity... \n");
#endif
}
#ifdef BOWEN_OPACITY
void refresh_bowen_opacity(void)
{
  int invalid_dust_temperatures = 0, zero_radiation_cells = 0;
  int first_invalid = NumGas, first_zero = NumGas;
#pragma omp parallel for reduction(+:invalid_dust_temperatures,zero_radiation_cells) reduction(min:first_invalid,first_zero)
  for(int i = 0; i < NumGas; i++)
    {
      SphP[i].KappaScatter = get_bowen_opacity(&SphP[i], 0);
      if(!isfinite(SphP[i].KappaScatter))
        {
          invalid_dust_temperatures++;
          first_invalid = i < first_invalid ? i : first_invalid;
        }
      else
        {
          SphP[i].tau = (SphP[i].KappaRoss + SphP[i].KappaScatter) *
              SphP[i].Density * All.UnitDensity_in_cgs * get_cell_radius(i);
#ifdef IDORT
          if(All.BowenTemperatureMode == 1 && SphP[i].Density > All.BowenOpacityRhomin && SphP[i].E_rad == 0.0)
            {
              zero_radiation_cells++;
              first_zero = i < first_zero ? i : first_zero;
            }
#endif
        }
    }
  if(invalid_dust_temperatures)
    {
      report_bowen_cell("refresh_bowen_opacity", first_invalid, invalid_dust_temperatures, 0, 0);
      terminate("BOWEN: cannot refresh changed checkpoint opacity: %d cells have invalid dust temperatures; "
                "see BOWEN DIAGNOSTIC.", invalid_dust_temperatures);
    }
  if(zero_radiation_cells)
    report_bowen_cell("refresh_bowen_opacity", first_zero, zero_radiation_cells, 0, 1);
}
#endif
#undef fprintf

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
  cells[0].EOSTemperature = NAN;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0), All.BowenOpacityKappamax / 2));
  cells[0].E_rad = 0.0;
  cells[0].EOSTemperature = 1100.0;
  CHECK(close_enough(get_bowen_opacity(&cells[0], 1), All.BowenOpacityKappamax / 2));
  CHECK(close_enough(get_bowen_opacity(&cells[0], 0),
    All.BowenOpacityKappamax / (1.0 + exp(-1100.0/50.0))));
  {
    const double invalid_gas[] = {0.0, -1.0, NAN, INFINITY};
    for(i=0; i<4; i++) {
      cells[0].EOSTemperature = invalid_gas[i];
      CHECK(isnan(get_bowen_opacity(&cells[0], 1)));
      CHECK(close_enough(get_bowen_opacity(&cells[0], 0),
        All.BowenOpacityKappamax / (1.0 + exp(-1100.0/50.0))));
    }
  }
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
  cells[0].E_rad = DBL_MAX; /* E_rad/RAD_CONST would overflow. */
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
  cells[0].Intensity[0] = 0.0;
  cells[0].Intensity[1] = -2.0;
  cells[0].Intensity[2] = NAN;
  cells[0].Intensity[3] = 3.0;
  cells[0].IV[0] = INFINITY;
  update_opacity(0);
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(calls_at_termination == 3); /* Errors aggregated after all cells. */
  CHECK(!terminated_in_parallel);
  CHECK(!diagnostic_in_parallel);
  CHECK(diagnostic_calls == 2);
  CHECK(strstr(diagnostic_text, "context=update_opacity task=29 local_count=2") != NULL);
  CHECK(strstr(diagnostic_text, "RestartFlag=1 bootstrap=0 i=0 ID=10000000000 type=0 mass=123") != NULL);
  CHECK(strstr(diagnostic_text, "pos=(1,2,3)") != NULL);
  CHECK(strstr(diagnostic_text, "rho=") != NULL && strstr(diagnostic_text, "Tgas=2000") != NULL);
  CHECK(strstr(diagnostic_text, "E_rad=-1 Volume=16") != NULL);
  CHECK(strstr(diagnostic_text, "stored_I_min_finite=-2 stored_I_max_finite=3 angles=4") != NULL);
  CHECK(strstr(diagnostic_text, "zero_angles=1 negative_angles=1 nonfinite_angles=1 nonfinite_IV=1") != NULL);
  CHECK(strstr(terminate_text, "finite and nonnegative") != NULL);
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
  for(i=0; i<IDORT_NUM_ANGLES; i++) cells[0].Intensity[i] = NAN;
  refresh_bowen_opacity();
  CHECK(terminate_calls == 1);
  CHECK(invalid_count == 2);
  CHECK(close_enough(cells[2].KappaScatter, 6.0));
  CHECK(!terminated_in_parallel);
  CHECK(!diagnostic_in_parallel);
  CHECK(diagnostic_calls == 2);
  CHECK(strstr(diagnostic_text, "context=refresh_bowen_opacity task=29 local_count=2") != NULL);
  CHECK(strstr(diagnostic_text, "i=0 ID=10000000000") != NULL);
  CHECK(strstr(diagnostic_text, "E_rad=-1 Volume=16") != NULL);
  CHECK(strstr(diagnostic_text, "zero_angles=0 negative_angles=0 nonfinite_angles=4 nonfinite_IV=0") != NULL);

  /* Excluded/bootstrap zeros must not consume the one-time warning. */
  reset();
  All.BowenTemperatureMode = 1;
  cells[0].E_rad = 0.0;
  cells[0].Density = All.BowenOpacityRhomin;
  update_opacity(0);
  CHECK(terminate_calls == 0 && diagnostic_calls == 0);
  CHECK(cells[0].KappaScatter == 0.0);
  cells[0].Density = 1e-12;
  cells[0].EOSTemperature = 1100.0;
  update_opacity(1);
  CHECK(terminate_calls == 0 && diagnostic_calls == 0);
  CHECK(close_enough(cells[0].KappaScatter, All.BowenOpacityKappamax / 2));
  All.BowenTemperatureMode = 0;
  update_opacity(0);
  CHECK(terminate_calls == 0 && diagnostic_calls == 0);

  reset();
  All.BowenTemperatureMode = 1;
  cells[1].E_rad = cells[2].E_rad = 0.0;
  for(i=0; i<IDORT_NUM_ANGLES; i++) cells[1].Intensity[i] = 0.0;
  update_opacity(0);
  CHECK(terminate_calls == 0 && diagnostic_calls == 2);
  CHECK(!diagnostic_in_parallel);
  CHECK(strstr(diagnostic_text, "finite zero radiation accepted; first report on this rank") != NULL);
  CHECK(strstr(diagnostic_text, "context=update_opacity task=29 local_count=2") != NULL);
  CHECK(strstr(diagnostic_text, "i=1 ID=10000000001") != NULL);
  CHECK(strstr(diagnostic_text, "E_rad=0 Volume=17") != NULL);
  CHECK(strstr(diagnostic_text, "zero_angles=4 negative_angles=0 nonfinite_angles=0") != NULL);
  update_opacity(0);
  CHECK(diagnostic_calls == 2); /* No duplicate warning on subsequent updates. */
  refresh_bowen_opacity();
  CHECK(diagnostic_calls == 2); /* Shared suppression across both callers. */
  CHECK(terminate_calls == 0);

  /* Invalid radiation remains diagnostic and fatal after a zero warning. */
  reset();
  All.BowenTemperatureMode = 1;
  cells[1].E_rad = INFINITY;
  update_opacity(0);
  CHECK(terminate_calls == 1 && invalid_count == 1 && diagnostic_calls == 2);
  CHECK(strstr(diagnostic_text, "invalid dust temperature") != NULL);
  CHECK(strstr(diagnostic_text, "i=1 ID=10000000001") != NULL);
  CHECK(!diagnostic_in_parallel && !terminated_in_parallel);
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
