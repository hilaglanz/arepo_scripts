/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/opal_opacity.h
 * \date        08/2023
 * \author      Jing-Ze Ma <jingze@mpa-garching.mpg.de>
 * \brief       source file for MESA opacity
 * \details     implementation of the MESA opacity in C using
 *    interpolation of the MESA tables or freely available data
 *    analogous to OPAL EOS implemented by Sebastian Ohlmann
 *
 *
 * \par Major modifications and contributions:
 *
 * - 10.08.2023 Integration in Arepo
 */

#include "arepoconfig.h"
#include "proto.h"

#ifdef OPACITY_MESA

#undef USE_ZLIB

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef USE_ZLIB
#include <zlib.h>
#endif
#ifdef USE_ZLIB
#define fclose gzclose
#define fgets(x, y, z) gzgets(z, x, y)
#endif

//#include "opal_opacity.h"


//#ifdef OPAL_GLOBAL
static struct opal_opacity_table *opal_table; /**< main opacity table object */
//#endif

/* function declarations */
static inline double interpolate_linear(double x, double *xs, double *ys);
static inline double interpolate_linear_slope(double slope, double *ys);
static inline double interpolate_quadratic(double x, double *xs, double *ys);
static inline double interpolate_quadratic_coeff(double A, double B, double C, double *ys);
static inline void binary_search(double x, double *list, int len, int points, int *low, int *m);
static inline int linear_search_asc(double x, double *list, int len);
int update_opacity_cache(struct opal_opacity_table *opal_table, double rho, double t);

/**
 * @brief initialization routine.
 * allocate table, open data file and read table from data file
 *
 * @param filename the name of the opacity data file
 *
 * @return the opacity table struct
 */
struct opal_opacity_table *opalopacity_init_local(const char *filename)
{
  /* allocate table */
  struct opal_opacity_table *opal_table = (struct opal_opacity_table *)malloc(sizeof(struct opal_opacity_table));

  if(!opal_table)
    {
#ifdef OPACITY_MESA_DEBUG
      fprintf(stderr, "Error allocating memory for opacity tables.\n");
#endif
      return NULL;
    }

  /* allocate cache */
  opal_table->cache = (struct opal_opacity_cache *)malloc(sizeof(struct opal_opacity_cache));

  if(!opal_table->cache)
    {
#ifdef OPACITY_MESA_DEBUG
      fprintf(stderr, "Error allocating memory for opacity cache.\n");
#endif
      return NULL;
    }
  /* initialize with negative values */
    opal_table->cache->rho  = -30.0;
    opal_table->cache->t    = -30.0;
    opal_table->cache->mint = -30.0;

  /* open data file */
#ifndef USE_ZLIB
  FILE *datafile = fopen(filename, "r");
#else
  gzFile datafile = gzopen(filename, "r");
#endif

  if(!datafile)
    {
#ifdef OPACITY_MESA_DEBUG
      fprintf(stderr, "Error opening opacity data file %s.\n", filename);
#endif
      return NULL;
    }

  /* read in opacity data */
  int ilogR, ilogT;
  int numlogR_read, numlogT_read;
  int n;
  int dumi;
  char line[400];
  double dum, minlogR_read, maxlogR_read, minlogT_read, maxlogT_read;
  int status;
//  float tempFloat;

  /* Reading opacity tables from data file */
  if(!fgets(line, sizeof(line), datafile))
    {
#ifdef OPACITY_MESA_DEBUG
      fprintf(stderr, "Error reading line.\n");
#endif
      fclose(datafile);
      return NULL;
    }
  /* 2 empty lines */
  fgets(line, sizeof(line), datafile);
  fgets(line, sizeof(line), datafile);
    
  /* read third line */
  n = sscanf(line, "%d %d %lf %lf %d %lf %lf %d %lf %lf", &dumi, &dumi, &opal_table->x, &opal_table->z, &numlogR_read, &minlogR_read, &maxlogR_read, &numlogT_read, &minlogT_read, &maxlogT_read);
  /* consistency checks */
  if(n < 3 || numlogR_read != numlogR || numlogT_read != numlogT)
    {
#ifdef OPACITY_MESA_DEBUG
      fprintf(stderr, "Error reading X line, read only %d elements.\n", n);
#endif
      fclose(datafile);
      return NULL;
    }
  /*printf("Reading values for X=%f.\n", opal_table->x[ix]); */
    
  /* 2 empty lines */
  fgets(line, sizeof(line), datafile);
  fgets(line, sizeof(line), datafile);
    
  /* read the line of logR */
  if(!fgets(line, sizeof(line), datafile))
    {
#ifdef OPACITY_MESA_DEBUG
      fprintf(stderr, "Error reading logR line.\n");
#endif
      fclose(datafile);
      return NULL;
    }
  double tempDouble;
  int numRead = 0;
  char *ptr = line;
  // Use sscanf to read doubles from the line separated by multiple spaces
  while (ptr != NULL && *ptr == ' ')
    ptr++; // Move past the space(s)
  while (sscanf(ptr, "%lf", &tempDouble) == 1) {
    opal_table->logR[numRead] = tempDouble;
    ptr = strchr(ptr, ' '); // Find the next space-separated value
    while (ptr != NULL && *ptr == ' ')
        ptr++; // Move past the space(s)
    if (ptr == NULL)
        break;
    else
        numRead++;
//    ptr++; // Move past the space
  }
  // Check if we read the expected number of floats
  if (numRead+1 != numlogR)
    {
#ifdef OPACITY_MESA_DEBUG
    fprintf(stderr, "Error reading logR line number. %d not equal to %d. Final number is %lf.\n", numRead+1, numlogR, tempDouble);
#endif
    fclose(datafile);
    return NULL;
    }
  opal_table->minlogR = opal_table->logR[0];
  opal_table->maxlogR = opal_table->logR[numlogR - 1];
  if (opal_table->minlogR != minlogR_read || opal_table->maxlogR != maxlogR_read)
    {
  #ifdef OPACITY_MESA_DEBUG
    fprintf(stderr, "Error reading logR. Reading range %lf %lf not equal to table range %lf %lf.\n", opal_table->minlogR, opal_table->maxlogR, minlogR_read, maxlogR_read);
  #endif
    fclose(datafile);
    return NULL;
    }
  
  /* 1 empty lines */
  fgets(line, sizeof(line), datafile);

  for(ilogT = 0; ilogT < numlogT; ilogT++)
    {
      /* read first line of logT block */
      if(!fgets(line, sizeof(line), datafile))
        {
#ifdef OPAL_DEBUG
          fprintf(stderr, "Error reading logT line.\n");
#endif
          fclose(datafile);
          return NULL;
        }
      double tempDouble;
      int numRead = 0;
      char *ptr = line;
      // Use sscanf to read doubles from the line separated by multiple spaces
      while (ptr != NULL && *ptr == ' ')
        ptr++; // Move past the space(s)
      while (sscanf(ptr, "%lf", &tempDouble) == 1) {
          if (numRead == 0)
              opal_table->logT[ilogT] = tempDouble;
          else
              opal_table->logkappa[ilogT][numRead-1] = tempDouble;
          ptr = strchr(ptr, ' '); // Find the next space-separated value
          while (ptr != NULL && *ptr == ' ')
              ptr++; // Move past the space(s)
          if (ptr == NULL)
              break;
          else
              numRead++;
//          ptr++; // Move past the space
      }
      // Check if we read the expected number of floats
      if (numRead != numlogR)
        {
      #ifdef OPACITY_MESA_DEBUG
        fprintf(stderr, "Error reading logT line number. %d not equal to %d. Final number is %lf.\n", numRead, numlogR, tempDouble);
      #endif
        fclose(datafile);
        return NULL;
        }
    }
    
  opal_table->minlogT = opal_table->logT[0];
  opal_table->maxlogT = opal_table->logT[numlogT - 1];

  /* close data file */
  fclose(datafile);

  printf("Read opacity tables from data file %s for metallicity %.3f.\n", filename, opal_table->z);
//  printf("logR range %lf %lf, logT range %lf %lf.\n", opal_table->minlogR, opal_table->maxlogR, opal_table->minlogT, opal_table->maxlogT);
    
//  opal_table->inter = gsl_interp2d_alloc(gsl_interp2d_bilinear, numlogR, numlogT);
//  opal_table->xacc  = gsl_interp_accel_alloc();
//  opal_table->yacc  = gsl_interp_accel_alloc();
//
//  status = gsl_interp2d_init(opal_table->inter, (double *)opal_table->logR, (double *)opal_table->logT, (double *)opal_table->logkappa, numlogR, numlogT);
    
//    opacityinterp->inter = gsl_interp2d_alloc(gsl_interp2d_bilinear, numlogR, numlogT);
//    opacityinterp->xacc  = gsl_interp_accel_alloc();
//    opacityinterp->yacc  = gsl_interp_accel_alloc();
//
//    if(gsl_interp2d_init(opacityinterp->inter, (double *)opal_table->logR, (double *)opal_table->logT, (double *)opal_table->logkappa, numlogR, numlogT) < 0)
//    {
//  #ifdef OPACITY_MESA_DEBUG
//      fprintf(stderr, "Error interpolating.\n");
//  #endif
//    }
    
    
  return opal_table;
}

/**
 * @brief deinitialization routine.
 *
 * free opacity table data
 */
void opalopacity_deinit_local(struct opal_opacity_table *opal_table)
{
  free(opal_table->cache);
  free(opal_table);
    
//  gsl_interp2d_free(opal_table->inter);
//  gsl_interp_accel_free(opal_table->xacc);
//  gsl_interp_accel_free(opal_table->yacc);
}

//void opacityinterp_deinit_local(opal_opacity_interp *opacityinterp)
//{
//  gsl_interp2d_free(opacityinterp->inter);
//  gsl_interp_accel_free(opacityinterp->xacc);
//  gsl_interp_accel_free(opacityinterp->yacc);
//}

//#ifdef OPAL_GLOBAL
/**
 * @brief Initialization routine if global table is used
 *
 * @param filename Path to table file
 *
 * @return  <0 if errror
 */
int opalopacity_init_global(const char *filename)
{
  /* read in file with task 0 and send to other tasks */
  if(ThisTask == 0)
    {
      opal_table = opalopacity_init_local(filename);
    }
  else
    {
      /* allocate table */
      opal_table = (struct opal_opacity_table *)malloc(sizeof(struct opal_opacity_table));

      if(!opal_table)
        {
#ifdef OPACITY_MESA_DEBUG
          fprintf(stderr, "Error allocating memory for opacity table.\n");
#endif
          return -1;
        }

      /* allocate cache */
      opal_table->cache = (struct opal_opacity_cache *)malloc(sizeof(struct opal_opacity_cache));

      if(!opal_table->cache)
        {
#ifdef OPACITY_MESA_DEBUG
          fprintf(stderr, "Error allocating memory for opacity cache.\n");
#endif
          return -1;
        }
      /* initialize with negative values */
        opal_table->cache->rho  = -30.0;
        opal_table->cache->t    = -30.0;
        opal_table->cache->mint = -30.0;
    }
  /* now broadcast table data */
  MPI_Bcast(&opal_table->x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->logR, numlogR, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->logT, numlogT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(opal_table->logkappa, numlogT * numlogR, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->minlogR, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->maxlogR, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->minlogT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&opal_table->maxlogT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  return 0;
}

/**
 * @brief deinitialization routine for global table
 */
void opalopacity_deinit_global(void) {
    opalopacity_deinit_local(opal_table);
    
}

/**
 * @brief get limits of table in density
 *
 * @param lower lower border
 * @param upper upper border
 */
void opalopacity_get_logT_limits(double *lower, double *upper)
{
  *lower = opal_table->minlogT;
  *upper = opal_table->maxlogT;
}

/**
 * @brief Updates values in cache for faster interpolation
 *
 * @param opal_table the opal table struct
 * @param rho,t opacity parameters
 *
 * @return -1 if error occurred, 0 otherwise
 */
int update_opacity_cache(struct opal_opacity_table *opal_table, double rho, double t)
{
  /* check if update needed */
  if(fabs(opal_table->cache->rho - rho) < 1e-15 && fabs(opal_table->cache->t - t) < 1e-15 &&
     opal_table->cache->valid)
    return 0;

  /* set to non valid first */
  opal_table->cache->valid = 0;

  int num_points = 0;
//#if INTERPOLATION_METHOD == IP_LINEAR
  num_points = 2;

  if(fabs(opal_table->cache->rho - rho) > 1e-15)
    {
      if((rho < opal_table->minlogR) || (rho > opal_table->maxlogR))
        {
#ifdef OPACITY_MESA_DEBUG
          fprintf(stderr, "Error: logR out of range. Value: %e, limits: %e to %e.\n", rho, opal_table->minlogR, opal_table->maxlogR);
#endif
          return -1;
        }
      /* table is in ascending order for rho */
      binary_search(rho, opal_table->logR, numlogR, num_points, &opal_table->cache->irho, &opal_table->cache->mrho);
      /*printf("irho: %d, rho: %f %f %f\n", irho, rho, opal_table->rho[irho], opal_table->rho[irho+1]); */
    }

  /* check also if temperature limits changed with rho! */
  double mint = opal_table->minlogT;
  // double mint = opal_table->mint[opal_table->cache->mrho];
  if((fabs(opal_table->cache->t - t) > 1e-15) || (fabs(opal_table->cache->mint - mint) > 1e-15))
    {
      if((t < mint) || (t > opal_table->maxlogT))
        {
#ifdef OPACITY_MESA_DEBUG
          fprintf(stderr,
                  "Error: logT out of range. Value: %e, limits: %e to %e \
for rho=%e.\n",
                  t, mint, opal_table->maxlogT, rho);
#endif
          opal_table->cache->mint = mint;
          return -1;
        }
      /* table is in descending order for T */
      int minnum = numlogT;
      // int minnum = opal_table->numt[opal_table->cache->mrho];
      binary_search(t, opal_table->logT, minnum, num_points, &opal_table->cache->it, &opal_table->cache->mt);
      if(fabs(t - opal_table->maxlogT) < 1e-15)
        opal_table->cache->it = opal_table->cache->mt = 0;
    }

    /* pre-compute interpolation quantities */
//#if INTERPOLATION_METHOD == IP_LINEAR
  /* set slopes */
  
  if(fabs(opal_table->cache->rho - rho) > 1e-15)
    opal_table->cache->sloperho = (rho - opal_table->logR[opal_table->cache->irho]) /
                                  (opal_table->logR[opal_table->cache->irho + 1] - opal_table->logR[opal_table->cache->irho]);
  if((fabs(opal_table->cache->t - t) > 1e-15) || (fabs(opal_table->cache->mint - mint) > 1e-15))
    opal_table->cache->slopet =
        (t - opal_table->logT[opal_table->cache->it]) / (opal_table->logT[opal_table->cache->it + 1] - opal_table->logT[opal_table->cache->it]);


  opal_table->cache->rho  = rho;
  opal_table->cache->t    = t;
  opal_table->cache->mint = mint;

  /* for debugging
     printf("irho: %d, mrho: %d; it: %d, mt: %d; numt[mrho]: %d, numt[mrho+1]: %d, numt[mrho+2]: %d;\n",
     opal_table->cache->irho, opal_table->cache->mrho,
     opal_table->cache->it, opal_table->cache->mt,
     opal_table->numt[opal_table->cache->mrho],
     opal_table->numt[opal_table->cache->mrho + 1],
     opal_table->numt[opal_table->cache->mrho + 2]);
     printf("t: %e, t[mt]: %e, t[mt+1]: %e, t[mt+2]: %e, mint[mrho+1]: %e\n",
     opal_table->cache->t, opal_table->t[opal_table->cache->mt],
     opal_table->t[opal_table->cache->mt + 1],
     opal_table->t[opal_table->cache->mt + 2],
     opal_table->mint[opal_table->cache->mrho + 1]);
   */

  /* everything fine, set to valid */
  opal_table->cache->valid = 1;
  return 0;
}


/**
 * @brief retrieve value of certain variable
 *
 * @param opal_table the opal table struct
 * @param var opacity variable to interpolate
 * @param rho density in g/cm^3
 * @param t temperature in K
 * @param result interpolated value
 *
 * @result negative if error occured, otherwise zero
 */
int get_opacitytable_value(struct opal_opacity_table *opal_table, opacitytable var, double rho, double t, double *result)
{
  int i,j;

  /* update cache */
  if(update_opacity_cache(opal_table, rho, t) < 0)
    return -1;

//#if INTERPOLATION_METHOD == IP_LINEAR
//  /* interpolate in x first */
  double varrt[2][2];
  /* use cached slopes for faster interpolation */
  for(i = 0; i < 2; i++)
    for(j = 0; j < 2; j++)
//      varrt[i][j] = interpolate_linear_slope(opal_table->cache->slopex,
//                                             var[opal_table->cache->it + i][opal_table->cache->irho + j] + opal_table->cache->ix);
        varrt[i][j] = var[opal_table->cache->it + i][opal_table->cache->irho + j];

  /* interpolate in rho */
  double vart[2];
  for(i = 0; i < 2; i++)
    vart[i] = interpolate_linear_slope(opal_table->cache->sloperho, varrt[i]);

  /* interpolate in T */
  //*result = interpolate_linear(t, opal_table->t+it, vart);
  *result = interpolate_linear_slope(opal_table->cache->slopet, vart);

  return 0;
}


/**
 * @brief linear interpolation
 *
 * @param x value to interpolate at
 * @param xs pointer to array x values for interpolation
 * @param ys pointer to array y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_linear(double x, double *xs, double *ys) { return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0]); }

/**
 * @brief linear interpolation using given slope
 *
 * @param slope given slope
 * @param ys pointer to array y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_linear_slope(double slope, double *ys) { return ys[0] + slope * (ys[1] - ys[0]); }


/**
 * @brief quadratic interpolation
 *
 * @param x value to interpolate at
 * @param xs pointer to array x values for interpolation
 * @param ys pointer to array y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_quadratic(double x, double *xs, double *ys)
{
  double A = (x - xs[1]) * (x - xs[2]) / (xs[0] - xs[1]) / (xs[0] - xs[2]);
  double B = (x - xs[0]) * (x - xs[2]) / (xs[1] - xs[0]) / (xs[1] - xs[2]);
  double C = (x - xs[0]) * (x - xs[1]) / (xs[2] - xs[0]) / (xs[2] - xs[1]);
  return A * ys[0] + B * ys[1] + C * ys[2];
}

/**
 * @brief quadratic interpolation using given coefficients
 *
 * @param A,B,C given coefficients
 * @param ys y values for interpolation
 *
 * @return interpolated value
 */
double interpolate_quadratic_coeff(double A, double B, double C, double *ys) { return A * ys[0] + B * ys[1] + C * ys[2]; }

/**
 * @brief binary search in lists in arbitrary order
 *
 * @param x value to search for
 * @param list list to search in
 * @param len length of list
 * @param points number of points for interpolating polynomials
 *        (2 for linear, 3 for quadratic etc)
 * @param low lowest index in array with x centered in a subrange of points points
 * @param m index with list[index] <= x < list[index+1] (ascending)
 *  or with list[index] < x <= list[index+1] (descending)
 */
void binary_search(double x, double *list, int len, int points, int *low, int *m)
{
  int ilow = 0, ihigh = len - 1, imid;
  int ascend = (list[0] <= list[len - 1]);
  while(ilow < ihigh - 1)
    {
      /* imid = (ilow + ihigh) / 2; */
      imid = (ilow + ihigh) >> 1;
      if((list[imid] <= x) == ascend)
        ilow = imid;
      else
        ihigh = imid;
    }
  /* now ilow is index with list[index] <= x < list[index+1] (ascending)
   *  or with list[index] < x <= list[index+1] (descending) */
  *m = ilow;

  /* low: min 0, max len-points, else m-floor((points-2)/2) */
  *low = ilow - ((points - 2) >> 1);
  *low = fmax(0, fmin(len - points, *low));
}

/**
 * @brief linear search in lists in ascending order
 *
 * @param x value to search for
 * @param list list to search in
 * @param len length of list
 *
 * @return index with list[index] <= x < list[index+1]
 */
int linear_search_asc(double x, double *list, int len)
{
  int ilow = 0;
  while((list[ilow] > x) && (ilow < len))
    ilow++;
  if(ilow < len)
    return ilow;
  else
    return -1;
}


/**
 * @brief main driver routine for the OPAL opacity
 *
 * @param opal_table the opal table struct
 * @param t temperature in K
 * @param rho density in g/cm^3
 * @param result pointer to double array where results are stored
 *
 * @result negative if error occured, otherwise zero
 */

int get_opal_opacity_local(struct opal_opacity_table *opal_table, double rho, double t, struct opal_opacity_result *result)
{
  int err;
  double res;
  int status;
  double logt, t6, logr;
  err = 0;

  logt = log10(t);
  t6   = t / 1e6;
  logr = log10(rho) - 3 * log10(t6);

  if(logr <= opal_table->logR[0] || logr >= opal_table->logR[numlogR - 1] || logt <= opal_table->logT[0] || logt >= opal_table->logT[numlogT - 1])
      res = 1e-30;
  else
    {
//      err    = gsl_interp2d_eval_e(opacityinterp->inter, (double *)opal_table->logR, (double *)opal_table->logT, (double *)opal_table->logkappa, logr, logt,
//                                   opacityinterp->xacc, opacityinterp->yacc, &result);
        
        /* update cache and return with error if value not in table */
      err = update_opacity_cache(opal_table, logr, logt);
      err += get_opacitytable_value(opal_table, opal_table->logkappa, logr, logt, &res);
      res = pow(10, res);
    }

  if(res < 1e-5)
      res = 1e-5;
  result->kappa_ross = res;
    
//  opacityinterp_deinit_local();
//#ifdef OPACITY_MESA_DEBUG
//  mpi_printf("IDORT: Updating opacity... %lf\n",result);
//#endif
  /* check for errors */
  if(err < 0)
    {
      fprintf(stderr, "Error obtaining opacity values for rho=%e, T=%e.\n", rho, t);
      return -1;
    }
  else
      return 0;
}

int get_opal_opacity_global(double rho, double t, struct opal_opacity_result *result)
{
    return get_opal_opacity_local(opal_table, rho, t, result);
}

#ifdef BOWEN_OPACITY
/* Returns NaN on invalid input so errors can be reported outside OpenMP. */
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

/* Call outside OpenMP. Zero-radiation reports are limited to one per MPI rank;
 * invalid inputs are always reported before termination. No checkpoint state. */
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

#ifdef BOWEN_OPACITY
void set_opacity(struct sph_particle_data *SphP, int i, int allow_radiation_bootstrap)
#else
void set_opacity(struct sph_particle_data *SphP, int i)
#endif
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

#ifdef BOWEN_OPACITY
/* Bootstrap is allowed only before the first radiation solve of a fresh IC. */
void update_opacity(int allow_radiation_bootstrap)
#else
void update_opacity(void)
#endif
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
#ifdef BOWEN_OPACITY
        set_opacity(SphP, i, allow_radiation_bootstrap);
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
#else
        set_opacity(SphP, i);
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
/* A changed checkpoint prescription must also update inactive cells before
 * the next RT solve. Baseline table opacities and gas state are untouched. */
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

//void IDORT(void)
//{
//#ifdef OPACITY_MESA_DEBUG
//  mpi_printf("IDORT: Start IDORT... \n");
//#endif
//  TIMER_START(CPU_IDORT);
//    
//  update_opacity();
//
//  TIMER_STOP(CPU_IDORT);
//    
////  opalopacity_deinit();
//    
//#ifdef OPACITY_MESA_DEBUG
//    mpi_printf("IDORT: Finish IDORT... \n");
//#endif
//}


#endif /* OPACITY_MESA */
