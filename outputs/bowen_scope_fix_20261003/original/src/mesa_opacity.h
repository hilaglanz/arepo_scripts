/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/opal_opacity.h
 * \date        08/2023
 * \author      Jing-Ze Ma <jingze@mpa-garching.mpg.de>
 * \brief       header file for MESA opacity
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

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>

#ifdef OPACITY_MESA

/* #define OPACITY_MESA_DEBUG */
//#undef OPACITY_MESA_DEBUG

/* opacity table definitions */
#define numlogR 37
#define numlogT 180

//#define IP_LINEAR 1
//#define IP_QUADRATIC 2
//#define INTERPOLATION_METHOD IP_QUADRATIC
//#define INTERPOLATION_METHOD IP_LINEAR

/* global table object? */
//#define OPAL_GLOBAL

/* maximum number of Newton iterations */
//#define OPAL_NEWTON_MAX 70
//#define OPAL_NEWTON_ACCURACY 1e-10

typedef double opacitytable[numlogT][numlogR]; /**< table for opacity variables */

//typedef struct
//{
//  gsl_interp2d *inter;
//  gsl_interp_accel *xacc;
//  gsl_interp_accel *yacc;
//} opal_opacity_interp;

/** @brief opacity result struct.
 * It contains variables for all opacity variables.
 */
struct opal_opacity_result
{
  double kappa_ross;           /**< Rosseland opacity [cm^2/g]                  */
};

/**
 * @brief cache structure for faster computation
 */
struct opal_opacity_cache
{
  double rho, t; /**< values for evaluation */
  int irho, it; /**< corresponding low indices */
  int mrho, mt; /**< indices directly near values */
//#if INTERPOLATION_METHOD == IP_LINEAR
  double sloperho, slopet; /**< pre-computed slopes for linear interpolation */
//#elif INTERPOLATION_METHOD == IP_QUADRATIC
//  double Arho1, Brho1, Crho1, At1, Bt1, Ct1; /**< pre-computed constants for quadratic interpolation 1 */
//  double Arho2, Brho2, Crho2, At2, Bt2, Ct2; /**< pre-computed constants for quadratic interpolation 2 */
//  double sloperho, slopet;                          /**< pre-computed slopes for overlapping quadratics */
//#endif
  double mint; /**< minimum T for current rho */
  int valid;   /**< flag, set to 0 if current values in range of table */
};

/**
 * @brief Main opacity struct.
 * It contains the data tables and the logR and logT values.
 *
 * opacity tables have the following dimensions:
 * - 37 points in logR
 * - 138 points in logT
 */
/* data structure for opal opacities */
struct opal_opacity_table
{
  opacitytable logkappa;
  double x;                    /**< X values (H mass fraction)*/
  double z;                    /**< metallicity*/
  double logR[numlogR];
  double logT[numlogT];
  double minlogR;              /**< minimum logR*/
  double maxlogR;              /**< maximum logR */
  double minlogT;              /**< minimum logT */
  double maxlogT;              /**< maximum logT */
    
//  gsl_interp2d *inter;
//  gsl_interp_accel *xacc;
//  gsl_interp_accel *yacc;

  struct opal_opacity_cache *cache; /**< opacity cache object */
};

/* function declarations */
struct opal_opacity_table *opalopacity_init_local(const char *filename);
void opalopacity_deinit_local(struct opal_opacity_table *opal_table);
//void opacityinterp_deinit_local(opal_opacity_interp *opacityinterp);
//#ifdef OPAL_GLOBAL
int opalopacity_init_global(const char *filename);
void opalopacity_deinit_global(void);
void opalopacity_get_logT_limits(double *lower, double *upper);
int get_opacitytable_value(struct opal_opacity_table *opal_table, opacitytable var, double rho, double t, double *result);
int get_opal_opacity_local(struct opal_opacity_table *opal_table, double rho, double t, struct opal_opacity_result *result);
int get_opal_opacity_global(double rho, double t, struct opal_opacity_result *result);
void set_opacity(struct sph_particle_data *SphP, int i, int allow_radiation_bootstrap);
void update_opacity(int allow_radiation_bootstrap);
#ifdef BOWEN_OPACITY
void refresh_bowen_opacity(void);
#endif
//void IDORT(void);
//#endif
//int opaleos(struct opal_eos_table *opal_table, double x, double rho, double t, int mode, struct opal_eos_result *result);
//int opaleos_all_tgiven_local(struct opal_eos_table *opal_table, double x, double rho, double t, struct opal_eos_result *result);
//int opaleos_all_egiven_local(struct opal_eos_table *opal_table, double x, double rho, double e, double *tguess,
//                             struct opal_eos_result *result);
//int opaleos_egiven_local(struct opal_eos_table *opal_table, double x, double rho, double e, double *tguess, int mode,
//                         struct opal_eos_result *result);
//int opaleos_all_pgiven_local(struct opal_eos_table *opal_table, double x, double rho, double p, double *tguess,
//                             struct opal_eos_result *result);
//int opaleos_pgiven_local(struct opal_eos_table *opal_table, double x, double rho, double p, double *tguess, int mode,
//                         struct opal_eos_result *result);
//#ifdef OPAL_GLOBAL
//int opaleos_all_tgiven_global(double x, double rho, double t, struct opal_eos_result *result);
//int opaleos_all_egiven_global(double x, double rho, double e, double *tguess, struct opal_eos_result *result);
//int opaleos_egiven_global(double x, double rho, double e, double *tguess, int mode, struct opal_eos_result *result);
//int opaleos_all_pgiven_global(double x, double rho, double p, double *tguess, struct opal_eos_result *result);
//int opaleos_pgiven_global(double x, double rho, double p, double *tguess, int mode, struct opal_eos_result *result);
//#endif
//int get_eostable_value(struct opal_eos_table *opal_table, eostable var, double x, double rho, double t, double *result);

//#ifdef OPAL_GLOBAL
#define opalopacity_init opalopacity_init_global
#define opalopacity_deinit opalopacity_deinit_global
#define get_opal_opacity get_opal_opacity_global
//#define opaleos_all_tgiven opaleos_all_tgiven_global
//#define opaleos_all_egiven opaleos_all_egiven_global
//#define opaleos_egiven opaleos_egiven_global
//#define opaleos_all_pgiven opaleos_all_pgiven_global
//#define opaleos_pgiven opaleos_pgiven_global
//#else
//#define opaleos_init opaleos_init_local
//#define opaleos_deinit opaleos_deinit_local
//#define opaleos_all_tgiven opaleos_all_tgiven_local
//#define opaleos_all_egiven opaleos_all_egiven_local
//#define opaleos_egiven opaleos_egiven_local
//#define opaleos_all_pgiven opaleos_all_pgiven_local
//#define opaleos_pgiven opaleos_pgiven_local
//#endif

#ifdef OPACITY_PLANCK

typedef double planckopacitytable[numlogT][numlogR]; /**< table for opacity variables */

struct planck_opacity_result
{
  double kappa_planck;           /**< Planck opacity [cm^2/g]                  */
};

/**
 * @brief cache structure for faster computation
 */
struct planck_opacity_cache
{
  double rho, t; /**< values for evaluation */
  int irho, it; /**< corresponding low indices */
  int mrho, mt; /**< indices directly near values */
//#if INTERPOLATION_METHOD == IP_LINEAR
  double sloperho, slopet; /**< pre-computed slopes for linear interpolation */
//#elif INTERPOLATION_METHOD == IP_QUADRATIC
//  double Arho1, Brho1, Crho1, At1, Bt1, Ct1; /**< pre-computed constants for quadratic interpolation 1 */
//  double Arho2, Brho2, Crho2, At2, Bt2, Ct2; /**< pre-computed constants for quadratic interpolation 2 */
//  double sloperho, slopet;                          /**< pre-computed slopes for overlapping quadratics */
//#endif
  double mint; /**< minimum T for current rho */
  int valid;   /**< flag, set to 0 if current values in range of table */
};

/**
 * @brief Main opacity struct.
 * It contains the data tables and the logR and logT values.
 *
 * opacity tables have the following dimensions:
 * - 37 points in logR
 * - 138 points in logT
 */
/* data structure for opal opacities */
struct planck_opacity_table
{
  planckopacitytable logkappa;
  double x;                    /**< X values (H mass fraction)*/
  double z;                    /**< metallicity*/
  double logR[numlogR];
  double logT[numlogT];
  double minlogR;              /**< minimum logR*/
  double maxlogR;              /**< maximum logR */
  double minlogT;              /**< minimum logT */
  double maxlogT;              /**< maximum logT */
    
//  gsl_interp2d *inter;
//  gsl_interp_accel *xacc;
//  gsl_interp_accel *yacc;

  struct planck_opacity_cache *cache; /**< opacity cache object */
};

/* function declarations */
struct planck_opacity_table *planckopacity_init_local(const char *filename);
void planckopacity_deinit_local(struct planck_opacity_table *planck_table);
//void opacityinterp_deinit_local(opal_opacity_interp *opacityinterp);
//#ifdef OPAL_GLOBAL
int planckopacity_init_global(const char *filename);
void planckopacity_deinit_global(void);
void planckopacity_get_logT_limits(double *lower, double *upper);
int get_planck_opacitytable_value(struct planck_opacity_table *planck_table, planckopacitytable var, double rho, double t, double *result);
int get_planck_opacity_local(struct planck_opacity_table *planck_table, double rho, double t, struct planck_opacity_result *result);
int get_planck_opacity_global(double rho, double t, struct planck_opacity_result *result);

#define planckopacity_init planckopacity_init_global
#define planckopacity_deinit planckopacity_deinit_global
#define get_planck_opacity get_planck_opacity_global

#endif



#endif /* OPACITY_MESA */
