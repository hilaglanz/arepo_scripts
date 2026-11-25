/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/IDORT/IDORT_run.c
 * \date        08/2023
 * \author      Jing-Ze Ma
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"
#include "IDORT_proto.h"


#if defined(EOS_SUBTRACT_RAD) || defined(EOS_NO_RAD) || defined(IRT_STAR_PH_EOS_SWITCH)
void subtract_radiation(void)
{
    
    if(RestartFlag != RESTART_RESTART && RestartFlag != RESTART_SNAPSHOT)
      {
#if defined(EOS_SUBTRACT_RAD) || defined(IRT_STAR_PH_EOS_SWITCH)
    int i,idx;
#ifdef EOS_DEGENERATE
    struct eos_result res;
#else
#ifdef EOS_OPAL
    struct opal_eos_result res;
#endif
#endif
    
//    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
//      {
//        i = TimeBinsHydro.ActiveParticleList[idx];
//        if(i < 0)
//          continue;
    for(int i = 0; i < NumGas; i++)
      {
        set_pressure_of_cell_internal(P, SphP, i); /* calculate the pressure from Density and Utherm (and composition) */
      }
#ifdef EOS_SUBTRACT_RAD
    EosNoRad = 1;
#endif
#ifdef IRT_STAR_PH_EOS_SWITCH
    EosCoreRad = 1;
#endif
//    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
//      {
//        i = TimeBinsHydro.ActiveParticleList[idx];
//        if(i < 0)
//          continue;
    for(int i = 0; i < NumGas; i++)
      {
          
          SphP[i].Energy -= SphP[i].Utherm * P[i].Mass;
#ifdef EOS_DEGENERATE
          if(eos_calc_tgiven(SphP[i].Density, SphP[i].Composition, SphP[i].EOSTemperature, &res) < 0)
          {
          mpi_printf("Warning: Error in getting E from T...");
          }
          SphP[i].Utherm = res.e.v;
#else
#ifdef EOS_OPAL
          if(opaleos_all_tgiven(SphP[i].Composition[0], SphP[i].Density, SphP[i].EOSTemperature, &res) < 0)
          {
          mpi_printf("Warning: Error in getting E from T...");
          }
          SphP[i].Utherm = res.e;
#endif
#endif
          SphP[i].Energy += SphP[i].Utherm * P[i].Mass;
          set_pressure_of_cell_internal(P, SphP, i); /* calculate the pressure from Density and Utherm (and composition) */
      }
    mpi_printf("Finish subtracting radiation....\n");
#endif
#if defined(EOS_NO_RAD)
    EosNoRad = 1;
#endif
      }
else
    {
#if defined(EOS_SUBTRACT_RAD) || defined(EOS_NO_RAD)
    EosNoRad = 1;
#endif
#ifdef IRT_STAR_PH_EOS_SWITCH
    EosCoreRad = 1;
#endif
    }
}
#endif


#if defined(IDORT) || defined(OUTPUT_OPACITY)
void update_primitive_variables_IRT(int ic)
{
    
    int idx, i;
    double rad_const, clight;
#ifdef IDORT_TEST
    rad_const = All.IdortRadConst;
#else
    rad_const = RAD_CONST;
#endif
#ifdef IDORT_TEST
    clight = All.IdortClight;
#else
    clight = CLIGHT_REAL;
#endif
    
#ifdef IDORT_TEST
#pragma omp parallel for private(idx, i)
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        SphP[i].EOSTemperature = SphP[i].Utherm * GAMMA_MINUS1;
#ifdef IDORT_TEST
        if(IDORT_TEST == 8)
            SphP[i].EOSTemperature *= GAMMA;
        if(IDORT_TEST == 9)
        {
            double mu = 2.33;
            SphP[i].EOSTemperature /= BOLTZMANN / PROTONMASS / mu;
        }
#endif
    }
#endif
    
#ifdef OPACITY_MESA
      /* Snapshot/checkpoint restarts already have radiation. Only fresh ICs
       * may bootstrap cells whose radiation has not yet been initialized. */
      update_opacity(ic == 1 && RestartFlag == RESTART_IC);
#ifdef OPACITY_RP_OPTHIN
    if(ic != 1)
    {
        for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
          {
            i = TimeBinsHydro.ActiveParticleList[idx];
            if(i < 0)
              continue;
//            double Tedd_r = SphP[i].P_rad_rr/SphP[i].E_rad;
//            double x = 10.0*fabs(Tedd_r - 1.0/3.0); //1000.0*
            double x = 1.0/pow(10.0, (log10(SphP[i].tau)/1.0));
            double f = 0.0;
            if(x > 10.0)
            {
                f = 1.0;
            }
            else if(x < 0.1)
            {
                f = x;
            }
            else
            {
                f = (1.0-exp(-x))/(1.0+exp(-x));
            }
            SphP[i].KappaRoss = pow(10.0, log10(SphP[i].KappaRoss) * (1.0 - f) + log10(SphP[i].KappaPlanck) * f);
          }
    }
#endif
#else

#ifdef OPACITY_INDEX_TEST
#pragma omp parallel for private(idx, i)
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        SphP[i].KappaRoss = SphP[i].Density * pow(SphP[i].EOSTemperature, OPACITY_INDEX_TEST);
        SphP[i].KappaPlanck = SphP[i].KappaRoss;
    }
#endif
    
#ifdef IDORT_TEST
    if(IDORT_TEST == 9)
    {
#pragma omp parallel for private(idx, i)
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;

        SphP[i].KappaRoss = 0.0316 * pow(SphP[i].EOSTemperature/10.0, 2);
        SphP[i].KappaPlanck = 0.1 * pow(SphP[i].EOSTemperature/10.0, 2);
      }
    }
#endif
    
#endif

    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
          
      double rho = SphP[i].Density * All.UnitDensity_in_cgs;
      double kappa = SphP[i].KappaRoss + SphP[i].KappaScatter;
      double dx = get_cell_radius(i);

      SphP[i].tau = kappa * rho * dx;
        
#ifdef IDORT
        for(int j = 0; j < IDORT_NUM_ANGLES; j++)
        {
            SphP[i].Iperrho[j] = SphP[i].IV[j] / P[i].Mass;
#ifdef IDORT_GLOBAL_TIMESTEPS
            if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
            {
                SphP[i].Intensity[j] = SphP[i].IV[j] / SphP[i].Volume;
            }
        }
#endif
          
    }
}
#endif

#ifdef IDORT

double ray_dir[IDORT_NUM_ANGLES][3];
double ray_weight[IDORT_NUM_ANGLES];
double rad_const, clight;

static struct Dir * Dir;
static struct Fval * Fval;
static struct Cval * Cval;

#ifdef ONEDIMS_SPHERICAL
    double coszeta_v[IDORT_NUM_ANGLES], coszeta_f[IDORT_NUM_ANGLES+1], len_zeta[IDORT_NUM_ANGLES];
#endif

static void iter_IDORT(tessellation *T, int ic, int stage, double *dI_top, double *dI_bot, double *dIdivI, double *dTdivT);
static void compute_and_check_residuals_IDORT(int ic);
static void calculate_Ic(tessellation *T, int ic, int update_all, int correct_iso);
static int calculate_I(tessellation *T, int ic, int switch_iter, int stage);
static void get_NR_temperature(int i, double k1, double b, double T_lastiter);
static void initialize_directions(void);
static void initialize_IDORT(tessellation *T, double ray_dir[IDORT_NUM_ANGLES][3], int ic, int stage);
static void get_io_IDORT(int ic, int update_rad);
static void couple_IDORT(tessellation *T, int ic);
static void boundaries_IDORT(tessellation *T);



void run_IDORT(tessellation *T, int ic)
{
    
    TIMER_START(CPU_IDORT);
    
    Fval = (struct Fval *)mymalloc_movable(&Fval, "Fval", Mesh.Nvf * sizeof(struct Fval));
    Cval = (struct Cval *)mymalloc_movable(&Cval, "Cval", All.MaxPartSph * sizeof(struct Cval));
    
    TIMER_START(CPU_IDORT_ITER);
    
    int switch_iter = 1;
#ifndef IDORT_RELAX_INI
    if(ic == 1)
        switch_iter = 0;
#endif
    
    int iter_error = 0;
#ifdef IDORT_ONLY_ROBUST_CONVERGE
    iter_error = calculate_I(T, ic, switch_iter, 1);
#else
    iter_error = calculate_I(T, ic, switch_iter, 0);
    if(iter_error < 0)
    {
        iter_error = calculate_I(T, ic, switch_iter, 1);
        if(iter_error < 0)
        {
            iter_error = calculate_I(T, ic, switch_iter, 2);
            if(iter_error < 0)
              iter_error = calculate_I(T, ic, switch_iter, 1);
        }
    }
#endif
    mpi_printf("IDORT IMPLICIT: Finish iteration....\n");
    TIMER_STOP(CPU_IDORT_ITER);
    
    int switch_flux = 1;
#if defined(FORCE_EQUAL_TIMESTEPS) || defined(IDORT_GLOBAL_TIMESTEPS)
    switch_flux = 0;
#endif

    /* apply fluxes to inactive cells with neighboring active cells */
    if(switch_flux == 1)
    {
        exchange_primitive_variables_IDORT(2, Cval);
        initialize_finite_volume_coeff_vc(T, ray_dir, ic, Cval, Fval);
        boundaries_IDORT(T);
    }
        
        
    int update_rad = 1;
    get_io_IDORT(ic, update_rad);
        
    
#ifdef IDORT_COUPLE
    if(ic != 1)
    {
        couple_IDORT(T, ic);
    }
#endif
    
        
//    boundaries_IDORT(T);
    
    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
        
        for(int j = 0; j < IDORT_NUM_ANGLES; j++)
        {
            SphP[i].IV[j] = SphP[i].Intensity[j] * SphP[i].Volume;
        }
    }
    
    mpi_printf("IDORT: Finish....\n");
        
    myfree(Cval);
    myfree(Fval);
        
TIMER_STOP(CPU_IDORT);
}



int calculate_I(tessellation *T, int ic, int switch_iter, int stage)
{
    double dI_glob = 1e3, dIdivI_glob = 1e3, dTdivT_glob = 1e3, Iw;
    double b[3], n[3];
    int i_iter = 0;
    int i_nt = 0, i_nv = 0;
    int i, j, idx, k1, k2, l;
    double dt;
    double cdt_inv = 0.0;
    double Vdt_inv = 0.0;
    
    initialize_IDORT(T, ray_dir, ic, stage);

//    TIMER_START(CPU_IDORT_ITER);

    int max_iter_num = All.IdortMaxIterNum;
    int min_iter_num = All.IdortMinIterNum;
    if(ic == 1)
    {
#ifdef IDORT_MAX_ITER_NUM_INI
        max_iter_num = IDORT_MAX_ITER_NUM_INI;
#endif
#ifdef IDORT_MIN_ITER_NUM_INI
        min_iter_num = IDORT_MIN_ITER_NUM_INI;
#endif
    }

    double dI_lastiter = 1e5, dIdivI_lastiter = 1e5, dTdivT_lastiter = 1e5;
    if(switch_iter == 1)
    {
        mpi_printf("IDORT IMPLICIT: Iterating....\n");
        while((dI_glob > All.IdortTolerance || i_iter < min_iter_num) && i_iter < max_iter_num)
        {
            int update_all = 1;
            calculate_Ic(T, ic, update_all, 0);
            
            double dI_top, dI_bot, dI_top_glob, dI_bot_glob, dIdivI, dTdivT;
            iter_IDORT(T, ic, stage, &dI_top, &dI_bot, &dIdivI, &dTdivT);
            
            MPI_Allreduce(&dI_top, &dI_top_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&dI_bot, &dI_bot_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            dI_glob = dI_top_glob / dI_bot_glob;
            
            MPI_Allreduce(&dIdivI, &dIdivI_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&dTdivT, &dTdivT_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
            if(ic == 0)
            {
#ifndef IDORT_NEGLECT_TIME
                compute_and_check_residuals_IDORT(ic);
#endif
            }
            

            if(((dI_glob > dI_lastiter && ic != 1) || isnan(dI_glob)))// && (dI_glob < All.IdortTolerance || stage == 1))// || i_iter > min_iter_num))
            {
                mpi_printf("IDORT IMPLICIT: diverge at %d stage. switch to iteration at next stage.\n", stage);
                return -1;
            }
            dI_lastiter = dI_glob;
            mpi_printf("IDORT IMPLICIT: iter=%d, dI = %lf, logdI = %lf\n", i_iter, dI_glob, log10(dI_glob));
            
            boundaries_IDORT(T);

            for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
              {
                i = TimeBinsHydro.ActiveParticleList[idx];
                if(i < 0)
                  continue;

                  for(int j = 0; j < IDORT_NUM_ANGLES; j++)
                  {
                      Cval[i].I_lastiter[j] = SphP[i].Intensity[j];
                  }
              }
                
            i_iter += 1;
                
                if(i_iter == max_iter_num)
                {
                        mpi_printf("IDORT IMPLICIT: reached maximum iteration limit %d stage. switch to iteration at next stage.\n", stage);
                        return -1;
                }
        }
    }
//    TIMER_STOP(CPU_IDORT_ITER);
    return 0;
}



void calculate_Ic(tessellation *T, int ic, int update_all, int correct_iso)
{
    int i, j, k, idx, s, l;


    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
    
          double bv[3], b2;
          for(l = 0; l < 3; l++)
          {
#ifdef IDORT_NEGLECT_VEL
              bv[l] = 0.0;
#else
              if(ic != 1)
                  bv[l] = P[i].Vel[l] * All.UnitVelocity_in_cm_per_s / clight;
#endif
              b2 += bv[l] * bv[l];
          }

          for(j = 0; j < IDORT_NUM_ANGLES; j++)
          {
              Cval[i].Ic[j] = 0.0;
#ifdef IRT_STAR_PH_EOS_SWITCH
              Cval[i].Ic_corr[j] = 0.0;
#endif
#ifdef ONEDIMS_SPHERICAL
              double dnx;
              dnx = - 2.0 / (IDORT_NUM_ANGLES-1);
              double mu = ray_dir[j][0];
              double mu_l = coszeta_f[j]; //mu - 0.5*dnx;
              double mu_r = coszeta_f[j+1]; //mu + 0.5*dnx;
              double I_l;
              if(j != 0)
                  I_l = SphP[i].Intensity[j-1];
              else
                  I_l = SphP[i].Intensity[j];
              double I_r = SphP[i].Intensity[j];
              double dIdmu = ((1-mu_l*mu_l)*I_l - (1-mu_r*mu_r)*I_r) / (mu_l - mu_r);

              double gadj= -(1-mu_r*mu_r) / (mu_l - mu_r) / SphP[i].Center[0] * SphP[i].Volume;
              double gc   = (1-mu_l*mu_l) / (mu_l - mu_r) / SphP[i].Center[0] * SphP[i].Volume;
              if(j == IDORT_NUM_ANGLES - 1)
              {
                  gc += gadj;
                  gadj = 0.0;
              }
              else
              {
                  Cval[i].Ic[j] = - gadj * SphP[i].Intensity[j+1];
              }
              
#ifndef IDORT_NEGLECT_TIME
        if(ic == 0)
        {
            double dt;
#ifdef IDORT_GLOBAL_TIMESTEPS
            dt = (((integertime)1) << All.HighestActiveTimeBin) * All.Timebase_interval;
#else
            dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
#endif
            Cval[i].Ic[j] *= dt;
        }
#endif
#endif
              
          }
      }
        
    exchange_primitive_variables_IDORT(correct_iso, Cval);
    initialize_finite_volume_coeff_vb(T, ray_dir, ic, Cval, Fval, update_all, correct_iso);
}



void iter_IDORT(tessellation *T, int ic, int stage, double *dI_top, double *dI_bot, double *dIdivI, double *dTdivT)
{
    double JB, denom, denom_inv, f1, f2, k1, b, W_ni, I0_ni, Ib_ni, J0_i, Utherm, I0_b1[IDORT_NUM_ANGLES], I0_b2[IDORT_NUM_ANGLES], I0_k[IDORT_NUM_ANGLES];
    double bv[3], b2;
    double Ic_c, gcp, gcm;
    double T_lastiter;
    *dI_bot = 0.0;
    *dI_top = 0.0;
    *dIdivI = 0.0;
    *dTdivT = 0.0;
    double dT_bot = 0.0;
    double dT_top = 0.0;
    double dI = 0.0;
    double dT = 0.0;
    double dt = 0.0;
    int i, j, k, idx, s, l;
   
    for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
          
#ifdef IRT_STAR_CORE
          if(get_distance_from_core(i) <= All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
              continue;
#endif
        double cdt_inv = 0.0;
        double Vdt_inv = 0.0;
          
        Vdt_inv = 1.0 / SphP[i].Volume;
        if(ic == 0)
        {
#ifdef IDORT_GLOBAL_TIMESTEPS
            dt = (((integertime)1) << All.HighestActiveTimeBin) * All.Timebase_interval;
#else
            dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
#endif
//                dt *= 0.5;
#ifndef IDORT_NEGLECT_TIME
                cdt_inv = 1.0 / (clight * dt * All.UnitTime_in_s);
                Vdt_inv = 1.0 / (SphP[i].Volume * dt);
#endif
        }
          T_lastiter = Cval[i].EOSTemperature;
          JB = rad_const*pow(Cval[i].EOSTemperature,4)/4.0/M_PI;
          double sum_top_i = 0.0;
          double sum_bot_i = 0.0;
          double sum_Ic = 0.0;
          double A = 0.0;
          double B = 0.0;
          for(j = 0; j < IDORT_NUM_ANGLES; j++)
          {
              double gm = 0.0;
              if(stage == 1)
                  gm = Cval[i].gm[j];
              Ic_c = Cval[i].Ic[j] * Vdt_inv - gm * Vdt_inv * SphP[i].Intensity[j] + cdt_inv * SphP[i].IV[j] / SphP[i].Volume;
              I0_b1[j] = Cval[i].f[j] * Ic_c;
              sum_Ic += Cval[i].fsi[j] * Ic_c;
          }
          k1 = Cval[i].k;
          b = Cval[i].fbs * sum_Ic + SphP[i].Utherm;

          if(ic == 0)
          {
              int IT_COEVOL = 0;
#if defined(IDORT_IT_COEVOL)// && !defined(IDORT_NEGLECT_TIME)
#ifdef IDORT_NEGLECT_TIME
              if(TimeBinSynchronized[P[i].TimeBinHydro] && P[i].Type == PTYPE_GAS && P[i].ID != 0 && P[i].Mass > 0)
                  IT_COEVOL = 1;
#else
              IT_COEVOL = 1;
#endif
#ifdef IRT_NO_IT_COEVOL_OPTHIN
              if(SphP[i].tau < 1.0)
                  IT_COEVOL = 0;
#endif
#ifdef IRT_NO_IT_COEVOL_OPTHICK
              if(SphP[i].tau > 1.0)
                  IT_COEVOL = 0;
#endif
#endif
#if defined(IRT_STAR_CORE)
              if(get_distance_from_core(i) < All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
              {
                  IT_COEVOL = 0;
              }
#endif
        if(IT_COEVOL == 1)
        {
            get_NR_temperature(i, k1, b, T_lastiter);
        }
          
          dT_bot += T_lastiter;
          dT_top += fabs(Cval[i].EOSTemperature - T_lastiter);
          if(T_lastiter > 0)
              *dTdivT += fabs(Cval[i].EOSTemperature - T_lastiter)/T_lastiter;
          
         } /* ic != 1 */

         JB = rad_const*pow(Cval[i].EOSTemperature,4)/4.0/M_PI;
          
          for(j = 0; j < IDORT_NUM_ANGLES; j++)
          {
              Ib_ni = I0_b1[j] + Cval[i].fis[j] * sum_Ic + Cval[i].fib[j] * JB;
              if(Ib_ni < 0.0)
                  Ib_ni = 0.0;
             *dI_bot += fabs(SphP[i].Intensity[j]);
             *dI_top += fabs(Ib_ni - SphP[i].Intensity[j]);
             if(SphP[i].Intensity[j] > 0)
                 *dIdivI += fabs(Ib_ni - SphP[i].Intensity[j])/SphP[i].Intensity[j];
             SphP[i].Intensity[j] = Ib_ni;
         }
          
         
      }

}
    
    

void compute_and_check_residuals_IDORT(int ic)
{
    double dEV_top = 0.0;
    double dEV_bot = 0.0;
    double dFpercV_top = 0.0;
    double dFpercV_bot = 0.0;
    double maxresEV = 0.0;
    int maxresEV_index = -1;
    double maxresFpercV = 0.0;
    int maxresFpercV_index = -1;
    
    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
          
#if defined(IRT_STAR_CORE)
          if(get_distance_from_core(i) <= All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
              continue;
#endif
          
        double Vdt_inv = 1.0 / SphP[i].Volume;
        double cdt_inv = 0.0;
        if(ic == 0)
        {
#ifdef IDORT_GLOBAL_TIMESTEPS
            double dt = (((integertime)1) << All.HighestActiveTimeBin) * All.Timebase_interval;
#else
            double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
#endif
//                dt *= 0.5;
#ifndef IDORT_NEGLECT_TIME
                cdt_inv = 1.0 / (clight * dt * All.UnitTime_in_s);
                Vdt_inv = 1.0 / (SphP[i].Volume * dt);
#endif
        }
        double b2 = 0.0, b[3];
        for(int q = 0; q < 3; q++)
        {
#ifdef IDORT_NEGLECT_VEL
            b[q] = 0.0;
#else
            if(ic != 1)
            {
                b[q] = P[i].Vel[q] * All.UnitVelocity_in_cm_per_s / clight;
            }
            else
            {
                b[q] = 0.0;
            }
#endif
            b2 += b[q] * b[q];
        }
        double gamma = 1.0 / sqrt(1.0 - b2);
        double partition = 0.0;
        double JB = rad_const*pow(Cval[i].EOSTemperature,4)/4.0/M_PI;
        double Jp = 0.0;
        for(int n = 0; n < IDORT_NUM_ANGLES; n++)
        {
            
            double Gamma = 1.0;
            for(int q = 0; q < 3; q++)
            {
                Gamma -= b[q] * ray_dir[n][q];
            }
            Gamma *= gamma;
            double invG = 1.0 / Gamma;
            double wp = ray_weight[n] * invG * invG; /* missing 1/partition */
            partition += wp;
            
            double Ip = Gamma*Gamma*Gamma*Gamma*SphP[i].Intensity[n];
            Jp += Ip*wp; /* missing 1/partition */
        }
        Jp /= partition;
        double resEFVr[4], EFVr[4];
        for(int q = 0; q < 4; q++)
        {
            resEFVr[q] = 0.0;
            EFVr[q] = 0.0;
        }
        for(int n = 0; n < IDORT_NUM_ANGLES; n++)
        {
            double Gamma = 1.0;
            for(int q = 0; q < 3; q++)
            {
                Gamma -= b[q] * ray_dir[n][q];
            }
            Gamma *= gamma;
            double invG = 1.0 / Gamma;
            
            double nablacnI = clight * ((Cval[i].gp[n] + Cval[i].gm[n]) * SphP[i].Intensity[n] - Cval[i].Ic[n]) * Vdt_inv;
            double dIdt = clight * cdt_inv * (SphP[i].Intensity[n] - SphP[i].IV[n] / SphP[i].Volume);
            double Ip = Gamma*Gamma*Gamma*Gamma*SphP[i].Intensity[n];
            double SI = invG*invG*invG * SphP[i].Density * (SphP[i].KappaPlanck * JB + (SphP[i].KappaRoss + SphP[i].KappaScatter - SphP[i].KappaPlanck) * Jp - (SphP[i].KappaRoss + SphP[i].KappaScatter) * Ip);
            double resIV = (dIdt + nablacnI - clight * SI) / Vdt_inv;
            double IV = SphP[i].IV[n];
            resIV /= SphP[i].Volume;
            IV /= SphP[i].Volume;
            for(int q = 0; q < 4; q++)
            {
                if(q == 0)
                {
                    resEFVr[q] += 4.0 * M_PI * resIV * ray_weight[n];
                    EFVr[q] += 4.0 * M_PI * IV * ray_weight[n];
                }
                else
                {
                    resEFVr[q] += 4.0 * M_PI * resIV * ray_dir[n][q-1] * ray_weight[n];
                    EFVr[q] += 4.0 * M_PI * IV * ray_dir[n][q-1] * ray_weight[n];
                }
            }
        }
          dEV_bot += fabs(EFVr[0]);
          dEV_top += fabs(resEFVr[0]);
          double FpercVsq_old = 0.0;
          double resFpercVsq = 0.0;
          for(int q = 0; q < 3; q++)
          {
              FpercVsq_old += EFVr[q+1]*EFVr[q+1];
              resFpercVsq += resEFVr[q+1]*resEFVr[q+1];
          }
          dFpercV_bot += sqrt(FpercVsq_old);
          dFpercV_top += sqrt(resFpercVsq);
          
          double resEVp = 0.0;
          if(EFVr[0] > 0.0)
              resEVp = resEFVr[0]/EFVr[0];
          if(resEVp > maxresEV)
          {
              maxresEV = resEVp;
              maxresEV_index = i;
          }
          double resFpercVp = 0.0;
          if(sqrt(FpercVsq_old) > 0.0)
              resFpercVp = sqrt(resFpercVsq)/sqrt(FpercVsq_old);
          if(resFpercVp > maxresFpercV)
          {
              maxresFpercV = resFpercVp;
              maxresFpercV_index = i;
          }
      }
        
        double resEVall, EVall;
        double resFpercVall, FpercVall;
        MPI_Allreduce(&dEV_top, &resEVall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&dEV_bot, &EVall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&dFpercV_top, &resFpercVall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&dFpercV_bot, &FpercVall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double resEV = 0.0;
        if(EVall > 0.0)
            resEV = resEVall / EVall;
        double resFpercV = 0.0;
        if(FpercVall > 0.0)
            resFpercV = resFpercVall / FpercVall;
            mpi_printf("IDORT IMPLICIT: Energy residual: %g, Flux residual: %g\n", resEV, resFpercV);

        struct
        {
          double val;
          int rank;
        } inE, outE, inF, outF;

        inE.val  = maxresEV;
        inE.rank = ThisTask;

        MPI_Allreduce(&inE, &outE, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        // Broadcast the maximum index to all tasks
        int global_max_index = -1;
        if (ThisTask == outE.rank) {
            global_max_index = maxresEV_index;
        }
        MPI_Bcast(&global_max_index, 1, MPI_INT, outE.rank, MPI_COMM_WORLD);

        // Output only from the task with the maximum value
        if (ThisTask == outE.rank && global_max_index != -1) {

            int p = global_max_index;//maxresEV_index;
              double E_rad_old = 0.0;
              double E_rad_new = 0.0;
              for(int n = 0; n < IDORT_NUM_ANGLES; n++)
              {
                  E_rad_old += 4.0 * M_PI * SphP[p].IV[n]/SphP[p].Volume * ray_weight[n];
                  E_rad_new += 4.0 * M_PI * SphP[p].Intensity[n] * ray_weight[n];
              }
              printf(
                "IDORT IMPLICIT: Particle %d on task %d has the highest energy residual res=%g, rho=%g, T=%g, Tradold=%g, Trad=%g.\n",
                p, ThisTask, maxresEV, SphP[p].Density, SphP[p].EOSTemperature, pow(E_rad_old/rad_const,0.25), pow(E_rad_new/rad_const,0.25));
#ifdef IRT_STAR_CORE
              double r = get_distance_from_core(p);
              printf("Radius=%g.\n",r/6.957e10);
#endif
          }
        
        inF.val  = maxresFpercV;
        inF.rank = ThisTask;

        MPI_Allreduce(&inF, &outF, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        if(ThisTask == outF.rank && maxresFpercV_index != -1)// && ic != 1)
          {
            int p = maxresEV_index;
              double E_rad_old = 0.0;
              double E_rad_new = 0.0;
              for(int n = 0; n < IDORT_NUM_ANGLES; n++)
              {
                  E_rad_old += 4.0 * M_PI * SphP[p].IV[n]/SphP[p].Volume * ray_weight[n];
                  E_rad_new += 4.0 * M_PI * SphP[p].Intensity[n] * ray_weight[n];
              }
              printf(
                "IDORT IMPLICIT: Particle %d on task %d has the highest flux residual res=%g, rho=%g, T=%g, Tradold=%g, Trad=%g.\n",
                p, ThisTask, maxresFpercV, SphP[p].Density, SphP[p].EOSTemperature, pow(E_rad_old/rad_const,0.25), pow(E_rad_new/rad_const,0.25));
#ifdef IRT_STAR_CORE
              double r = get_distance_from_core(p);
              printf("Radius=%g.\n",r/6.957e10);
#endif
          }
}
    
    

void get_NR_temperature(int i, double k1, double b, double T_lastiter)
{
    double Utherm = SphP[i].Utherm;
#ifdef EOS_DEGENERATE
      struct eos_result res;

      if(eos_calc_et_relation_given(SphP[i].Density, SphP[i].Composition, Utherm, &Cval[i].EOSTemperature, &res, k1, b) ==
         -2)
        {
          Cval[i].EOSTemperature = -1.;
        }
    if(eos_calc_et_relation_given(SphP[i].Density, SphP[i].Composition, Utherm, &Cval[i].EOSTemperature,
                       &res, k1, b) != 0)
      {
        if(res.temp < 1000.)
        {
          Cval[i].EOSTemperature = 1000;
        }
        else if(res.temp > 1e10)
        {
          Cval[i].EOSTemperature = 1e10;
        }
#ifdef IRT_IT_COEVOL_DEBUG
        mpi_printf("Warning: Error in I and T coevol... last iter %lf, new %lf\n", T_lastiter, Cval[i].EOSTemperature);
#endif
      }
//          SphP[i].Utherm = res.e;
#else
#ifdef EOS_OPAL
      struct opal_eos_result res;
      if(opaleos_et_relation_given(SphP[i].Composition[0], SphP[i].Density, Utherm, &Cval[i].EOSTemperature, EOS_ENERGY | EOS_CV | EOS_RADIATION, &res, k1, b) < 0)
      {
#ifdef IRT_IT_COEVOL_DEBUG
      mpi_printf("Warning: Error in I and T coevol... last iter %lf, new %lf\n", T_lastiter, Cval[i].EOSTemperature);
#endif
      }
#else
     double t, de, told, e, f;
     t = Cval[i].EOSTemperature;
     f = 1 / GAMMA_MINUS1;
#ifdef IDORT_TEST
     if(IDORT_TEST == 8)
         f /= GAMMA;
    if(IDORT_TEST == 9)
    {
        double mu = 2.33;
        f *= BOLTZMANN / PROTONMASS / mu;
    }
#endif
//         e = SphP[i].Utherm;
     for(int s = 0; s < 1000; s++)
       {
         e = t * f;
         de   = (e - k1*rad_const*t*t*t*t - b);
         told = t;
         t -= de / (f - 4.0*k1*rad_const*t*t*t);
         /* converged to desired accuracy? */
         if(fabs(de) <= 1e-8 * e)
           break;
        }
     Cval[i].EOSTemperature = t;
//         SphP[i].Utherm = e;
#endif
#endif
}
   

        
void get_io_IDORT(int ic, int update_rad)
{
    double dir_r[3], x[3], r, pos_rgcore[3];

    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
        
#ifdef ID_RGCORE
            get_pos_rgcore(pos_rgcore);
#else
            pos_rgcore[0] = boxHalf_X;
            pos_rgcore[1] = boxHalf_Y;
            pos_rgcore[2] = boxHalf_Z;
#endif
            x[0] = P[i].Pos[0] - pos_rgcore[0];
            r = x[0]*x[0];
            if(NUMDIMS > 1)
            {
                x[1] = P[i].Pos[1] - pos_rgcore[1];
                r += x[1]*x[1];
                if(NUMDIMS > 2)
                {
                    x[2] = P[i].Pos[2] - pos_rgcore[2];
                    r += x[2]*x[2];
                }
            }
            r = sqrt(r);          
          
          if(update_rad == 1)
          {
              SphP[i].E_rad = 0.0;
          }
          SphP[i].F_rad_r = 0.0;
          SphP[i].F_rad_r_com = 0.0;
          SphP[i].P_rad_rr = 0.0;
          for(int k1 = 0; k1< 3; k1++)
          {
              if(update_rad == 1)
              {
                  SphP[i].F_rad[k1] = 0.0;
              }
              for(int k2 = 0; k2< 3; k2++)
              {
                  SphP[i].P_rad[k1][k2] = 0.0;
              }
          }
          
          for(int l = 0; l< 3; l++)
          {
              if(r > 0.0)
              {
                  dir_r[l] = x[l] / r;
              }
              else
              {
                  dir_r[l] = 0.0;
              }
          }
          for(int j = 0; j < IDORT_NUM_ANGLES; j++)
          {
              double Iw = SphP[i].Intensity[j] * ray_weight[j];
              if(update_rad == 1)
              {
                  SphP[i].E_rad += 4.0 * M_PI * Iw;
              }
              
              for(int k1 = 0; k1< 3; k1++)
              {
                  if(update_rad == 1)
                  {
                      SphP[i].F_rad[k1] += clight * 4.0 * M_PI * Iw * ray_dir[j][k1];
                  }
                  
                  for(int k2 = 0; k2< 3; k2++)
                  {
                      SphP[i].P_rad[k1][k2] += 4.0 * M_PI * Iw * ray_dir[j][k1] * ray_dir[j][k2];
                  }
              }
          }
          for(int k1 = 0; k1 < 3; k1++)
          {
              SphP[i].F_rad_com[k1] = SphP[i].F_rad[k1];
              if(ic != 1)
              {
#if defined(IRT_STAR_CORE) && defined(IDORT_CORE_CONST_L)
                  if(r > All.IrtStarCoreRadius)
                  {
#endif
#ifndef IDORT_NEGLECT_VEL
                  SphP[i].F_rad_com[k1] -= P[i].Vel[k1] * SphP[i].E_rad;
                  for(int k2 = 0; k2< 3; k2++)
                  {
                      SphP[i].F_rad_com[k1] -= P[i].Vel[k2] * SphP[i].P_rad[k1][k2];
                  }
#endif
#if defined(IRT_STAR_CORE) && defined(IDORT_CORE_CONST_L)
                  }
#endif
              }
          }
          for(int k1 = 0; k1< 3; k1++)
          {
              SphP[i].F_rad_r += SphP[i].F_rad[k1] * dir_r[k1];
              SphP[i].F_rad_r_com += SphP[i].F_rad_com[k1] * dir_r[k1];
              for(int k2 = 0; k2< 3; k2++)
              {
                  SphP[i].P_rad_rr += SphP[i].P_rad[k1][k2] * dir_r[k1] * dir_r[k2];
              }
          }
          if(ic == 1)
          {
              SphP[i].E_rad_com = SphP[i].E_rad;
          }
      }
}
        
        
        
void couple_IDORT(tessellation *T, int ic)
{

    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
          
#ifdef IRT_STAR_CORE
              double r = get_distance_from_core(i);
#endif
          
#if defined(IRT_STAR_CORE) && !defined(EOS_SUBTRACT_RAD) && !defined(EOS_NO_RAD)
    if(r < All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
        continue;
#endif
          
        if(P[i].Type == PTYPE_GAS && P[i].ID != 0 && P[i].Mass > 0)
            {
                
                double dt;
#ifdef IDORT_GLOBAL_TIMESTEPS
                dt = (((integertime)1) << All.HighestActiveTimeBin) * All.Timebase_interval;
#else
                dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
#endif

                double Vdt_inv = 1.0 / SphP[i].Volume;
#ifndef IDORT_NEGLECT_TIME
                Vdt_inv = 1.0 / (SphP[i].Volume * dt);
#endif
                
                double dE_rad, dF_rad[3];
                double Sr, Srv;
                double nablacnIw;
                
                  dE_rad = 0.0;
                  
                  for(int k1 = 0; k1< 3; k1++)
                  {
                      dF_rad[k1] = 0.0;
                  }
                  
                  for(int j = 0; j < IDORT_NUM_ANGLES; j++)
                  {
                      double I_in_flux = SphP[i].Intensity[j];//Cval[i].I_lastiter[j];
                      nablacnIw = clight * ((Cval[i].gp[j] + Cval[i].gm[j]) * I_in_flux - Cval[i].Ic[j]) * Vdt_inv * ray_weight[j];

                      dE_rad += 4.0 * M_PI * nablacnIw;
                      
                      for(int k1 = 0; k1< 3; k1++)
                      {
#ifdef IRT_STAR_PH_EOS_SWITCH
                          double rho = SphP[i].Density;
                          double factor = 0.0;
                          if(rho > 10.0 * All.IrtStarPhDensity)
                          {
                              factor = 1.0;
                          }
                          else if(rho > 0.1 * All.IrtStarPhDensity && rho <= 10.0 * All.IrtStarPhDensity)
                          {
                              factor = (rho - 0.1 * All.IrtStarPhDensity) / (10.0 * All.IrtStarPhDensity - 0.1 * All.IrtStarPhDensity);
                          }
                          dF_rad[k1] += clight * 4.0 * M_PI * (clight * ((Cval[i].gp[j] + Cval[i].gm[j]) * I_in_flux * (1.0-factor) - Cval[i].Ic_corr[j]) * Vdt_inv * ray_weight[j]) * ray_dir[j][k1];
#else
                          dF_rad[k1] += clight * 4.0 * M_PI * nablacnIw * ray_dir[j][k1];
#endif
                      }
                  }
                
                double KE_old = 0.0;
                double KE_new = 0.0;
                
                dE_rad *= dt;
                
                int non_conservative = 0;
#ifdef IDORT_COUPLE_NON_CONSERVATIVE_FORM
                non_conservative = 1;
#endif
                
            double EB = rad_const*pow(Cval[i].EOSTemperature,4);
            double b[3];
            if(non_conservative == 1)
            {
                double bF = 0.0;
                double bF_com = 0.0;
                double b2 = 0.0;
                for(int k1 = 0; k1< 3; k1++)
                {
#ifdef IDORT_NEGLECT_VEL
                    b[k1] = 0.0;
#else
                    b[k1] = P[i].Vel[k1] * All.UnitVelocity_in_cm_per_s / clight;
#endif
                    b2 += b[k1] * b[k1];
                    bF += b[k1] * SphP[i].F_rad[k1];
                    bF_com += b[k1] * SphP[i].F_rad_com[k1];
                }
                Sr = clight * dt * SphP[i].Density * All.UnitDensity_in_cgs * (SphP[i].KappaPlanck * ((EB - SphP[i].E_rad) + bF/clight) - (SphP[i].KappaScatter + SphP[i].KappaRoss - SphP[i].KappaPlanck) * bF_com/clight);//clight * dt * SphP[i].Density * All.UnitDensity_in_cgs * (SphP[i].KappaPlanck * ((EB - SphP[i].E_rad)*(1.0+0.5*b2) + bF/clight) - (SphP[i].KappaScatter + SphP[i].KappaRoss - SphP[i].KappaPlanck) * bF_com/clight);//
#ifdef IRT_STAR_CORE
              if(r < All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
              {
                  Sr = 0.0;
              }
#endif
            }
            else
            {
#if defined(IDORT_NEGLECT_VEL) && defined(IDORT_NEGLECT_TIME)
                Sr = dE_rad;
#else
                Sr = SphP[i].E_rad - Cval[i].E_rad_post_adv_gas + dE_rad;
#endif
#ifdef IRT_STAR_CORE
              if(r < All.IrtStarCoreRadius)
              {
                  if(r > All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
                    {
                        Sr *= (r/All.IrtStarCoreRadius - All.IrtStarCoreRadiusFactor)/(1.0-All.IrtStarCoreRadiusFactor);
                    }
                  else
                  {
                      Sr = 0.0;
                  }
              }
#endif
            }
                
                for(int k1 = 0; k1< 3; k1++)
                {
                    KE_old += SphP[i].Momentum[k1] * SphP[i].Momentum[k1] / P[i].Mass / 2.0;
                    dF_rad[k1] *= dt;
                    
                    if(non_conservative == 1)
                    {
                    Srv = dt * SphP[i].Density * All.UnitDensity_in_cgs * (SphP[i].KappaPlanck * (EB - SphP[i].E_rad) * b[k1] - (SphP[i].KappaScatter + SphP[i].KappaRoss) * SphP[i].F_rad_com[k1]/clight);
#ifdef IRT_STAR_CORE
                          if(r < All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
                          {
                              Srv = dF_rad[k1] / clight / clight;
#if !defined(EOS_SUBTRACT_RAD) && !defined(EOS_NO_RAD)
                              Srv = 0.0;
#endif
                          }
#endif
                    }
                    else
                    {
#if defined(IDORT_NEGLECT_VEL) && defined(IDORT_NEGLECT_TIME)
                    Srv = 0.0;
#ifndef IDORT_NEGLECT_TIME_NO_RAD_FORCE_OPTHIN
#ifdef IRT_STAR_PH_EOS_SWITCH
//                      double rho = SphP[i].Density;
//                      double factor = 0.0;
//                      if(rho > 10.0 * All.IrtStarPhDensity)
//                      {
//                          factor = 1.0;
//                      }
//                      else if(rho > 0.1 * All.IrtStarPhDensity && rho <= 10.0 * All.IrtStarPhDensity)
//                      {
//                          factor = (rho - 0.1 * All.IrtStarPhDensity) / (10.0 * All.IrtStarPhDensity - 0.1 * All.IrtStarPhDensity);
//                      }
//                    Srv = (1.0 - factor) * dF_rad[k1] / clight / clight;
                        Srv = dF_rad[k1] / clight / clight;
#ifdef IDORT_NEGLECT_TIME_RAD_FORCE_OPTHIN_PLANCK
                    if(SphP[i].tau < 1.0)
                        Srv *= (SphP[i].KappaPlanck+SphP[i].KappaScatter)/SphP[i].KappaRoss;
#endif
#else
                    if(SphP[i].tau < 1.0)
                        Srv = dF_rad[k1] / clight / clight;
#endif
#endif
#else
                    Srv = (SphP[i].F_rad[k1] - Cval[i].F_rad_post_adv_gas[k1] + dF_rad[k1]) / clight / clight;
#endif
                    }

                    SphP[i].Momentum[k1] += - Srv * SphP[i].Volume;
                    KE_new += SphP[i].Momentum[k1] * SphP[i].Momentum[k1] / P[i].Mass / 2.0;
                }
                SphP[i].Energy += - Sr * SphP[i].Volume;
            }
        }
}



void initialize_IDORT(tessellation *T, double ray_dir[IDORT_NUM_ANGLES][3], int ic, int stage)
{

//  TIMER_START(CPU_IDORT_INIT);
#ifdef IDORT_ONLY_ROBUST_CONVERGE
if(stage == 0 || stage == 1)
#else
if(stage == 0)
#endif
{
  initialize_directions();
    
#ifdef IDORT_TEST
    rad_const = All.IdortRadConst;
#else
    rad_const = RAD_CONST;
#endif
#ifdef IDORT_TEST
    clight = All.IdortClight;
#else
    clight = CLIGHT_REAL;
#endif
    
    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
          
          double b2 = 0.0, b[3];
          for(int q = 0; q < 3; q++)
          {
#ifdef IDORT_NEGLECT_VEL
              b[q] = 0.0;
#else
              if(ic != 1)
              {
                  b[q] = P[i].Vel[q] * All.UnitVelocity_in_cm_per_s / clight;
              }
              else
              {
                  b[q] = 0.0;
              }
#endif
              b2 += b[q] * b[q];
          }
          double gamma = 1.0 / sqrt(1.0 - b2);
          for(int n = 0; n < IDORT_NUM_ANGLES; n++)
            {
                double Gamma = 1.0;
                for(int q = 0; q < 3; q++)
                {
                    Gamma -= b[q] * ray_dir[n][q];
                }
                Gamma *= gamma;
                
                Cval[i].gp[n] = 0.0;
                Cval[i].gm[n] = 0.0;
                
#ifdef ONEDIMS_SPHERICAL
                double dnx;
                dnx = - 2.0 / (IDORT_NUM_ANGLES-1);
                double mu = ray_dir[n][0];
                double mu_l = coszeta_f[n];
                double mu_r = coszeta_f[n+1];
                double I_l;
                if(n != 0)
                    I_l = SphP[i].Intensity[n-1];
                else
                    I_l = SphP[i].Intensity[n];
                double I_r = SphP[i].Intensity[n];
                double dIdmu = ((1-mu_l*mu_l)*I_l - (1-mu_r*mu_r)*I_r) / (mu_l - mu_r);

                double gadj= -(1-mu_r*mu_r) / (mu_l - mu_r) / SphP[i].Center[0] * SphP[i].Volume;
                double gc   = (1-mu_l*mu_l) / (mu_l - mu_r) / SphP[i].Center[0] * SphP[i].Volume;
                if(n == IDORT_NUM_ANGLES - 1)
                {
                    gc += gadj;
                    gadj = 0.0;
                }
                if(gc > 0)
                    Cval[i].gp[n]  = gc;
                else
                    Cval[i].gm[n]  = gc;
                
#ifndef IDORT_NEGLECT_TIME
              if(ic == 0)
              {
                  double dt;
#ifdef IDORT_GLOBAL_TIMESTEPS
                  dt = (((integertime)1) << All.HighestActiveTimeBin) * All.Timebase_interval;
#else
                  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
#endif
                  Cval[i].gp[n] *= dt;
                  Cval[i].gm[n] *= dt;
              }
#endif
#endif
        }
    }
    initialize_finite_volume_coeff_va(T, ray_dir, ic, Cval, Fval);
        
}
        

    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
      {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
          
          Cval[i].E_rad_post_adv_gas = 0.0;
          for(int k1 = 0; k1 < 3; k1++)
          {
              Cval[i].F_rad_post_adv_gas[k1] = 0.0;
          }
  
          double dt;
          double cdt_inv = 0.0;
          double Vdt_inv = 1.0 / SphP[i].Volume;

          if(ic == 0)
          {
#ifdef IDORT_GLOBAL_TIMESTEPS
              dt = (((integertime)1) << All.HighestActiveTimeBin) * All.Timebase_interval;
#else
              dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
#endif
#ifndef IDORT_NEGLECT_TIME
                  cdt_inv = 1.0 / (clight * dt * All.UnitTime_in_s);
                  Vdt_inv = 1.0 / (SphP[i].Volume * dt);
#endif
          }

          double b2 = 0.0, b[3];
          for(int q = 0; q < 3; q++)
          {
#ifdef IDORT_NEGLECT_VEL
              b[q] = 0.0;
#else
              if(ic != 1)
              {
                  b[q] = P[i].Vel[q] * All.UnitVelocity_in_cm_per_s / clight;
              }
              else
              {
                  b[q] = 0.0;
              }
#endif
              b2 += b[q] * b[q];
          }
          double gamma = 1.0 / sqrt(1.0 - b2);
          double partition = 0.0;
          double sum_wp_f_G = 0.0;
          for(int n = 0; n < IDORT_NUM_ANGLES; n++)
          {
              
              double Gamma = 1.0;
              for(int q = 0; q < 3; q++)
              {
                  Gamma -= b[q] * ray_dir[n][q];
              }
              Gamma *= gamma;
              double invG = 1.0 / Gamma;
              double wp = ray_weight[n] * invG * invG; /* missing 1/partition */
              partition += wp;
              double gp = Cval[i].gp[n] + Cval[i].gm[n];
              if(stage == 1)
                  gp = Cval[i].gp[n];
              Cval[i].f[n] = 1.0 / (cdt_inv + gp * Vdt_inv + SphP[i].Density * (SphP[i].KappaRoss + SphP[i].KappaScatter) * Gamma);
              Cval[i].fis[n] = SphP[i].Density * (SphP[i].KappaRoss + SphP[i].KappaScatter - SphP[i].KappaPlanck) * Cval[i].f[n] * invG * invG * invG; /* missing h */
              Cval[i].fib[n] = SphP[i].Density * SphP[i].KappaPlanck * Cval[i].f[n] * invG * invG * invG; /* missing h */
              Cval[i].fsi[n] = wp * Cval[i].f[n] * Gamma * Gamma * Gamma * Gamma; /* missing 1/partition */
              
              sum_wp_f_G += wp * Cval[i].f[n] * Gamma; /* missing 1/partition */
          }
          Cval[i].fbs = 4.0 * M_PI * clight * dt * SphP[i].KappaPlanck; /* missing h */
          double invP = 1.0 / partition;
          sum_wp_f_G *= invP;
          
          double h = 1.0 / (1.0 - (SphP[i].Density * (SphP[i].KappaRoss + SphP[i].KappaScatter - SphP[i].KappaPlanck) * sum_wp_f_G));
          Cval[i].k = clight * dt * SphP[i].KappaPlanck * ((h * SphP[i].Density * SphP[i].KappaPlanck * sum_wp_f_G) - 1);
          Cval[i].fbs *= h;
          
          double Ip, Im;
#ifdef IDORT_TEST
          if(IDORT_TEST == 9)
          {
              double T0 = 82.0;
              double Er0 = rad_const * pow(T0,4);
              double F0 = clight * rad_const * pow(T0,4);
//              Ip = F0 / clight * 6.0 / (4.0 * M_PI);
              double fE = 0.0;
              double fF = 0.0;
              for(int n = 0; n < IDORT_NUM_ANGLES; n++)
              {
                  if(ray_dir[n][1]>0 && ray_dir[n][2]<1e-3)
                  {
                      fE += ray_weight[n];
                      fF += ray_dir[n][1] * ray_weight[n];
                  }
              }
              fE *= 2.0;
              fF *= 2.0;
              Ip = Er0 / fE + F0 / clight / fF;
              Im = Er0 / fE - F0 / clight / fF;
              Ip /= 4.0 * M_PI;
              Im /= 4.0 * M_PI;
//              Ip = F0 / clight / fF;
//              Ip /= 4.0 * M_PI;
          }
#endif
          if(stage < 2)
              Cval[i].EOSTemperature = SphP[i].EOSTemperature;
          for(int n = 0; n < IDORT_NUM_ANGLES; n++)
          {
              Cval[i].fis[n] *= h;
              Cval[i].fib[n] *= h;
              Cval[i].fsi[n] *= invP;
              
//              SphP[i].Intensity[n] = SphP[i].Iperrho[n] * P[i].Mass / SphP[i].Volume;
              if(ic == 1)
              {
#ifdef RAD_INI_CON
                  SphP[i].Intensity[n] = SphP[i].E_rad/4.0/M_PI;
#else
                  SphP[i].Intensity[n] = rad_const*pow(SphP[i].EOSTemperature,4)/4.0/M_PI;
#endif
                  SphP[i].IV[n] = SphP[i].Intensity[n] * SphP[i].Volume;
              }
              if(stage < 2)
              {
#ifdef IDORT_GLOBAL_TIMESTEPS
                  SphP[i].Intensity[n] = SphP[i].IV[n] / SphP[i].Volume;
#endif
                  Cval[i].I_lastiter[n] = SphP[i].Intensity[n];
              }
              if(ic != 1)
              {
                  double Iw = SphP[i].IV[n] / SphP[i].Volume * ray_weight[n];
                  Cval[i].E_rad_post_adv_gas += 4.0 * M_PI * Iw;
                  for(int k1 = 0; k1 < 3; k1++)
                  {
                      Cval[i].F_rad_post_adv_gas[k1] += clight * 4.0 * M_PI * Iw * ray_dir[n][k1];
                  }
              }
          }
      }
    boundaries_IDORT(T);
}



void initialize_directions(void)
{
    double ii,iii;
    double dir_temp[3];
    double *dir_ptr = dir_temp;
    int i,j,k,l,m;
    int Nside, Nangles;
    
#ifdef USE_BRULS
//    AngularGrid(All.IdortBrulsNum);
//#if defined(BRULS_NUM) && defined(IDORT_NUM_ANGLES)
    int N = 2*All.IdortBrulsNum;
    Nangles = N*(N+2);
    if(Nangles == IDORT_NUM_ANGLES)
    {
        double musq[All.IdortBrulsNum];
        double delta = 2.0/(N-1.0);
        int sigma = All.IdortBrulsNum - 1;
        for(j = 0; j < All.IdortBrulsNum; j++)
            musq[j] = 1.0/3.0/(N-1.0) + delta*j;
        l = 0;
        for(i = 0; i <= sigma; i++)
        {
            int jj = sigma-i;
            for(j = 0; j <= jj; j++)
            {
                k = sigma-i-j;
                for(int m1 = 0; m1 < 2; m1++)
                {
                    double s1 = m1*2.0-1.0;
                    for(int m2 = 0; m2 < 2; m2++)
                    {
                        double s2 = m2*2.0-1.0;
                        for(int m3 = 0; m3 < 2; m3++)
                        {
                            double s3 = m3*2.0-1.0;
                            double q1 = s1*sqrt(musq[i]);
                            double q2 = s2*sqrt(musq[j]);
                            double q3 = s3*sqrt(musq[k]);
                            ray_dir[l][0] = q1;
                            ray_dir[l][1] = q2;
                            ray_dir[l][2] = q3;
                            ray_weight[l] = 1.0 / IDORT_NUM_ANGLES;
                            l++;
//                            mpi_printf("IDORT: Using BRULS. Ray vector = (%lf,%lf,%lf)\n", q1, q2, q3);
                        }
                    }
                }
            }
        }
        if(l != Nangles)
            mpi_printf("Warning: Error using BRULS. BRULS_NUM=%d, Nloop = %d not matching IDORT_NUM_ANGLES = %d\n", All.IdortBrulsNum, l, IDORT_NUM_ANGLES);
        mpi_printf("IDORT: Using BRULS. BRULS_NUM=%d, IDORT_NUM_ANGLES = %d\n", All.IdortBrulsNum, IDORT_NUM_ANGLES);
    }
    else
    {
        mpi_printf("Warning: Error using BRULS. BRULS_NUM=%d, Nangles = %d not matching IDORT_NUM_ANGLES = %d\n", All.IdortBrulsNum, Nangles, IDORT_NUM_ANGLES);
    }
#else
    
    if(IDORT_RAY_NUMDIMS == 2)
    {
        double dtheta, theta;
#ifdef IDORT_TEST
        if(IDORT_TEST == 2)
        {
            int N = IDORT_NUM_ANGLES;
            dtheta = 2.0 * M_PI / N;
            l = 0;
            for(j = 0; j < N; j++)
            {
                theta = dtheta * (j+0.5);
                double s;
                s = 1.0;
                ray_dir[l][0] = s*cos(theta);
                ray_dir[l][1] = s*sin(theta);
                ray_dir[l][2] = sqrt(1-s*s);
                ray_weight[l] = 1.0 / IDORT_NUM_ANGLES;
                l++;
            }
            if(l != IDORT_NUM_ANGLES)
                mpi_printf("Warning: Error using rays in 2D. Nloop = %d not matching IDORT_NUM_ANGLES = %d\n", l, IDORT_NUM_ANGLES);
        }
        else
#endif
        {
            int N = IDORT_NUM_ANGLES / 2;
            if(IDORT_NUM_ANGLES % 2 > 0)
                mpi_printf("Warning: Error using rays in 2D. Angle number not even, thus not consistent!\n");
            dtheta = 2.0 * M_PI / N;
            l = 0;
            for(j = 0; j < N; j++)
            {
                theta = dtheta * (j+0.5);
                for(i = 0; i < 2; i++)
                {
                    double s;
                    if(i == 0)
                        s = 1.0;
                    else
                        s = 1.0/sqrt(3.0);
                    ray_dir[l][0] = s*cos(theta);
                    ray_dir[l][1] = s*sin(theta);
                    ray_dir[l][2] = sqrt(1-s*s);
                    ray_weight[l] = 1.0 / IDORT_NUM_ANGLES;
                    l++;
                }
            }
            if(l != IDORT_NUM_ANGLES)
                mpi_printf("Warning: Error using rays in 2D. Nloop = %d not matching IDORT_NUM_ANGLES = %d\n", l, IDORT_NUM_ANGLES);
        }
    }
    else
    {
        if(IDORT_RAY_NUMDIMS == 1)
        {
            double dnx;
            dnx = 2.0 / (IDORT_NUM_ANGLES-1);
//            if(IDORT_NUM_ANGLES!=2)
//                mpi_printf("Warning: Error using rays in 1D. Angle number not consistent!\n");
            for(j = 0; j < IDORT_NUM_ANGLES; j++)
            {
                ray_dir[j][0] = 1.0 - dnx*j;
                ray_dir[j][1] = sqrt(1 - ray_dir[j][0] * ray_dir[j][0]);
                ray_dir[j][2] = 0.0;
                ray_weight[j] = 1.0 / IDORT_NUM_ANGLES;
            }
        }
        else
        {
#ifdef USE_HEALPIX
//#if defined(HEALPIX_NUM) && defined(IDORT_NUM_ANGLES)
            Nside = 1;
            if(All.IdortHealpixNum > 0)
            {
                for(i = 0; i < All.IdortHealpixNum; i++)
                {
                    Nside *= 2;
                }
            }
            Nangles = 12 * Nside * Nside;
            if(Nangles == IDORT_NUM_ANGLES)
            {
                for(i = 0; i < IDORT_NUM_ANGLES; i++)
                {
                    pix2vec_nest(Nside, i, dir_ptr);
                    for(j = 0; j < 3; j++)
                    {
                        ray_dir[i][j] = dir_temp[j];
                    }
                    ray_weight[i] = 1.0 / IDORT_NUM_ANGLES;
                }
                mpi_printf("IDORT: Using HEALPIX. HEALPIX_NUM=%d, IDORT_NUM_ANGLES = %d\n", All.IdortHealpixNum, IDORT_NUM_ANGLES);
            }
            else
            {
                mpi_printf("Warning: Error using HEALPIX. HEALPIX_NUM=%d, Nangles = %d not matching IDORT_NUM_ANGLES = %d\n", All.IdortHealpixNum, Nangles, IDORT_NUM_ANGLES);
            }
//#else
//            mpi_printf("Warning: Error using HEALPIX. HEALPIX_NUM or IDORT_NUM_ANGLES not defined!\n");
//#endif
#else
            
#ifdef ONEDIMS_SPHERICAL
            int nang = IDORT_NUM_ANGLES;
            int nzeta = (int)(IDORT_NUM_ANGLES / 2);
//            double coszeta_v[2 * nzeta], coszeta_f[2 * nzeta + 1], len_zeta[2 * nzeta];
            int zs = 0; // ze = 2*nzeta - 1;

                double dcoszeta = 1.0 / nzeta;
                coszeta_f[zs] = 1.0;
                for (int i = 1; i < nzeta; ++i) {
                    coszeta_f[i + zs] = 1.0 - i * dcoszeta;
                }
                coszeta_f[nzeta + zs] = 0.0;
                for (int i = nzeta + 1; i < 2 * nzeta + 1; ++i)
                    coszeta_f[i + zs] = -coszeta_f[2 * nzeta - i + zs];

                for (int i = 0; i < nzeta; ++i) {
                    coszeta_v[i + zs] = 0.5 * (coszeta_f[i + zs] + coszeta_f[i + zs + 1]);
                }
                // re-normalize
                double normalization = 2 * nzeta / sqrt(4 * nzeta * nzeta - 1);

                for (int i = 0; i < nzeta; ++i) {
                    coszeta_v[i + zs] *= normalization;
                }
                for (int i = nzeta; i < 2 * nzeta; ++i) {
                    coszeta_v[i + zs] = -coszeta_v[2 * nzeta - i - 1 + zs];
                }

                for (int i = 0; i < nzeta * 2; ++i) {
                    len_zeta[i] = coszeta_f[i] - coszeta_f[i + 1];
                }
            
            for (int n=0; n<2*nzeta; ++n) {
               ray_dir[n][0] = coszeta_v[n];
            }
            // equal weight for all cases
            for (int n=0; n<nang; ++n) {
              ray_weight[n] = 1.0/nang;
            }
#endif
#endif
        }
    }
#endif /* USE_BRULS */
}
        
        

void boundaries_IDORT(tessellation *T)
{
#ifdef IRT_STAR_CORE
    double x[3], r, pos_rgcore[3];

    for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
        int i = TimeBinsHydro.ActiveParticleList[idx];
        if(i < 0)
          continue;
        
          for(int l = 0; l< 3; l++)
          {
              x[l] = 0.0;
          }

#ifdef ID_RGCORE
            get_pos_rgcore(pos_rgcore);
#else
            pos_rgcore[0] = boxHalf_X;
            pos_rgcore[1] = boxHalf_Y;
            pos_rgcore[2] = boxHalf_Z;
#endif
            x[0] = P[i].Pos[0] - pos_rgcore[0];
            r = x[0]*x[0];
            if(NUMDIMS > 1)
            {
                x[1] = P[i].Pos[1] - pos_rgcore[1];
                r += x[1]*x[1];
                if(NUMDIMS > 2)
                {
                    x[2] = P[i].Pos[2] - pos_rgcore[2];
                    r += x[2]*x[2];
                }
            }
            r = sqrt(r);

          if(r <= All.IrtStarCoreRadius * All.IrtStarCoreRadiusFactor)
          {
//              mpi_printf("Boundaries");
              double r0 = 1e-5*All.IrtStarCoreRadius;
              double F_rad_ideal;
              if(r > r0)
              {
                  F_rad_ideal = All.IrtStarCoreLuminosity / (4.0 * M_PI * r * r);
              }
              else
              {
                  F_rad_ideal = All.IrtStarCoreLuminosity / (4.0 * M_PI * r0 * r0);
              }
              double b2 = 0.0, b[3];
              for(int q = 0; q < 3; q++)
              {
#ifdef IDORT_NEGLECT_VEL
                  b[q] = 0.0;
#else
                  b[q] = P[i].Vel[q] * All.UnitVelocity_in_cm_per_s / clight;
#endif
                  b2 += b[q] * b[q];
              }
              double gamma = 1.0 / sqrt(1.0 - b2);
              for(int j = 0; j < IDORT_NUM_ANGLES; j++)
              {
                  double Gamma = 1.0;
                  double nr = 0.0;
                  double dir_r;
                  for(int l = 0; l< 3; l++)
                  {
                      if(r > 0.0)
                      {
                          dir_r = x[l] / r;
                      }
                      else
                      {
                          dir_r = 0.0;
                      }
                      nr += dir_r * ray_dir[j][l];
                      Gamma -= b[l] * ray_dir[j][l];
                  }
                  Gamma *= gamma;
#ifdef IDORT_CORE_CONST_L
                  SphP[i].Intensity[j] = (RAD_CONST*pow(SphP[i].EOSTemperature,4)/4.0/M_PI + 3.0/(4.0*M_PI*CLIGHT_REAL) * F_rad_ideal * nr);
#else
                  SphP[i].Intensity[j] = (RAD_CONST*pow(SphP[i].EOSTemperature,4)/4.0/M_PI + 3.0/(4.0*M_PI*CLIGHT_REAL) * F_rad_ideal * nr) / pow(Gamma,4);
#endif
                  Cval[i].I_lastiter[j] = SphP[i].Intensity[j];
              }
          }
      }
#endif
    
    
    // we go through all the interfaces and look for boundary points
    int iface, j;
    face *VF  = T->VF;
    point *DP = T->DP;
    for(iface = 0; iface < T->Nvf; iface++)
      {
        int p1 = VF[iface].p1;
        int p2 = VF[iface].p2;

#ifndef ONEDIMS
        if(p1 < 0 || p2 < 0)
          continue;
#endif

        if(DP[p1].task == ThisTask && DP[p2].task == ThisTask)
          {
            int part1 = DP[p1].index;
            int part2 = DP[p2].index;

            if((part1 >= NumGas || part2 >= NumGas) && (part1 < NumGas || part2 < NumGas))
              {
                // we found a boundary point
                int pCell, pGhost;
                if(part1 >= NumGas)
                {
                  pCell = part2;
                  pGhost = part1 - NumGas;
                }
                if(part2 >= NumGas)
                {
                  pCell = part1;
                  pGhost = part2 - NumGas;
                }
                int p = DP[pCell].index;
                int pg = DP[pGhost].index;
                  /* normal vector pointing to "right" state */
                  double nx = DP[pGhost].x - DP[pCell].x;
                  double ny = DP[pGhost].y - DP[pCell].y;
                  double nz = DP[pGhost].z - DP[pCell].z;
                  double nn = sqrt(nx * nx + ny * ny + nz * nz);
                  nx /= nn;
                  ny /= nn;
                  nz /= nn;
                double x[3];
                x[0] = P[p].Pos[0]-boxHalf_X;
                x[1] = P[p].Pos[1]-boxHalf_Y;
                x[2] = P[p].Pos[2]-boxHalf_Z;
#ifdef ONEDIMS_SPHERICAL
                x[0] = P[p].Pos[0];
#endif
#ifdef IDORT_OUTFLOW
                  for(j = 0; j < IDORT_NUM_ANGLES; j++)
                  {
                      SphP[p].Intensity[j] = 0.0;
                      SphP[pg].Intensity[j] = 0.0;
                  }
#endif
#ifdef IDORT_TEST
                  if(IDORT_TEST == 2)
                  {
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if(x[1]<-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0])
                          {
                              if((fabs(ray_dir[j][0] - cos(45.0/180.0*M_PI))<1e-3 && fabs(ray_dir[j][1] - sin(45.0/180.0*M_PI))<1e-3 && fabs(x[0]-0.1)<1e-2) || (fabs(ray_dir[j][0] - cos(135.0/180.0*M_PI))<1e-3 && fabs(ray_dir[j][1] - sin(135.0/180.0*M_PI))<1e-3 && fabs(x[0]+0.1)<1e-2))
                              {
                                  SphP[pg].Intensity[j] = 0.8;
                              }
                              else
                              {
                                  SphP[pg].Intensity[j] = 0.0;
                              }
                          }
                          if(x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]>boxHalf_Y/boxHalf_X*x[0])
                          {
                              if(ray_dir[j][0]>0)
                                  SphP[pg].Intensity[j] = SphP[p].Intensity[j];
                              else
                                  SphP[pg].Intensity[j] = 0.0;
                          }
                      }
                  }
                  if(IDORT_TEST == 3)
                  {
                      if(x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0])
                      {
                          SphP[p].Utherm = 1.5;
                          SphP[p].Density = 1.0;
                          P[p].Mass = SphP[p].Density * SphP[p].Volume;
                          SphP[p].Energy = P[p].Mass * SphP[p].Utherm;
                          for(int l = 0; l< 3; l++)
                          {
                              SphP[p].Momentum[l] = 0.0;
                          }
                      }
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if(x[1]<-boxHalf_Y/boxHalf_X*x[0] && x[1]>boxHalf_Y/boxHalf_X*x[0])
                          {
                              double r = sqrt(ray_dir[j][0]*ray_dir[j][0] + ray_dir[j][1]*ray_dir[j][1]);
                              double cosine = ray_dir[j][0]/r;
                              double sine = ray_dir[j][1]/r;
                              if(((fabs(cosine - cos(15.0/180.0*M_PI))<1e-3 && fabs(sine - sin(15.0/180.0*M_PI))<1e-3) || (fabs(cosine - cos(-15.0/180.0*M_PI))<1e-3 && fabs(sine - sin(-15.0/180.0*M_PI))<1e-3)) && ray_dir[j][2]<1e-3)
                              {
                                  SphP[pg].Intensity[j] = rad_const * 1031.3 * IDORT_NUM_ANGLES / 24.0;
                              }
                              else
                              {
                                  SphP[pg].Intensity[j] = 0.0;
                              }
                          }
                          if(x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0])
                          {
                              if(ray_dir[j][0]>0)
                                  SphP[pg].Intensity[j] = SphP[p].Intensity[j];
                              else
                                  SphP[pg].Intensity[j] = 0.0;
                          }
                      }
                  }
                  if(IDORT_TEST == 5)
                  {
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if(x[0]<0)
                          {
//                              SphP[p].Intensity[j] = rad_const*1.0/4.0/M_PI;
                              SphP[pg].Intensity[j] = rad_const*1.0/4.0/M_PI;
                          }
                          if(x[0]>0)
                          {
                              if(ray_dir[j][0]>0)
                                  SphP[pg].Intensity[j] = SphP[p].Intensity[j];
                              else
                                  SphP[pg].Intensity[j] = 0.0;
                          }
                      }
                  }
                  if(IDORT_TEST == 6)
                  {
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if(P[p].Pos[0]<0.1)
                          {
//                              SphP[p].Intensity[j] = rad_const*1.0/4.0/M_PI;
                              SphP[pg].Intensity[j] = rad_const*1.0/4.0/M_PI;
                          }
                          if(P[p].Pos[0]>6.9)
                          {
                              if(ray_dir[j][0]>0)
                              {
                                  SphP[pg].Intensity[j] = SphP[p].Intensity[j];
                              }
                              else
                              {
                                  SphP[pg].Intensity[j] = 0.0;
                              }
                          }
                      }
                  }
                  if(IDORT_TEST == 8)
                  {
                          if(x[0]<0.0)
                          {
                              for(j = 0; j < IDORT_NUM_ANGLES; j++)
                              {
                                  SphP[pg].Intensity[j] = rad_const*1.0/4.0/M_PI;// / pow(SphP[p].Gamma_ni[j],4);
                              }
                              if(strcmp(All.InitCondFile, "ics_m1.2") == 0)
                              {
                                  P[pg].Vel[0] = 1.2;
                              }
                              if(strcmp(All.InitCondFile, "ics_m2") == 0)
                              {
                                  P[pg].Vel[0] = 2.0;
                              }
                              if(strcmp(All.InitCondFile, "ics_m3") == 0)
                              {
                                  P[pg].Vel[0] = 3.0;
                              }
                                  SphP[pg].Utherm = 9e-1;
                                  SphP[pg].Density = 1.0;
                                  P[pg].Mass = SphP[pg].Density * SphP[pg].Volume;
                                  SphP[pg].Momentum[0] = P[pg].Mass * P[pg].Vel[0];
                                  SphP[pg].Energy = 0.5 * P[pg].Mass * P[pg].Vel[0] * P[pg].Vel[0] + P[pg].Mass * SphP[pg].Utherm;
                          }
                          if(x[0]>0.0)
                          {
                              for(j = 0; j < IDORT_NUM_ANGLES; j++)
                              {
                                  SphP[pg].Intensity[j] = (SphP[p].Intensity[j]-SphP[p-1].Intensity[j])/(P[p].Pos[0]-P[p-1].Pos[0])*(P[pg].Pos[0]-P[p].Pos[0]) + SphP[p].Intensity[j];
                              }
                          }
                  }
                  if(IDORT_TEST == 9)
                  {
//                      double mu = 2.33;
//                      double T0 = 82.0;
//                      double Tu = 1.0 * T0;
//                      double rho0 = 7.02e-16;
//                      double rhou = 1e-10 * rho0;
//                      double g0 = 1.45e-6;
//                      double k0 = 0.0316 * pow(T0/10.0,2);
//                      double Fr_ref = 0.5 * clight * g0 / k0;
//                      double Er_ref = 9.26 * rad_const * pow(T0,4); //Fr_ref / clight - SphP[p].F_rad[1] / clight + SphP[p].E_rad; //9.26 * rad_const * pow(T0,4);
//                      double factor = 0.0;
//                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
//                      {
//                          if(ray_dir[j][1]>0)
//                              factor += ray_dir[j][1] * ray_weight[j];
//                      }
//                      factor *= 2.0;
//                      double Ip = Er_ref + Fr_ref / clight / factor;
//                      double Im = Er_ref - Fr_ref / clight / factor;
//                      Ip /= 4.0 * M_PI;
//                      Im /= 4.0 * M_PI;
                      double tau0 = 3.0;
                      double fE0 = 0.5;
                      double mu = 2.33;
                      double T0 = 82.0;
                      double F0 = clight*rad_const*pow(T0,4);
                      double cs20 = BOLTZMANN/PROTONMASS*T0/mu;
                      double k0 = 0.0316*pow(T0/10.0,2);
                      double sigma0 = tau0/k0;
                      double g0 = k0*F0/(fE0*clight);
                      double L0 = cs20/g0;
                      double rho0 = sigma0/L0;
                      double Er0 = rad_const * pow(T0,4);
//                      double Iy0 = Er0/(2.0*M_PI);
                      double Tu = 1e3 * T0;
                      double rhou = 1e-13 * rho0;
                      double fE = 0.0;
                      double fF = 0.0;
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if(ray_dir[j][1]>0)// && ray_dir[j][2]<1e-3)
                          {
                              fE += ray_weight[j];
                              fF += ray_dir[j][1] * ray_weight[j];
                          }
                      }
                      fE *= 2.0;
                      fF *= 2.0;
                      double Ip = 9.26*Er0 / fE + F0 / clight / fF;
                      double Im = 9.26*Er0 / fE - F0 / clight / fF;
                      Ip /= 4.0 * M_PI;
                      Im /= 4.0 * M_PI;
//                      double Ip = F0 / clight / fF;
                      
                      if(x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]>boxHalf_Y/boxHalf_X*x[0])
                      {
                          SphP[pg].Energy -= P[pg].Mass * SphP[pg].Utherm;
                          for(int l = 0; l< 3; l++)
                          {
                              SphP[pg].Energy -= SphP[pg].Momentum[l]*SphP[pg].Momentum[l]/2.0/P[pg].Mass;
                          }
                          SphP[pg].Utherm = BOLTZMANN / PROTONMASS / GAMMA_MINUS1 / mu * Tu;
                          SphP[pg].Density = rhou;
                          P[pg].Mass = SphP[pg].Density * SphP[pg].Volume;
                          SphP[pg].Energy += P[pg].Mass * SphP[pg].Utherm;
                          for(int l = 0; l< 3; l++)
                          {
                              SphP[pg].Momentum[l] = 0.0;
                          }
                      }
                        if(x[1]<-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0])
                        {
                            SphP[p].Energy -= P[p].Mass * SphP[p].Utherm;
                            for(int l = 0; l< 3; l++)
                            {
                                SphP[p].Energy -= SphP[p].Momentum[l]*SphP[p].Momentum[l]/2.0/P[p].Mass;
                            }
                            SphP[p].Utherm = BOLTZMANN / PROTONMASS / GAMMA_MINUS1 / mu * 10.0 * T0;
                            SphP[p].Density = rho0;
                            P[p].Mass = SphP[p].Density * SphP[p].Volume;
                            SphP[p].Energy += P[p].Mass * SphP[p].Utherm;
                            for(int l = 0; l< 3; l++)
                            {
                                SphP[p].Momentum[l] = 0.0;
                            }
                        }
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if(x[1]<-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0])
                          {
                              SphP[pg].Intensity[j] = Er0 / (4.0 * M_PI) * (1.0 + 3.0 * ray_dir[j][1]);
//                              if(fabs(ray_dir[j][1]-1.0)<1e-3 && ray_dir[j][2]<1e-3)
//                              {
//                                  SphP[pg].Intensity[j] = F0 / (4.0 * M_PI * clight) * IDORT_NUM_ANGLES;
//                              }
//                              else
//                              {
//                                  SphP[pg].Intensity[j] = 0.0;
//                              }
//                              if(ray_dir[j][2]<1e-3)
//                              {
//                                  if(ray_dir[j][1]>0)
//                                  {
//                                      SphP[pg].Intensity[j] = Ip;//Iy0;//Ip;
//                                  }
//                                  else
//                                  {
//                                      SphP[pg].Intensity[j] = Im;
//                                  }
//                              }
//                              else
//                                  SphP[pg].Intensity[j] = 0.0;
//                              SphP[p].Intensity[j] = SphP[pg].Intensity[j];
                          }
//                          if(x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]>boxHalf_Y/boxHalf_X*x[0])
//                          {
//                              if(ray_dir[j][1]>0)
//                                  SphP[pg].Intensity[j] = SphP[p].Intensity[j];
//                              else
//                                  SphP[pg].Intensity[j] = 0.0;
//////                              if(ray_dir[j][1]>0)
//////                              {
//////                                  SphP[pg].Intensity[j] = Ip;//Iy0;//Ip;
//////                              }
//////                              else
//////                              {
//////                                  SphP[pg].Intensity[j] = Im;//0.0;//Im;
//////                              }
//////                              SphP[p].Intensity[j] = SphP[pg].Intensity[j];
//                          }
                      }
                  }
                  if(IDORT_TEST == 10)
                  {
                      for(j = 0; j < IDORT_NUM_ANGLES; j++)
                      {
                          if((x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0] && ray_dir[j][0]>0) || (x[1]<-boxHalf_Y/boxHalf_X*x[0] && x[1]>boxHalf_Y/boxHalf_X*x[0] && ray_dir[j][0]<0) || (x[1]>-boxHalf_Y/boxHalf_X*x[0] && x[1]>boxHalf_Y/boxHalf_X*x[0] && ray_dir[j][1]>0) || (x[1]<-boxHalf_Y/boxHalf_X*x[0] && x[1]<boxHalf_Y/boxHalf_X*x[0] && ray_dir[j][1]<0))
                          {
                              SphP[pg].Intensity[j] = SphP[p].Intensity[j];
                          }
                          else
                          {
                              SphP[pg].Intensity[j] = 0.0;
                              
                          }
                      }
                  }
#endif
              }
          }
      }
}
#endif /* IDORT */
