//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_FLUXES_H
#define CHAMR_3D_CHARM_FLUXES_H

#include "charm_globals.h"

#ifdef CHARM_FLUX_RIM
    void rim_orig(  double* RI, double* EI, double* PI, double* UI, double* VI, double* WI,
                    double RB, double PB, double UB, double VB, double WB,
                    double RE, double PE, double UE, double VE, double WE, double gam);
#endif

void charm_calc_flux_godunov(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double* qc, double n[3]);
void charm_calc_flux_lf(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double* qc, double n[3]);
void charm_calc_flux_cd(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double* qc, double n[3]);
void charm_calc_flux_hllc(p4est_t *p4est, charm_prim_t prim[2], double* qu, double* qv, double* qw, double* qe, double qc[], double n[3]);


#endif //CHAMR_3D_CHARM_FLUXES_H
