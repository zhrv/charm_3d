//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_FLUXES_H
#define CHAMR_3D_CHARM_FLUXES_H

#include "charm_globals.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CHARM_FLUX_RIM
void rim_orig(  charm_real_t* RI, charm_real_t* EI, charm_real_t* PI, charm_real_t* UI, charm_real_t* VI, charm_real_t* WI,
                charm_real_t RB, charm_real_t PB, charm_real_t UB, charm_real_t VB, charm_real_t WB,
                charm_real_t RE, charm_real_t PE, charm_real_t UE, charm_real_t VE, charm_real_t WE, charm_real_t gam);
#endif

void charm_calc_flux_godunov(p4est_t *p4est, charm_prim_t prim[2], charm_real_t *qu, charm_real_t *qv, charm_real_t *qw, charm_real_t *qe,
                             charm_real_t *qc, charm_real_t n[3]);

void
charm_calc_flux_lf(p4est_t *p4est, charm_prim_t prim[2], charm_real_t *qu, charm_real_t *qv, charm_real_t *qw, charm_real_t *qe, charm_real_t *qc,
                   charm_real_t n[3]);

void
charm_calc_flux_cd(p4est_t *p4est, charm_prim_t prim[2], charm_real_t *qu, charm_real_t *qv, charm_real_t *qw, charm_real_t *qe, charm_real_t *qc,
                   charm_real_t n[3]);

void
charm_calc_flux_hllc(p4est_t *p4est, charm_prim_t prim[2], charm_real_t *qu, charm_real_t *qv, charm_real_t *qw, charm_real_t *qe, charm_real_t qc[],
                     charm_real_t n[3]);

#ifdef __cplusplus
}
#endif

#endif //CHAMR_3D_CHARM_FLUXES_H
