//
// Created by zhrv on 05.10.18.
//

#ifndef CHARM_DG_CHARM_EOS_H
#define CHARM_DG_CHARM_EOS_H

#include "charm_globals.h"


void charm_mat_eos_ideal            (p4est_t * p4est, charm_prim_t * p, int flag);
void charm_mat_eos_mix              (p4est_t * p4est, charm_prim_t * p, int flag);
void charm_mat_eos_ideal_low_mach   (p4est_t * p4est, charm_prim_t * p, int flag);
void charm_mat_eos_mix_low_mach     (p4est_t * p4est, charm_prim_t * p, int flag);
void charm_mat_eos_table            (p4est_t * p4est, charm_prim_t * p, int flag);


#endif //CHARM_DG_CHARM_EOS_H
