//
// Created by zhrv on 05.10.18.
//

#ifndef CHARM_DG_CHARM_EOS_H
#define CHARM_DG_CHARM_EOS_H

#include "charm_globals.h"

#ifdef __cplusplus
extern "C" {
#endif


void charm_mat_eos_ideal(p4est_t *p4est, charm_prim_t *p, charm_eos_flag_t flag);

void charm_mat_eos_mix(p4est_t *p4est, charm_prim_t *p, charm_eos_flag_t flag);

void charm_mat_eos_table(p4est_t *p4est, charm_prim_t *p, charm_eos_flag_t flag);

charm_real_t charm_eos_get_r();

#ifdef __cplusplus
}
#endif

#endif //CHARM_DG_CHARM_EOS_H
