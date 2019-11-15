//
// Created by zhrv on 15.11.19.
//

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"
#include "charm_xml.h"
#include "charm_eos.h"


void charm_model_ns_timestep_diffusion(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    size_t c_count = charm_get_comp_count(p4est);
    if (c_count < 2) return;
}