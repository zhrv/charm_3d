//
// Created by zhrv on 19.10.17.
//

#ifndef CHAMR_3D_CHARM_CONNECTIVITY_H
#define CHAMR_3D_CHARM_CONNECTIVITY_H

#define P4EST_ENABLE_DEBUG

#include <p4est_to_p8est.h>

#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>

#include "charm_globals.h"


p4est_connectivity_t* charm_conn_create(charm_ctx_t *ctx);


#endif //CHAMR_3D_CHARM_CONNECTIVITY_H
