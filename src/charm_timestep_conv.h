//
// Created by zhrv on 15.09.18.
//

#ifndef CHARM_DG_CHARM_TIMESTEP_CONV_H
#define CHARM_DG_CHARM_TIMESTEP_CONV_H
#include "charm_globals.h"

void charm_convect_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
void charm_convect_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data);

#endif //CHARM_DG_CHARM_TIMESTEP_CONV_H
