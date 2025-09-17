#include "charm_globals.h"

#ifndef CHARM_3D_COMPONENT_H
#define CHARM_3D_COMPONENT_H
#define M_H 0.00100797
#define M_O 0.0159994

#endif

typedef enum {
    H, O, H2, O2, OH, H2O, HO2
} Component;

typedef struct {
    Component component;
    charm_real_t M;
} ComponentInfo;

ComponentInfo getComponentInfo(Component component) {
    ComponentInfo info;
    switch (component) {
        case H:
            info.component = component;
            info.M = M_H;
            return info;
        case O:
            info.component = component;
            info.M = M_O;
            return info;
        case H2:
            info.component = component;
            info.M = 2*M_H;
            return info;
        case O2:
            info.component = component;
            info.M = 2*M_O;
            return info;
        case OH:
            info.component = component;
            info.M = M_O + M_H;
            return info;
        case H2O:
            info.component = component;
            info.M = 2*M_H + M_O;
            return info;
        case HO2:
            info.component = component;
            info.M = M_H + 2*M_O;
            return info;
    }
}

// nu_stage_table
const int nst[7][7] =
        {
                { -1, 1, 1, -1, 1, -1, 1}, //H
                { 1, -1, -1, 1, 0, 0, 0}, //O
                { 0, 0, -1, 1, -1, 1, -1}, //H2
                { -1, 1, 0, 0, 0, 0, -1}, //O2
                { 1, -1, 1, -1, -1, 1, 0}, //OH
                { 0, 0, 0, 0, 1, -1, 0}, //H2O
                { 0, 0, 0, 0, 0, 0, 1} //HO2
        };

const int nst_H[7] = {-1, 1, 1, -1, 1, -1, 1}; //H
const int nst_O[7] = {1, -1, -1, 1, 0, 0, 0}; //O
const int nst_H2[7] = {0, 0, -1, 1, -1, 1, -1}; //H2
const int nst_O2[7] = {-1, 1, 0, 0, 0, 0, -1}; //O2
const int nst_OH[7] = {1, -1, 1, -1, -1, 1, 0}; //OH
const int nst_H2O[7] = {0, 0, 0, 0, 1, -1, 0}; //H2O
const int nst_HO2[7] = {0, 0, 0, 0, 0, 0, 1}; //HO2