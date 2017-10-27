//
// Created by zhrv on 26.10.17.
//

#include "charm_fluxes.h"



void rim_orig(  double* RI, double* EI, double* PI, double* UI, double* VI, double* WI,
                double RB, double PB, double UB, double VB, double WB,
                double RE, double PE, double UE, double VE, double WE, double gam) {

    double AGAM = (gam - 1.0);
    double BGAM = (2.0 * sqrt(gam / AGAM));
    double CGAM = (1.0 / gam);
    double DGAM = (2.0 / AGAM);
    double EGAM = (AGAM / (gam + 1.0));
    double GGAM = (sqrt(gam * AGAM));
    double HGAM = (AGAM / 2.0);
    double FGAM = (3.0 * gam - 1.0);
    double OGAM = (AGAM / (2.0 * gam));
    double QGAM = (gam + 1.0);
    double PGAM = (QGAM / (2.0 * gam));
    double RGAM = (4.0 * gam);
    double SGAM = (gam * AGAM);
    double TGAM = (QGAM / 2.0);
    double UGAM = (sqrt(AGAM / gam));

    double RF, RS, EF, ES, SBL, SFL, SSL, SEL, D, FS1, F1, ZNB, PKB, ZFB, PKE, ZNE, F2, FS2, ZFE, DP, UBD, RUBD, UF, UED, RUED, US, PPE, PPB, P;
    double eps = RIM_EPS;
    double CB = sqrt(gam * PB / RB);
    double CE = sqrt(gam * PE / RE);
    double EB = CB * CB / SGAM;
    double EE = CE * CE / SGAM;
    double RCB = RB * CB;
    double RCE = RE * CE;
    double DU = UB - UE;
    if (DU < -2.0 * (CB + CE) / AGAM) {
        printf(" ATTENTION!!!  VACUUM \n");
        RF = 0.0;
        RS = 0.0;
        EF = 0.0;
        ES = 0.0;
        SBL = UB - CB;
        SFL = UB + 2.0 * CB / AGAM;
        SSL = UE - 2.0 * CE / AGAM;
        SEL = UE + CE;
    } else {
        P = (PB * RCE + PE * RCB + DU * RCB * RCE) / (RCB + RCE);
        do {

            if (P < eps) P = eps;

            PPB = P / PB;
            if (PB <= P) {
                PKB = PGAM * PPB + OGAM;
                ZNB = RCB * sqrt(PKB);
                F1 = (P - PB) / ZNB;
                FS1 = (QGAM * PPB + FGAM) / (RGAM * ZNB * PKB);
            }
            else {
                ZFB = CB * exp(log(PPB) * OGAM);
                F1 = DGAM * (ZFB - CB);
                FS1 = ZFB / (gam * P);
            }
            PPE = P / PE;
            if (PE <= P) {
                PKE = PGAM * PPE + OGAM;
                ZNE = RCE * sqrt(PKE);
                F2 = (P - PE) / ZNE;
                FS2 = (QGAM * PPE + FGAM) / (RGAM * ZNE * PKE);
            }
            else {
                ZFE = CE * exp(log(PPE) * OGAM);
                F2 = DGAM * (ZFE - CE);
                FS2 = ZFE / (gam * P);
            }
            DP = (DU - F1 - F2) / (FS1 + FS2);
            P = P + DP;
        } while (fabs(DU - F1 - F2) > eps);


        PPB = P / PB;
        PPE = P / PE;

//       ZFB=CB*PPB**OGAM;
//       ZFE=CE*PPE**OGAM;
        ZFB = CB * exp(log(PPB) * OGAM);
        ZFE = CE * exp(log(PPE) * OGAM);
        if (PB <= P) {
            D = UB - sqrt((TGAM * P + HGAM * PB) / RB);
            UBD = UB - D;
            RUBD = RB * UBD;
            RF = RUBD * RUBD / (PB - P + RUBD * UBD);
            UF = D + RUBD / RF;
            EF = P / (AGAM * RF);
            SBL = D;
            SFL = D;
        }
        else {
            EF = ZFB * ZFB / SGAM;
            UF = UB + DGAM * (CB - ZFB);
            RF = P / (AGAM * EF);
            SBL = UB - CB;
            SFL = UF - ZFB;
        }
        if (PE <= P) {
            D = UE + sqrt((TGAM * P + HGAM * PE) / RE);
            UED = UE - D;
            RUED = RE * UED;
            RS = RUED * RUED / (PE - P + RUED * UED);
            US = D + RUED / RS;
            ES = P / (AGAM * RS);
            SEL = D;
            SSL = D;
        }
        else {
            ES = ZFE * ZFE / SGAM;
            US = UE - DGAM * (CE - ZFE);
            RS = P / (AGAM * ES);
            SSL = US + ZFE;
            SEL = UE + CE;
        }
    }
//
// C     compute the interpolation value
    if (SEL<=0.0) {
        *RI= RE;
        *EI= EE;
        *UI= UE;
        *VI= VE;
        *WI= WE;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SBL>=0.0) {
        *RI= RB;
        *EI= EB;
        *UI= UB;
        *VI= VB;
        *WI= WB;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if ((SSL>=0.0)&&(SFL<=0.0)) {
        if (US>=0.0) {
            *RI= RF;
            *EI= EF;
            *UI= UF;
            *VI= VB;
            *WI= WB;
        } else {
            *RI= RS;
            *EI= ES;
            *UI= US;
            *VI= VE;
            *WI= WE;
        }
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SFL>0.0) {
        *UI= (UB+DGAM*GGAM*sqrt(EB))/(1+DGAM);
        *VI= VB;
        *WI= WB;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RB*exp(log((*EI)/EB)*(1/AGAM));
    } else {
        *UI= (UE-DGAM*GGAM*sqrt(EE))/(1+DGAM);
        *VI= VE;
        *WI= WE;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RE*exp(log((*EI)/EE)*(1/AGAM)) ;
    }

    *PI= AGAM*(*EI)*(*RI);

    return;
}


void calc_flux(double r_[2], double u_[2], double v_[2], double w_[2], double p_[2], double* qr, double* qu, double* qv, double* qw, double* qe, double n[3], int bnd)
{
#ifdef FLUX_RIM
    double ri, ei, pi, ui, vi, wi;
    switch (which_face) {
        case 0:                      /* -x side */
            if (bnd) {
                rim_orig(  &ri, &ei, &pi, &ui, &vi, &wi,
                           r_[1], p_[1], u_[1], v_[1], w_[1],
                           r_[0], p_[0], u_[0], v_[0], w_[0], GAM);
                *qr = ri*ui;
                *qu = (*qr)*ui+pi;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*ui;

            }
            else {
                rim_orig(  &ri, &ei, &pi, &ui, &vi, &wi,
                           r_[0], p_[0], u_[0], v_[0], w_[0],
                           r_[1], p_[1], u_[1], v_[1], w_[1], GAM);
                *qr = ri*ui;
                *qu = (*qr)*ui+pi;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*ui;
            }
            break;
        case 1:                      /* +x side */

            rim_orig(  &ri, &ei, &pi, &ui, &vi, &wi,
                       r_[0], p_[0], u_[0], v_[0], w_[0],
                       r_[1], p_[1], u_[1], v_[1], w_[1], GAM);
            *qr = ri*ui;
            *qu = (*qr)*ui+pi;
            *qv = (*qr)*vi;
            *qw = (*qr)*wi;
            *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*ui;
            break;
        case 2:                      /* -y side */
            if (bnd) {
                rim_orig(  &ri, &ei, &pi, &vi, &wi, &ui,
                           r_[1], p_[1], v_[1], w_[1], u_[1],
                           r_[0], p_[0], v_[0], w_[0], u_[0], GAM);
                *qr = ri*vi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi+pi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*vi;
            }
            else {
                rim_orig(  &ri, &ei, &pi, &vi, &wi, &ui,
                           r_[0], p_[0], v_[0], w_[0], u_[0],
                           r_[1], p_[1], v_[1], w_[1], u_[1], GAM);
                *qr = ri*vi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi+pi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*vi;
            }
            break;
        case 3:                      /* +y side */

            rim_orig(  &ri, &ei, &pi, &vi, &wi, &ui,
                       r_[0], p_[0], v_[0], w_[0], u_[0],
                       r_[1], p_[1], v_[1], w_[1], u_[1], GAM);
            *qr = ri*vi;
            *qu = (*qr)*ui;
            *qv = (*qr)*vi+pi;
            *qw = (*qr)*wi;
            *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*vi;
            break;
        case 4:                      /* -z side */
            if (bnd) {
                rim_orig(  &ri, &ei, &pi, &wi, &ui, &vi,
                           r_[1], p_[1], w_[1], u_[1], v_[1],
                           r_[0], p_[0], w_[0], u_[0], v_[0], GAM);
                *qr = ri*wi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi+pi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*wi;
            }
            else {
                rim_orig(  &ri, &ei, &pi, &wi, &ui, &vi,
                           r_[0], p_[0], w_[0], u_[0], v_[0],
                           r_[1], p_[1], w_[1], u_[1], v_[1], GAM);
                *qr = ri*wi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi+pi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*wi;
            }
            break;
        case 5:                      /* +z side */

            rim_orig(  &ri, &ei, &pi, &wi, &ui, &vi,
                       r_[0], p_[0], w_[0], u_[0], v_[0],
                       r_[1], p_[1], w_[1], u_[1], v_[1], GAM);
            *qr = ri*wi;
            *qu = (*qr)*ui;
            *qv = (*qr)*vi;
            *qw = (*qr)*wi+pi;
            *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*wi;
            break;
    }
#else
#ifdef FLUX_LF
    double fr[2], fu[2], fv[2], fw[2], fe[2];
    double ro[2], ru[2], rv[2], rw[2], re[2];
    double alpha;
    double v[2][3] = {{u_[0], v_[0], w_[0]}, {u_[1], v_[1], w_[1]}};
    double vn[2];
    vn[0] = scalar_prod(v[0], n);
    vn[1] = scalar_prod(v[1], n);
    int i;

    for (i = 0; i < 2; i++) {
        ro[i] = r_[i];
        ru[i] = r_[i] * u_[i];
        rv[i] = r_[i] * v_[i];
        rw[i] = r_[i] * w_[i];
        re[i] = p_[i] / (GAM - 1.0) + 0.5 * r_[i] * (u_[i] * u_[i] + v_[i] * v_[i] + w_[i] * w_[i]);

        fr[i] = r_[i]*vn[i];
        fu[i] = fr[i]*u_[i]+p_[i]*n[0];
        fv[i] = fr[i]*v_[i]+p_[i]*n[1];
        fw[i] = fr[i]*w_[i]+p_[i]*n[2];
        fe[i] = (re[i]+p_[i])*vn[i];
    }
    alpha = _MAX_(fabs(vn[0])+sqrt(GAM*p_[0]/r_[0]), fabs(vn[1])+sqrt(GAM*p_[1]/r_[1]));


    *qr = 0.5*(fr[0]+fr[1]-alpha*(ro[1]-ro[0]));
    *qu = 0.5*(fu[0]+fu[1]-alpha*(ru[1]-ru[0]));
    *qv = 0.5*(fv[0]+fv[1]-alpha*(rv[1]-rv[0]));
    *qw = 0.5*(fw[0]+fw[1]-alpha*(rw[1]-rw[0]));
    *qe = 0.5*(fe[0]+fe[1]-alpha*(re[1]-re[0]));
#endif // FLUX_LF
#endif // FLUX_RIM
}

